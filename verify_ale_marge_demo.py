import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path

import cobra
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
from marge_runtime import activate, reframed
from marge_sources import validate_marge_inputs
from reframed.core.transformation import gpr_transform


ROOT = Path(__file__).resolve().parent


def check(condition, message):
    if not condition:
        raise ValueError(message)


def close(actual, expected, tolerance=3e-6):
    check(math.isfinite(float(actual)) and abs(float(actual) - expected) <= tolerance * max(1, abs(expected)),
        f'Value mismatch: {actual} != {expected}')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path):
    with path.open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle))


def matrix_model(model, constraints):
    ids = list(model.reactions)
    positions = {rid: i for i, rid in enumerate(ids)}
    data, rows, columns = [], [], []
    for i, met in enumerate(model.metabolites):
        for rid, coefficient in model.metabolite_reaction_lookup()[met].items():
            rows.append(i); columns.append(positions[rid]); data.append(coefficient)
    matrix = coo_matrix((data, (rows, columns)), shape=(len(model.metabolites), len(ids))).tocsr()
    bounds = [constraints.get(rid, (model.reactions[rid].lb, model.reactions[rid].ub)) for rid in ids]
    bounds = [(None if a is None or not math.isfinite(a) else a, None if b is None or not math.isfinite(b) else b) for a, b in bounds]
    return ids, matrix, bounds


def reference_optima(model, constraints):
    ids, matrix, bounds = matrix_model(model, constraints)
    objective = np.zeros(len(ids)); objective[ids.index(model.biomass_reaction)] = -1
    result = linprog(objective, A_eq=matrix, b_eq=np.zeros(matrix.shape[0]), bounds=bounds, method='highs')
    check(result.success, 'HiGHS reference growth failed')
    growth = -result.fun
    bounds[ids.index(model.biomass_reaction)] = (growth, growth)
    objective[:] = 0; objective[ids.index('proteome_synth')] = 1
    result = linprog(objective, A_eq=matrix, b_eq=np.zeros(matrix.shape[0]), bounds=bounds, method='highs')
    check(result.success, 'HiGHS reference proteome failed')
    return growth, float(result.fun)


def lp_arrays(lp):
    check(lp['status'] == 'Optimal' and lp['direction'] == 'min', 'Nonoptimal or incorrect LP direction')
    ids = [v['id'] for v in lp['variables']]
    check(len(set(ids)) == len(ids), 'Duplicate LP variable')
    check(len({r['id'] for r in lp['constraints']}) == len(lp['constraints']), 'Duplicate LP constraint')
    positions = {rid: i for i, rid in enumerate(ids)}
    equality, eq_rhs, inequality, ub_rhs = [], [], [], []
    for row in lp['constraints']:
        check(len({x[0] for x in row['terms']}) == len(row['terms']), 'Duplicate LP coefficient')
        check(all(name in positions and math.isfinite(value) for name, value in row['terms']), 'Invalid LP term')
        lo, hi = row['lb'], row['ub']
        if lo is not None and hi is not None and lo == hi:
            equality.append(row['terms']); eq_rhs.append(lo)
        else:
            if hi is not None:
                inequality.append(row['terms']); ub_rhs.append(hi)
            if lo is not None:
                inequality.append([(name, -value) for name, value in row['terms']]); ub_rhs.append(-lo)
    def sparse(entries):
        rr, cc, vv = [], [], []
        for i, terms in enumerate(entries):
            for name, value in terms:
                rr.append(i); cc.append(positions[name]); vv.append(value)
        return coo_matrix((vv, (rr, cc)), shape=(len(entries), len(ids))).tocsr()
    arrays = {'A_eq': sparse(equality), 'b_eq': np.array(eq_rhs), 'A_ub': sparse(inequality),
        'b_ub': np.array(ub_rhs), 'bounds': [(v['lb'], v['ub']) for v in lp['variables']]}
    values = np.array([v['value'] for v in lp['variables']])
    check(np.isfinite(values).all(), 'Nonfinite LP solution')
    for value, (lo, hi) in zip(values, arrays['bounds']):
        check((lo is None or value >= lo - 3e-7) and (hi is None or value <= hi + 3e-7), 'LP variable bound violated')
    check(np.max(np.abs(arrays['A_eq'] @ values - arrays['b_eq']), initial=0) < 2e-6, 'LP equality balance failed')
    check(np.max(arrays['A_ub'] @ values - arrays['b_ub'], initial=0) < 2e-6, 'LP inequality violated')
    objective = np.array([v['objective'] for v in lp['variables']])
    close(float(objective @ values), lp['objective_value'])
    result = linprog(objective, **arrays, method='highs')
    check(result.success, 'HiGHS joint LP solve failed')
    close(lp['objective_value'], float(result.fun))
    return ids, arrays, values


def model_for_case(protocol, design):
    path = ROOT / protocol['model_path']
    cobra_model = cobra.io.read_sbml_model(str(path))
    genes = defaultdict(set)
    for line in (ROOT / protocol['gene_name_table']).read_text().splitlines():
        columns = line.split('\t')
        if len(columns) == 2 and columns[1]:
            genes[columns[1]].add(columns[0])
    for symbol in protocol['knockouts'][design]:
        check(len(genes[symbol]) == 1, 'Unresolved knockout identity')
        cobra_model.genes.get_by_id(next(iter(genes[symbol]))).knock_out()
    source = reframed.load_cbmodel(str(path), flavor='fbc2', use_infinity=False)
    check(set(source.reactions) == {'R_' + r.id for r in cobra_model.reactions}, 'Model loader identity mismatch')
    for reaction in cobra_model.reactions:
        source.reactions['R_' + reaction.id].lb, source.reactions['R_' + reaction.id].ub = reaction.bounds
    check(cobra_model.reactions.get_by_id('ATPM').bounds == (1, 1), 'Source maintenance changed')
    close(cobra_model.reactions.get_by_id('EX_glc_e_').lower_bound, -1.15)
    return source, gpr_transform(source, inplace=False, add_proteome=True)


def expected_inputs(protocol, source, layer, design):
    rows = [r for r in read_csv(ROOT / protocol['abundance_path']) if r['omics'] == layer and r['design'] == design]
    counts = Counter(r['systematic_id'] for r in rows if r['systematic_id'])
    valid, excluded = {}, []
    for row in rows:
        gid = 'G_' + row['systematic_id']
        reason = ''
        if not row['systematic_id'] or gid not in source.genes or gid in ['G_UNKNOWN', 'G_S0001']:
            reason = 'no_model_gene_or_pseudo_gene'
        elif counts[row['systematic_id']] != 1:
            reason = 'multiple_released_features_map_to_same_gene'
        elif not row['released_log2_ratio_evolved_over_parental']:
            reason = 'fewer_than_two_positive_observations_in_a_group'
        else:
            try:
                value = 2 ** float(row['released_log2_ratio_evolved_over_parental'])
            except OverflowError:
                value = math.inf
            if not math.isfinite(value) or value <= 0:
                reason = 'invalid_positive_ratio'
        if reason:
            excluded.append((layer, design, row['feature_id'], row['systematic_id'], reason))
        else:
            valid[gid] = row
    return valid, excluded


def validate_formulation(lp, transformed, ratios, constraints, fractions, growths, proteomes, stage, cap,
                         keep_balance=False):
    variables = {v['id']: v for v in lp['variables']}
    rows = {r['id']: r for r in lp['constraints']}
    expected_vars = {rid + '_' + c for rid in transformed.reactions for c in ['a', 'b']}
    expected_vars.update(gid + suffix for gid in ratios for suffix in ['_+', '_-'])
    check(set(variables) == expected_vars, 'LP variable coverage differs from MARGE formulation')
    expected_rows = {}
    for condition in ['a', 'b']:
        for rid, reaction in transformed.reactions.items():
            lo, hi = constraints[condition].get(rid, (reaction.lb, reaction.ub))
            if rid == transformed.biomass_reaction:
                lo, hi = growths[condition] * fractions[condition], None
            variable = variables[rid + '_' + condition]
            if lo is None or not math.isfinite(lo):
                check(variable['lb'] is None, 'Unexpected MARGE lower bound')
            else:
                close(variable['lb'], lo, tolerance=3e-6 if rid == transformed.biomass_reaction else 1e-12)
            check(variable['ub'] == (hi if hi is not None and math.isfinite(hi) else None), 'Unexpected MARGE upper bound')
        for met, reactions in transformed.metabolite_reaction_lookup().items():
            expected_rows[met + '_' + condition] = ({rid + '_' + condition: value for rid, value in reactions.items()}, 0, 0)
        # upstream ReFramed adds the capacity row under the same id as the proteome balance, replacing it;
        # the diagnostic variant files it as proteome_cap_<condition> and keeps the balance
        balance_id = 'proteome_' + condition
        replaced = expected_rows[balance_id]
        check(len(replaced[0]) > 1 and any(name.startswith('u_') for name in replaced[0]),
            'Proteome balance row carries no gene usage to override')
        capacity_id = ('proteome_cap_' + condition) if keep_balance else balance_id
        expected_rows[capacity_id] = ({'proteome_synth_' + condition: 1}, None, proteomes[condition])
    for gid, ratio in ratios.items():
        terms = {'u_' + gid[2:] + '_b': ratio ** -0.5, 'u_' + gid[2:] + '_a': -(ratio ** 0.5)}
        expected_rows[gid + '_c+'] = ({**terms, gid + '_-': 1}, 0, None)
        expected_rows[gid + '_c-'] = ({**terms, gid + '_+': -1}, None, 0)
        for suffix in ['_+', '_-']:
            check((variables[gid + suffix]['lb'], variables[gid + suffix]['ub']) == (0, None), 'Error slack bounds changed')
    if stage == 2:
        expected_rows['obj1'] = ({gid + suffix: 1 for gid in ratios for suffix in ['_+', '_-']}, None, cap)
    check(set(rows) == set(expected_rows), 'LP constraint coverage differs from MARGE formulation')
    for name, (terms, lo, hi) in expected_rows.items():
        row = rows[name]
        actual = dict(row['terms'])
        check(set(actual) == set(terms), 'LP stoichiometry or relative-expression terms changed')
        for key in terms:
            close(actual[key], terms[key], tolerance=2e-10)
        for field, expected in [('lb', lo), ('ub', hi)]:
            if expected is None:
                check(row[field] is None, 'Unexpected LP constraint bound')
            else:
                close(row[field], expected, tolerance=3e-6 if name.startswith('proteome') else 2e-10)
    for name, variable in variables.items():
        expected = (1 if name.endswith(('_+', '_-')) else 0) if stage == 1 else (1 if name in ['proteome_synth_a', 'proteome_synth_b'] else 0)
        close(variable['objective'], expected, tolerance=1e-12)


def verify(output, keep_balance=False):
    path = ROOT / 'research/ale_marge_protocol.json'
    protocol = json.loads(path.read_text(encoding='utf-8'))
    report = json.loads((output / 'report.json').read_text(encoding='utf-8'))
    check(report['protocol_sha256'] == sha(path) and report['input_sha256'] == protocol['input_sha256'], 'Protocol fingerprint mismatch')
    input_validation = validate_marge_inputs(ROOT, protocol)
    if 'input_validation' in report:
        check(report['input_validation'] == input_validation, 'Input validation receipt mismatch')
        check(report['verified_backend'] == 'ReFramed 1.2.1 / OptLang / GLPK'
            and report['paper_cplex_behavior_verified'] is False, 'Backend scope misrepresented')
    for name, expected in report['artifact_sha256'].items():
        check(sha(output / name) == expected, 'Artifact checksum mismatch')
    check(report['actual_upstream_marge_executed'] is True and report['paper_exact_reproduction_completed'] is False,
        'Demo misrepresented as exact-paper reproduction')
    check(report['constraints_use_descriptive_ratios_without_q_filter'] is True and report['rna_and_protein_fitted_separately'] is True,
        'Undocumented q filter or omics fusion')
    check(report['intracellular_flux_measured'] is False and report['strain_optimization_gain_validated'] is False, 'Unsupported experimental claim')
    check(report['limitations'] == protocol['limitations'], 'Omitted limitations')
    if keep_balance:
        check(report['proteome_balance_kept'] is True and report['proteome_budget_effective'] is True
            and report['stage2_carries_information'] is True
            and report['proteome_balance_override_confirmed_cases'] == 0,
            'Balance-preserving variant misrepresented')
    else:
        check(report['proteome_budget_effective'] is False and report['stage2_carries_information'] is False,
            'Proteome budget effectiveness misrepresented')
        check(report['proteome_balance_override_confirmed_cases'] == len(report['cases']),
            'Proteome balance override coverage misrepresented')
    strict = json.loads((output / 'strict_reproduction_status.json').read_text(encoding='utf-8'))
    check(strict['strict_entry_allowed'] is False and strict['paper_exact_reproduction_completed'] is False
        and strict['missing_required_inputs'] == protocol['missing_paper_inputs'], 'Strict-paper gate incorrectly opened')
    check(json.loads((output / 'source_audit.json').read_text(encoding='utf-8')) == protocol['source_audit'], 'Source audit changed')
    check(report['software'] == activate(protocol, keep_proteome_balance=keep_balance),
        'Software or compatibility receipt mismatch')
    input_rows = read_csv(output / 'input_gene_ratios.csv')
    grouped = defaultdict(dict)
    for row in input_rows:
        key = (row['design'], row['layer'])
        check(row['gene_id'] not in grouped[key], 'Duplicate input gene')
        grouped[key][row['gene_id']] = row
    expected_cases = {f'{d}_{layer}_{mode}' for d in protocol['designs'] for layer in protocol['layers'] for mode in protocol['modes']}
    check(len(report['cases']) == len(expected_cases) and {r['case_id'] for r in report['cases']} == expected_cases, 'Incomplete MARGE comparisons')
    exclusions, receipts = [], []
    for design in protocol['designs']:
        source, transformed = model_for_case(protocol, design)
        constraints = {'a': transformed.convert_constraints({'R_EX_gly_e_': (-1.0, 1000.0)}),
            'b': transformed.convert_constraints({'R_EX_gly_e_': (0.0, 1000.0)})}
        growths, proteomes = {}, {}
        for c in ['a', 'b']:
            growths[c], proteomes[c] = reference_optima(transformed, constraints[c])
        fractions = {'a': protocol['growth_rates_h_inverse'][design]['parental'] / protocol['growth_rates_h_inverse']['WT'],
            'b': protocol['growth_rates_h_inverse'][design]['evolved'] / protocol['growth_rates_h_inverse']['WT']}
        for layer in protocol['layers']:
            valid, excluded = expected_inputs(protocol, source, layer, design); exclusions.extend(excluded)
            check(set(grouped[(design, layer)]) == set(valid), 'Input gene selection changed')
            for gid, recorded in valid.items():
                row = grouped[(design, layer)][gid]
                check(row['feature_id'] == recorded['feature_id'], 'Wrong source feature for gene')
                close(row['released_log2_ratio'], float(recorded['released_log2_ratio_evolved_over_parental']))
                close(row['relative_expression_b_over_a'], 2 ** float(recorded['released_log2_ratio_evolved_over_parental']))
                check(row['q_value_available'] == 'False' and row['paper_q_filter_applied'] == 'False', 'Fabricated q-value filter')
            for mode in protocol['modes']:
                key = f'{design}_{layer}_{mode}'
                folder = output / key
                case = json.loads((folder / 'case.json').read_text(encoding='utf-8'))
                recorded = next(r for r in report['cases'] if r['case_id'] == key)
                check({k: v for k, v in recorded.items() if k != 'case_id'} == case, 'Case summary differs from report')
                check(case['design'] == design and case['layer'] == layer and case['mode'] == mode and case['input_genes'] == len(valid), 'Case identity mismatch')
                ratios = {gid: (2 ** float(row['released_log2_ratio_evolved_over_parental']) if mode == 'released_ratio' else 1.0) for gid, row in valid.items()}
                for c in ['a', 'b']:
                    close(case['max_growth_' + c], growths[c]); close(case['reference_proteome_' + c], proteomes[c])
                    close(case['growth_fraction_' + c], fractions[c])
                stage1 = json.loads((folder / 'stage1_lp.json').read_text(encoding='utf-8'))
                stage2 = json.loads((folder / 'stage2_lp.json').read_text(encoding='utf-8'))
                close(case['stage1_error_optimum'], stage1['objective_value'])
                cap = stage1['objective_value'] * (1 + protocol['step2_relax'])
                close(case['stage1_error_cap'], cap, tolerance=1e-12)
                close(case['stage2_proteome_optimum'], stage2['objective_value'])
                rows2 = {r['id']: r for r in stage2['constraints']}
                vars2 = {v['id']: v for v in stage2['variables']}
                for c in ['a', 'b']:
                    usage = sum(v['value'] for name, v in vars2.items()
                        if name.startswith('u_') and name.endswith('_' + c))
                    synth = vars2['proteome_synth_' + c]['value']
                    close(case['sum_gene_usage_' + c], usage, tolerance=1e-6)
                    close(case['proteome_synth_' + c], synth)
                    check(usage > 1e-6, 'Gene usage collapsed, proteome accounting not demonstrable')
                    if keep_balance:
                        row = rows2['proteome_' + c]
                        check(row['lb'] == 0 and row['ub'] == 0, 'Restored proteome row is not a balance')
                        check(any(name.startswith('u_') for name in dict(row['terms'])),
                            'Restored proteome balance has no gene usage')
                        cap_row = rows2['proteome_cap_' + c]
                        check(set(dict(cap_row['terms'])) == {'proteome_synth_' + c} and cap_row['lb'] is None,
                            'Renamed capacity row is malformed')
                        close(case['proteome_capacity_slack_' + c], cap_row['ub'] - synth, tolerance=1e-6)
                        close(synth, usage, tolerance=1e-6)
                    else:
                        row = rows2['proteome_' + c]
                        check(row['lb'] is None and row['ub'] is not None, 'Proteome row is not a capacity row')
                        check(set(dict(row['terms'])) == {'proteome_synth_' + c}, 'Proteome row terms changed')
                        check(synth <= 1e-9, 'Proteome synth is not at its lower bound')
                        check('proteome_cap_' + c not in rows2, 'Capacity row was renamed unexpectedly')
                if keep_balance:
                    check(case['proteome_mass_balance_absent'] is False and case['proteome_budget_constrains_fluxes'] is True
                        and case['stage2_degenerate'] is False and case['proteome_balance_kept'] is True,
                        'Balance-preserving variant not recorded')
                    check(case['stage2_proteome_optimum'] > 1e-6, 'Stage2 optimum did not become informative')
                else:
                    check(case['proteome_mass_balance_absent'] is True and case['proteome_budget_constrains_fluxes'] is False
                        and case['stage2_degenerate'] is True, 'Upstream proteome override not recorded')
                    override = protocol.get('proteome_balance_override', {})
                    check(override.get('expected') is True and case['proteome_row_terms_a'] == override.get('expected_terms_per_row')
                        and case['proteome_row_terms_b'] == override.get('expected_terms_per_row'),
                        'Protocol does not declare the proteome balance override')
                for index, lp in enumerate([stage1, stage2], 1):
                    validate_formulation(lp, transformed, ratios, constraints, fractions, growths, proteomes, index, cap,
                        keep_balance)
                    ids, arrays, values = lp_arrays(lp)
                check(case['lp_variables'] == len(ids) and case['lp_constraints'] == len(stage2['constraints']) and case['statuses_all_optimal'] is True, 'Case LP size or status mismatch')
                point = dict(zip(ids, values))
                net_rows = read_csv(folder / 'net_fluxes.csv')
                check(len(net_rows) == 2 * len(source.reactions) == 2 * case['net_reactions_per_condition']
                    and {(r['condition'], r['reaction_id']) for r in net_rows} == {(c, rid) for c in ['a', 'b'] for rid in source.reactions}, 'Incomplete net flux vectors')
                for row in net_rows:
                    expr = transformed.convert_id_to_expr(row['reaction_id'])
                    expected = sum(value * point[rid + '_' + row['condition']] for rid, value in expr.items())
                    close(row['flux'], expected)
                    reaction = source.reactions[row['reaction_id']]
                    bounds = (reaction.lb, reaction.ub)
                    if row['reaction_id'] == 'R_EX_gly_e_':
                        bounds = (-1.0 if row['condition'] == 'a' else 0.0, 1000.0)
                    close(row['lower_bound'], bounds[0]); close(row['upper_bound'], bounds[1])
                    check(bounds[0] - 3e-7 <= float(row['flux']) <= bounds[1] + 3e-7, 'Original net flux bound violated')
                for c in ['a', 'b']:
                    close(case['growth_' + c], point[transformed.biomass_reaction + '_' + c])
                intervals = read_csv(folder / 'flux_ranges.csv')
                wanted = {(c, rid) for c in ['a', 'b'] for rid in protocol['tracked_reactions'][design]}
                check(len(intervals) == len(wanted) and {(r['condition'], r['reaction_id']) for r in intervals} == wanted, 'Incomplete MARGE flux ranges')
                for row in intervals:
                    expression = transformed.convert_id_to_expr(row['reaction_id'])
                    objective = np.zeros(len(ids)); positions = {rid: i for i, rid in enumerate(ids)}
                    for rid, value in expression.items():
                        objective[positions[rid + '_' + row['condition']]] = value
                    for field, direction in [('minimum', 1), ('maximum', -1)]:
                        result = linprog(direction * objective, **arrays, method='highs')
                        check(result.success, 'HiGHS MARGE range solve failed')
                        close(row[field], direction * float(result.fun), tolerance=1e-5)
                    close(row['stage1_relaxation_fixed'], 0.1)
                    check(row['stage2_proteome_optimum_fixed'] == 'False', 'Range semantics changed')
                receipts.append({'case_id': key, 'stage1_optimum_highs': stage1['objective_value'],
                    'stage2_optimum_highs': stage2['objective_value'], 'range_endpoints_resolved': 2 * len(intervals),
                    'sum_gene_usage_a': case['sum_gene_usage_a'], 'sum_gene_usage_b': case['sum_gene_usage_b'],
                    'proteome_synth_a': case['proteome_synth_a'], 'proteome_synth_b': case['proteome_synth_b'],
                    'proteome_capacity_binding': case['proteome_capacity_binding'],
                    'proteome_balance_override_confirmed': not keep_balance})
                print('verified ' + key, flush=True)
    recorded_exclusions = read_csv(output / 'excluded_features.csv')
    check(len(recorded_exclusions) == len(exclusions) and Counter(tuple(r[x] for x in ['layer', 'design', 'feature_id', 'systematic_id', 'reason']) for r in recorded_exclusions) == Counter(exclusions), 'Feature exclusion audit changed')
    check(set(grouped) == {(d, layer) for d in protocol['designs'] for layer in protocol['layers']}, 'Unexpected input layer')
    return {'passed': True, 'independent_solver': 'SciPy HiGHS', 'cases': receipts, 'joint_lps_resolved': 2 * len(receipts),
        'range_endpoints_resolved': sum(r['range_endpoints_resolved'] for r in receipts),
        'reference_optima_resolved': 4 * len(protocol['designs']), 'input_gene_rows_recomputed': len(input_rows),
        'proteome_balance_override_confirmed_cases': sum(r['proteome_balance_override_confirmed'] for r in receipts),
        'proteome_balance_kept': keep_balance, 'proteome_budget_effective': keep_balance,
        'stage2_carries_information': keep_balance,
        'proteome_capacity_binding_cases': sum(r['proteome_capacity_binding'] for r in receipts),
        'report_sha256': sha(output / 'report.json'), 'exact_paper_reproduction': False}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT / 'research/results/ale_marge_demo_v1')
    parser.add_argument('--keep-proteome-balance', action='store_true')
    args = parser.parse_args()
    result = verify(args.output.resolve(), keep_balance=args.keep_proteome_balance)
    (args.output / 'verification.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: v for k, v in result.items() if k != 'cases'}), flush=True)


if __name__ == '__main__':
    main()
