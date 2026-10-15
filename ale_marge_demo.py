import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import json
import math
from pathlib import Path

import cobra
from marge_sources import validate_marge_inputs
from marge_runtime import activate, INSTANCES, reframed, export_lp
from reframed.cobra.transcriptomics import marge


ROOT = Path(__file__).resolve().parent
PROTOCOL = ROOT / 'research/ale_marge_protocol.json'
DEFAULT_OUTPUT = ROOT / 'research/results/ale_marge_demo_v1'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def read_csv(path):
    with path.open(encoding='utf-8', newline='') as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    with path.open('w', encoding='utf-8', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def strict_status(protocol):
    return {'paper_exact_reproduction_completed': False, 'strict_entry_allowed': False,
        'missing_required_inputs': protocol['missing_paper_inputs'],
        'rna_q_filter_reproduced': False, 'protein_q_filter_reproduced': False,
        'extracellular_rate_constraints_reproduced': False,
        'version_selected_by_publication_date_not_confirmed_by_authors': True}


def prepared_model(protocol, design):
    path = ROOT / protocol['model_path']
    base = cobra.io.read_sbml_model(str(path))
    symbols = defaultdict(set)
    for line in (ROOT / protocol['gene_name_table']).read_text().splitlines():
        columns = line.split('\t')
        if len(columns) == 2 and columns[1]:
            symbols[columns[1]].add(columns[0])
    for symbol in protocol['knockouts'][design]:
        if len(symbols[symbol]) != 1:
            raise ValueError('Ambiguous genotype gene')
        base.genes.get_by_id(next(iter(symbols[symbol]))).knock_out()
    model = reframed.load_cbmodel(str(path), flavor='fbc2', use_infinity=False)
    if set(model.reactions) != {'R_' + r.id for r in base.reactions}:
        raise ValueError('Reaction identity differs between model loaders')
    for reaction in base.reactions:
        converted = model.reactions['R_' + reaction.id]
        converted.lb, converted.ub = reaction.bounds
    return model


def input_ratios(protocol, model, layer, design):
    rows = [r for r in read_csv(ROOT / protocol['abundance_path']) if r['omics'] == layer and r['design'] == design]
    counts = Counter(r['systematic_id'] for r in rows if r['systematic_id'])
    accepted, excluded = [], []
    for row in rows:
        gid = 'G_' + row['systematic_id']
        effect = row['released_log2_ratio_evolved_over_parental']
        reason = ''
        if not row['systematic_id'] or gid not in model.genes or gid in ['G_UNKNOWN', 'G_S0001']:
            reason = 'no_model_gene_or_pseudo_gene'
        elif counts[row['systematic_id']] != 1:
            reason = 'multiple_released_features_map_to_same_gene'
        elif not effect:
            reason = 'fewer_than_two_positive_observations_in_a_group'
        else:
            try:
                ratio = 2.0 ** float(effect)
            except OverflowError:
                ratio = math.inf
            if not math.isfinite(ratio) or ratio <= 0:
                reason = 'invalid_positive_ratio'
        if reason:
            excluded.append({'layer': layer, 'design': design, 'feature_id': row['feature_id'],
                'systematic_id': row['systematic_id'], 'reason': reason})
        else:
            accepted.append({'layer': layer, 'design': design, 'feature_id': row['feature_id'],
                'gene_id': gid, 'released_log2_ratio': float(effect), 'relative_expression_b_over_a': ratio,
                'parental_n_positive': int(row['parental_n_positive']), 'evolved_n_positive': int(row['evolved_n_positive']),
                'q_value_available': False, 'paper_q_filter_applied': False})
    return accepted, excluded


def proteome_row_audit(lp, protocol, keep_balance=False):
    """Upstream MARGE builds the proteome pseudo-metabolite balance as row 'proteome_<condition>'
    and afterwards adds a capacity row under the same name, which replaces the balance.
    With the diagnostic rename the second row is filed as 'proteome_cap_<condition>' and the balance survives.
    The audit records what the solved LP actually contains instead of assuming."""
    rows = {r['id']: r for r in lp['constraints']}
    variables = {v['id']: v for v in lp['variables']}
    audit = {}
    for condition in ['a', 'b']:
        usage = sum(value['value'] for name, value in variables.items()
            if name.startswith('u_') and name.endswith('_' + condition))
        synth = variables['proteome_synth_' + condition]['value']
        row = rows['proteome_' + condition]
        terms = dict(row['terms'])
        gene_terms = sorted(name for name in terms if name.startswith('u_'))
        ceiling = rows['proteome_cap_' + condition]['ub'] if keep_balance else None
        audit['sum_gene_usage_' + condition] = usage
        audit['proteome_synth_' + condition] = synth
        audit['proteome_row_terms_' + condition] = len(terms)
        if keep_balance:
            if row['lb'] != 0 or row['ub'] != 0 or not gene_terms:
                raise ValueError(f'Proteome balance row for condition {condition} missing after the rename')
            cap = rows['proteome_cap_' + condition]
            if set(dict(cap['terms'])) != {'proteome_synth_' + condition} or cap['lb'] is not None:
                raise ValueError('Renamed capacity row is malformed')
            if usage <= 1e-9:
                raise ValueError('Gene usage collapsed under the restored balance')
            if abs(synth - usage) > 1e-6 * max(1.0, usage):
                raise ValueError('Restored proteome mass balance is violated by the solved LP')
            audit['proteome_mass_balance_present_' + condition] = True
            audit['proteome_capacity_slack_' + condition] = ceiling - synth
        else:
            if set(terms) != {'proteome_synth_' + condition} or gene_terms:
                raise ValueError(f'Upstream proteome row for condition {condition} changed: '
                    f'{len(terms)} terms, {len(gene_terms)} gene-usage terms')
            if usage <= 1e-9 and synth <= 1e-9:
                raise ValueError('Gene usage is zero, proteome decoupling cannot be demonstrated')
            audit['proteome_mass_balance_present_' + condition] = False
    audit['proteome_capacity_binding'] = keep_balance and any(
        audit['proteome_capacity_slack_' + c] <= 1e-6 * max(1.0, audit['sum_gene_usage_' + c]) for c in ['a', 'b'])
    return audit


def run_case(protocol, design, layer, mode, accepted, folder, keep_balance=False):
    model = prepared_model(protocol, design)
    ratios = {r['gene_id']: r['relative_expression_b_over_a'] if mode == 'released_ratio' else 1.0 for r in accepted}
    fractions = {key: value / protocol['growth_rates_h_inverse']['WT']
        for key, value in protocol['growth_rates_h_inverse'][design].items()}
    constraints_a = {'R_EX_gly_e_': (-protocol['parental_glycine_uptake_cap'], 1000.0)}
    constraints_b = {'R_EX_gly_e_': (0.0, 1000.0)}
    INSTANCES.clear()
    flux_a, flux_b = marge(model, rel_expr=ratios, constraints_a=constraints_a, constraints_b=constraints_b,
        growth_frac_a=fractions['parental'], growth_frac_b=fractions['evolved'],
        step2_relax=protocol['step2_relax'])
    joint = next(x for x in INSTANCES if x.model is None)
    reference = next(x for x in INSTANCES if x.model is not None)
    if len(joint.joint_snapshots) != 2 or len(reference.calls) != 4:
        raise ValueError('Unexpected upstream optimization stages')
    for index, snapshot in enumerate(joint.joint_snapshots, 1):
        save_json(folder / f'stage{index}_lp.json', snapshot)
    transformed = reference.model
    flux_rows = []
    for condition, values in [('a', flux_a), ('b', flux_b)]:
        for rid in sorted(model.reactions):
            reaction = model.reactions[rid]
            bounds = (constraints_a if condition == 'a' else constraints_b).get(rid, (reaction.lb, reaction.ub))
            flux_rows.append({'condition': condition, 'reaction_id': rid, 'flux': float(values[rid]),
                'lower_bound': bounds[0], 'upper_bound': bounds[1]})
    write_csv(folder / 'net_fluxes.csv', flux_rows)
    ranges = []
    for rid in protocol['tracked_reactions'][design]:
        expression = transformed.convert_id_to_expr(rid)
        for condition in ['a', 'b']:
            objective = {key + '_' + condition: value for key, value in expression.items()}
            limits = [joint.solve(objective, minimize=minimize, get_values=False).fobj for minimize in [True, False]]
            ranges.append({'condition': condition, 'reaction_id': rid, 'minimum': float(limits[0]),
                'maximum': float(limits[1]), 'stage1_relaxation_fixed': protocol['step2_relax'],
                'stage2_proteome_optimum_fixed': False})
    write_csv(folder / 'flux_ranges.csv', ranges)
    pre = reference.calls
    stage1, stage2 = joint.joint_snapshots
    audit = proteome_row_audit(stage2, protocol, keep_balance)
    result = {'design': design, 'layer': layer, 'mode': mode, 'input_genes': len(ratios),
        'proteome_balance_kept': keep_balance, **audit,
        'proteome_mass_balance_absent': not keep_balance,
        'proteome_budget_constrains_fluxes': keep_balance,
        'stage2_degenerate': not keep_balance,
        'growth_fraction_a': fractions['parental'], 'growth_fraction_b': fractions['evolved'],
        'max_growth_a': pre[0]['objective_value'], 'max_growth_b': pre[2]['objective_value'],
        'reference_proteome_a': pre[1]['objective_value'], 'reference_proteome_b': pre[3]['objective_value'],
        'stage1_error_optimum': stage1['objective_value'], 'stage1_error_cap': stage1['objective_value'] * (1 + protocol['step2_relax']),
        'stage2_proteome_optimum': stage2['objective_value'], 'growth_a': flux_a[model.biomass_reaction],
        'growth_b': flux_b[model.biomass_reaction], 'statuses_all_optimal': True,
        'lp_variables': len(stage2['variables']), 'lp_constraints': len(stage2['constraints']),
        'net_reactions_per_condition': len(model.reactions)}
    if keep_balance:
        result['stage2_degeneracy_cause'] = None
        result['flux_basis'] = ('stage1 relaxation, restored proteome mass balance, capacity row kept as '
            'proteome_cap_<condition>, and stage2 minimises total proteome')
        result['range_semantics'] = 'stage1 relaxation cap, restored proteome balance, capacity row and reaction bounds'
    else:
        result['stage2_degeneracy_cause'] = ('upstream ReFramed adds the proteome pseudo-metabolite balance and the '
            'capacity row under the same name, so the balance is replaced and proteome_synth is decoupled '
            'from gene usage; the second-stage minimum therefore collapses to its lower bound of zero')
        result['flux_basis'] = 'stage1 relaxation only; no proteome budget and no second-stage proteome minimum'
        result['range_semantics'] = 'stage1 relaxation cap and reaction bounds only'
    save_json(folder / 'case.json', result)
    print(f'{design}/{layer}/{mode}: genes={len(ratios)}, error={stage1["objective_value"]:.6g}', flush=True)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument('--strict-paper', action='store_true')
    parser.add_argument('--keep-proteome-balance', action='store_true',
        help='diagnostic variant: file the capacity row as proteome_cap_<condition> so the proteome balance survives')
    args = parser.parse_args()
    protocol = json.loads(PROTOCOL.read_text(encoding='utf-8'))
    if args.strict_paper:
        raise RuntimeError('Strict paper reproduction unavailable: ' + '; '.join(protocol['missing_paper_inputs']))
    if args.keep_proteome_balance and args.output == DEFAULT_OUTPUT:
        raise ValueError('The balance-preserving variant must be written to a separate directory')
    input_validation = validate_marge_inputs(ROOT, protocol)
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError('Use a new empty output directory')
    software = activate(protocol, keep_proteome_balance=args.keep_proteome_balance)
    output.mkdir(parents=True, exist_ok=True)
    save_json(output / 'strict_reproduction_status.json', strict_status(protocol))
    cases, genes, exclusions = [], [], []
    for design in protocol['designs']:
        model = prepared_model(protocol, design)
        for layer in protocol['layers']:
            accepted, excluded = input_ratios(protocol, model, layer, design)
            if not accepted:
                raise ValueError('No mapped positive gene ratios')
            genes.extend(accepted); exclusions.extend(excluded)
            for mode in protocol['modes']:
                case_id = f'{design}_{layer}_{mode}'
                folder = output / case_id; folder.mkdir()
                result = run_case(protocol, design, layer, mode, accepted, folder, args.keep_proteome_balance)
                cases.append({'case_id': case_id, **result})
    write_csv(output / 'input_gene_ratios.csv', genes)
    write_csv(output / 'excluded_features.csv', exclusions)
    save_json(output / 'source_audit.json', protocol['source_audit'])
    report = {'protocol_sha256': sha(PROTOCOL), 'software': software, 'cases': cases,
        'actual_upstream_marge_executed': True, 'paper_exact_reproduction_completed': False,
        'constraints_use_descriptive_ratios_without_q_filter': True, 'rna_and_protein_fitted_separately': True,
        'intracellular_flux_measured': False, 'strain_optimization_gain_validated': False,
        'proteome_balance_kept': args.keep_proteome_balance,
        'proteome_budget_effective': args.keep_proteome_balance,
        'stage2_carries_information': args.keep_proteome_balance,
        'proteome_balance_override_confirmed_cases': 0 if args.keep_proteome_balance else len(cases),
        'proteome_capacity_binding_cases': sum(c['proteome_capacity_binding'] for c in cases),
        'limitations': protocol['limitations'], 'input_sha256': protocol['input_sha256'],
        'input_validation': input_validation,
        'verified_backend': 'ReFramed 1.2.1 / OptLang / GLPK',
        'paper_cplex_behavior_verified': False,
        'artifact_sha256': {p.relative_to(output).as_posix(): sha(p) for p in sorted(output.rglob('*')) if p.is_file()}}
    save_json(output / 'report.json', report)
    print(json.dumps({'output': str(output), 'cases': len(cases)}, ensure_ascii=False), flush=True)


if __name__ == '__main__':
    main()
