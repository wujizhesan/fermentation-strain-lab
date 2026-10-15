import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
VARIANTS = {'original': 'ale_marge_demo_v1', 'balanced': 'ale_marge_balancefix_v1'}
DESIGNS = {'F': '富马酸设计', 'Sz': '琥珀酸设计', 'Mz': '苹果酸设计'}
LAYERS = ['protein', 'rna']
MODES = ['released_ratio', 'neutral_ratio']


def sha(body):
    return hashlib.sha256(body).hexdigest()


def number(value):
    result = float(value)
    if not math.isfinite(result):
        raise ValueError('Nonfinite viewer value')
    return result


def verified_file(folder, report, name):
    path = (folder / name).resolve()
    if not path.is_relative_to(folder.resolve()):
        raise ValueError('Artifact path leaves the result directory')
    body = path.read_bytes()
    if sha(body) != report['artifact_sha256'].get(name):
        raise ValueError('Displayed artifact fingerprint changed: ' + name)
    return body


def table(body):
    return list(csv.DictReader(io.StringIO(body.decode('utf-8'))))


def load_viewer_data(results_root, snapshot_kind='verified_reference'):
    results_root = Path(results_root)
    expected_cases = {f'{design}_{layer}_{mode}' for design in DESIGNS for layer in LAYERS for mode in MODES}
    payload = {'schema_version': 1, 'snapshot_kind': snapshot_kind,
        'title': 'ALE多组学与代谢模型约束复核',
        'verified_backend': 'ReFramed 1.2.1 / OptLang / GLPK',
        'paper_cplex_behavior_verified': False, 'paper_exact_reproduction_completed': False,
        'strain_optimization_gain_validated': False, 'intracellular_flux_measured': False,
        'interval_semantics': '模型可行域区间；固定第一阶段偏差松弛上限，未固定第二阶段最小值。',
        'variants': {}, 'genes': [], 'summary': {'cases': 0, 'joint_lps': 0, 'range_endpoints': 0}}
    gene_rows = None
    for variant, name in VARIANTS.items():
        folder = results_root / name
        report_body = (folder / 'report.json').read_bytes()
        report = json.loads(report_body)
        verification = json.loads((folder / 'verification.json').read_text(encoding='utf-8'))
        if verification.get('passed') is not True or verification.get('report_sha256') != sha(report_body):
            raise ValueError('Current independently verified report required: ' + name)
        if report.get('paper_exact_reproduction_completed') is not False or report.get('intracellular_flux_measured') is not False:
            raise ValueError('Viewer evidence boundary changed')
        if report.get('proteome_balance_kept') is not (variant == 'balanced'):
            raise ValueError('Wrong result variant')
        cases = report['cases']
        ids = [case['case_id'] for case in cases]
        if len(ids) != len(set(ids)) or set(ids) != expected_cases:
            raise ValueError('Incomplete viewer case coverage')
        strict = json.loads(verified_file(folder, report, 'strict_reproduction_status.json'))
        if strict['strict_entry_allowed'] is not False:
            raise ValueError('Strict paper boundary changed')
        source = json.loads(verified_file(folder, report, 'source_audit.json'))
        genes = table(verified_file(folder, report, 'input_gene_ratios.csv'))
        if gene_rows is not None and genes != gene_rows:
            raise ValueError('Variant input genes differ')
        gene_rows = genes
        output_cases = []
        for case in cases:
            cid = case['case_id']
            ranges = table(verified_file(folder, report, cid + '/flux_ranges.csv'))
            pairs = [(row['condition'], row['reaction_id']) for row in ranges]
            if len(pairs) != 10 or len(set(pairs)) != 10 or {p[0] for p in pairs} != {'a', 'b'}:
                raise ValueError('Tracking range coverage changed')
            numeric_ranges = []
            for row in ranges:
                lower, upper = number(row['minimum']), number(row['maximum'])
                if lower > upper + 1e-7 or row['stage2_proteome_optimum_fixed'] != 'False':
                    raise ValueError('Invalid range semantics')
                numeric_ranges.append({'condition': row['condition'], 'reaction_id': row['reaction_id'],
                    'minimum': lower, 'maximum': upper, 'width': upper - lower})
            selected = {key: case[key] for key in ['case_id', 'design', 'layer', 'mode', 'input_genes',
                'stage1_error_optimum', 'stage2_proteome_optimum', 'reference_proteome_a', 'reference_proteome_b',
                'growth_fraction_a', 'growth_fraction_b', 'proteome_balance_kept', 'stage2_degenerate']}
            selected['ranges'] = numeric_ranges
            output_cases.append(selected)
        by_id = {case['case_id']: case for case in output_cases}
        counts = {'narrower': 0, 'same': 0, 'wider': 0}
        wider = []
        for design in DESIGNS:
            for layer in LAYERS:
                released = by_id[f'{design}_{layer}_released_ratio']['ranges']
                neutral = {(row['condition'], row['reaction_id']): row for row in by_id[f'{design}_{layer}_neutral_ratio']['ranges']}
                for row in released:
                    key = row['condition'], row['reaction_id']
                    if key not in neutral:
                        raise ValueError('Released and neutral tracking ranges differ')
                    baseline = neutral[key]['width']
                    delta = row['width'] - baseline
                    kind = 'narrower' if delta < -1e-7 else ('wider' if delta > 1e-7 else 'same')
                    counts[kind] += 1
                    if kind == 'wider':
                        wider.append({'design': design, 'layer': layer, 'condition': key[0], 'reaction_id': key[1],
                            'neutral_width': baseline, 'released_width': row['width']})
        payload['variants'][variant] = {'label': '保留平衡修订' if variant == 'balanced' else '原样实现',
            'report_sha256': sha(report_body), 'cases': output_cases, 'width_counts': counts,
            'wider_cases': wider, 'missing_paper_inputs': strict['missing_required_inputs'],
            'source_audit': source, 'independent_solver': verification['independent_solver']}
        payload['summary']['cases'] += len(cases)
        payload['summary']['joint_lps'] += verification['joint_lps_resolved']
        payload['summary']['range_endpoints'] += verification['range_endpoints_resolved']
    for row in gene_rows:
        if row['q_value_available'] != 'False' or row['paper_q_filter_applied'] != 'False':
            raise ValueError('Gene evidence boundary changed')
        payload['genes'].append({key: row[key] for key in ['design', 'layer', 'feature_id', 'gene_id']} | {
            'released_log2_ratio': number(row['released_log2_ratio']),
            'relative_expression_b_over_a': number(row['relative_expression_b_over_a']),
            'parental_n_positive': int(row['parental_n_positive']), 'evolved_n_positive': int(row['evolved_n_positive'])})
    for data in payload['variants'].values():
        for case in data['cases']:
            count = sum(row['design'] == case['design'] and row['layer'] == case['layer'] for row in payload['genes'])
            if count != case['input_genes']:
                raise ValueError('Viewer gene count differs from the report')
    payload['summary']['input_gene_ratios'] = len(payload['genes'])
    return payload


def render_viewer(payload):
    template = (ROOT / 'demo/marge_viewer.html').read_text(encoding='utf-8')
    encoded = json.dumps(payload, ensure_ascii=False, separators=(',', ':'), allow_nan=False)
    encoded = encoded.replace('<', '\\u003c').replace('>', '\\u003e').replace('&', '\\u0026')
    if template.count('__MARGE_PAYLOAD__') != 1:
        raise ValueError('Viewer payload marker must appear exactly once')
    return template.replace('__MARGE_PAYLOAD__', encoded)


def build_viewer(results_root, destination, snapshot_kind='verified_reference'):
    destination = Path(destination)
    if destination.exists():
        raise ValueError('Viewer already exists; choose a new output file')
    payload = load_viewer_data(results_root, snapshot_kind)
    body = render_viewer(payload).encode('utf-8')
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_bytes(body)
    return {'path': str(destination), 'bytes': len(body), 'sha256': sha(body), 'summary': payload['summary'],
        'snapshot_kind': snapshot_kind, 'display_source_report_sha256': {
            key: value['report_sha256'] for key, value in payload['variants'].items()}}


def main():
    parser = argparse.ArgumentParser(description='从已核验MARGE结果生成离线交互式报告。')
    parser.add_argument('--results-root', type=Path, default=ROOT / 'research/results')
    parser.add_argument('--output', type=Path, default=ROOT / 'research/results/ale_marge_viewer_v1/marge-demo.html')
    parser.add_argument('--snapshot-kind', choices=['verified_reference', 'run_result'], default='verified_reference')
    args = parser.parse_args()
    print(json.dumps(build_viewer(args.results_root, args.output, args.snapshot_kind), ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
