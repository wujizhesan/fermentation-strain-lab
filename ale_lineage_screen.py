import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics

from modeling import solve_linear


ROOT = Path(__file__).resolve().parent
PROTOCOL = ROOT / 'research/ale_lineage_screen_protocol.json'
DEFAULT_OUTPUT = ROOT / 'research/results/ale_lineage_screen_v1'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def read_csv(path):
    with Path(path).open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream))


def write_csv(path, rows):
    with Path(path).open('w', encoding='utf-8', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def build_features(observations, endpoints, protocol):
    selected = [row for row in endpoints if row['design'] == protocol['design']]
    expected = {f'Chassis_F-{group}-{clone}' for group in [1, 2, 3] for clone in ['S', 'M', 'L']}
    ids = [row['isolate_id'] for row in selected]
    if len(ids) != len(set(ids)) or set(ids) != expected:
        raise ValueError('Expected nine distinct fumarate isolates')
    result = []
    for endpoint in sorted(selected, key=lambda row: row['isolate_id']):
        identity = endpoint['isolate_id']
        group = identity.rsplit('-', 1)[0]
        if endpoint['lineage_group'] != group or endpoint['product'] != protocol['product']:
            raise ValueError('Source lineage or product identity changed')
        points = [row for row in observations if row['isolate_id'] == identity and float(row['recorded_time_h']) <= protocol['cutoff_h']]
        points.sort(key=lambda row: float(row['recorded_time_h']))
        times = [float(row['recorded_time_h']) for row in points]
        values = [float(row['od600']) for row in points]
        if not times or times[0] != 0 or len(times) != len(set(times)):
            raise ValueError('Initial OD required; duplicate early time rejected')
        if any(not math.isfinite(value) or value <= 0 for value in values):
            raise ValueError('Finite positive OD required')
        titer = float(endpoint['titer_g_L'])
        if not math.isfinite(titer) or titer < 0:
            raise ValueError('Finite nonnegative endpoint titer required')
        initial, latest = values[0], values[-1]
        result.append({'isolate_id': identity, 'lineage_group': group, 'initial_od600': initial,
            'early_od600': latest, 'early_od_fold': latest / initial,
            'log_early_od_fold': math.log(latest / initial), 'early_recorded_time_h': times[-1],
            'early_points': len(points), 'early_source_excel_rows': json.dumps([int(row['source_excel_row']) for row in points]),
            'titer_g_L': titer, 'endpoint_source_sheet': endpoint['source_sheet'],
            'endpoint_source_excel_row': int(endpoint['source_excel_row']),
            'assay_collection_time_h': endpoint['assay_collection_time_h'],
            'assay_return_time_verified': endpoint['assay_return_time_verified']})
    return result


def fit_ridge(train, features, penalty):
    if not train or not math.isfinite(penalty) or penalty <= 0:
        raise ValueError('Positive ridge penalty and training data required')
    means = [statistics.fmean(row[key] for row in train) for key in features]
    scales = [math.sqrt(statistics.fmean((row[key] - mean) ** 2 for row in train)) for key, mean in zip(features, means)]
    scales = [value if value > 1e-12 else 1.0 for value in scales]
    vectors = [[1.0, *[(row[key] - mean) / scale for key, mean, scale in zip(features, means, scales)]] for row in train]
    size = len(features) + 1
    gram = [[sum(vector[i] * vector[j] for vector in vectors) for j in range(size)] for i in range(size)]
    cross = [sum(vector[i] * row['titer_g_L'] for vector, row in zip(vectors, train)) for i in range(size)]
    for index in range(1, size):
        gram[index][index] += penalty
    return {'features': features, 'means': means, 'scales': scales, 'coefficients': solve_linear(gram, cross),
        'ridge_penalty': penalty, 'training_isolate_ids': [row['isolate_id'] for row in train],
        'training_lineage_groups': sorted({row['lineage_group'] for row in train})}


def predict(model, row):
    vector = [1.0, *[(row[key] - mean) / scale for key, mean, scale in zip(model['features'], model['means'], model['scales'])]]
    value = sum(coefficient * item for coefficient, item in zip(model['coefficients'], vector))
    if not math.isfinite(value):
        raise ValueError('Nonfinite prediction')
    return value


def selection_weights(scores, tolerance):
    if not scores or any(not math.isfinite(value) for value in scores):
        raise ValueError('Finite candidate scores required')
    maximum = max(scores)
    chosen = [index for index, value in enumerate(scores) if abs(value - maximum) <= tolerance]
    return [1.0 / len(chosen) if index in chosen else 0.0 for index in range(len(scores))]


def fold_decisions(records, heldout_group, protocol):
    train = [row for row in records if row['lineage_group'] != heldout_group]
    test = [row for row in records if row['lineage_group'] == heldout_group]
    if len(train) != 6 or len(test) != 3 or len({row['lineage_group'] for row in train}) != 2:
        raise ValueError('Each fold must train on two complete populations and hold out three isolates')
    features = protocol['features']
    model = fit_ridge(train, features, protocol['ridge_penalty'])
    scores = {'uniform_random_expectation': [0.0] * len(test),
        'early_density': [row['early_od600'] for row in test],
        'early_density_fold': [row['early_od_fold'] for row in test],
        'fixed_ridge': [predict(model, row) for row in test]}
    decisions = []
    for method in protocol['methods']:
        weights = selection_weights(scores[method], protocol['score_tie_absolute_tolerance'])
        for row, score, weight in zip(test, scores[method], weights):
            decisions.append({'heldout_lineage': heldout_group, 'method': method, 'isolate_id': row['isolate_id'],
                'lineage_group': row['lineage_group'], 'candidate_score': score, 'selection_probability': weight,
                'score_is_predicted_titer': method == 'fixed_ridge', 'target_assay_budget': 1,
                'common_historical_titer_labels': len(train), 'test_titer_used_to_score': False})
    return model, decisions


def evaluate_decisions(decisions, records, protocol):
    labels = {row['isolate_id']: row['titer_g_L'] for row in records}
    results = []
    for group in sorted({row['heldout_lineage'] for row in decisions}):
        for method in protocol['methods']:
            rows = [row for row in decisions if row['heldout_lineage'] == group and row['method'] == method]
            if not math.isclose(sum(row['selection_probability'] for row in rows), 1.0, abs_tol=1e-12):
                raise ValueError('Selection probabilities must sum to one')
            best = max(labels[row['isolate_id']] for row in rows)
            expected = sum(row['selection_probability'] * labels[row['isolate_id']] for row in rows)
            hit = sum(row['selection_probability'] for row in rows if math.isclose(labels[row['isolate_id']], best, rel_tol=0, abs_tol=1e-12))
            results.append({'heldout_lineage': group, 'method': method, 'candidate_count': len(rows),
                'target_assay_budget': 1, 'selected_isolates': json.dumps([row['isolate_id'] for row in rows if row['selection_probability'] > 0]),
                'expected_selected_titer_g_L': expected, 'heldout_best_titer_g_L': best,
                'fraction_of_heldout_best': expected / best if best > 0 else None,
                'absolute_regret_g_L': best - expected, 'top1_hit_probability': hit,
                'new_experimental_assays_performed': 0, 'prospective_assay_availability_verified': False})
    return results


def summarize(metrics, decisions, records, protocol):
    labels = {row['isolate_id']: row['titer_g_L'] for row in records}
    summary = {}
    for method in protocol['methods']:
        rows = [row for row in metrics if row['method'] == method]
        summary[method] = {key: statistics.fmean(row[key] for row in rows) for key in [
            'expected_selected_titer_g_L', 'fraction_of_heldout_best', 'absolute_regret_g_L', 'top1_hit_probability']}
        summary[method]['independent_population_groups'] = len(rows)
    ridge = [row for row in decisions if row['method'] == 'fixed_ridge']
    summary['fixed_ridge']['prediction_mae_g_L'] = statistics.fmean(abs(row['candidate_score'] - labels[row['isolate_id']]) for row in ridge)
    summary['fixed_ridge']['negative_predictions'] = sum(row['candidate_score'] < 0 for row in ridge)
    return summary


def validate_protocol_inputs(root, protocol):
    for name, expected in protocol['input_sha256'].items():
        if sha(Path(root) / name) != expected:
            raise ValueError('Frozen input changed: ' + name)
    for name, expected in protocol['implementation_sha256'].items():
        if sha(Path(root) / name) != expected:
            raise ValueError('Frozen implementation changed: ' + name)
    base = Path(root) / 'research/results/ale_chassis_audit'
    verified = json.loads((base / 'verification.json').read_text(encoding='utf-8'))
    if verified['verified'] is not True or verified['report_sha256'] != sha(base / 'report.json'):
        raise ValueError('Verified upstream screen required')


def run(output=DEFAULT_OUTPUT):
    protocol = json.loads(PROTOCOL.read_text(encoding='utf-8'))
    validate_protocol_inputs(ROOT, protocol)
    output = Path(output).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError('Choose a new empty output directory')
    source = ROOT / 'research/results/ale_chassis_audit'
    records = build_features(read_csv(source / 'isolate_observations.csv'), read_csv(source / 'isolate_endpoints.csv'), protocol)
    models, decisions = {}, []
    for group in sorted({row['lineage_group'] for row in records}):
        model, rows = fold_decisions(records, group, protocol)
        models[group] = model
        decisions.extend(rows)
    metrics = evaluate_decisions(decisions, records, protocol)
    summary = summarize(metrics, decisions, records, protocol)
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output / 'features.csv', records)
    write_csv(output / 'candidate_decisions.csv', decisions)
    write_csv(output / 'lineage_metrics.csv', metrics)
    save_json(output / 'fitted_models.json', models)
    save_json(output / 'summary.json', summary)
    report = {'protocol_sha256': sha(PROTOCOL), 'input_sha256': protocol['input_sha256'],
        'species': protocol['species'], 'design': 'F', 'product': 'fumarate',
        'endpoint': 'released titer_g_L', 'isolates': len(records), 'population_groups': len(models),
        'early_source_points_used': sum(row['early_points'] for row in records),
        'actual_ridge_fits': len(models), 'fold_method_comparisons': len(metrics),
        'test_candidate_score_rows': len(decisions), 'summary': summary,
        'dataset_seen_before_protocol_lock': True, 'new_blind_validation': False,
        'prospective_selection_validated': False, 'productivity_scored': False,
        'omics_used_as_isolate_features': False, 'strain_optimization_gain_validated': False,
        'selection_expectations_are_retrospective': True, 'limitations': protocol['limitations'],
        'comparison_scope': 'Within each heldout population of3; not comparable to prior global9-isolate selection ratios.',
        'artifact_sha256': {path.name: sha(path) for path in sorted(output.iterdir()) if path.is_file()}}
    save_json(output / 'report.json', report)
    return report


def main():
    parser = argparse.ArgumentParser(description='公开富马酸ALE数据的整谱系留出选株回放。')
    parser.add_argument('--output', type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument('--strict-prospective', action='store_true')
    args = parser.parse_args()
    if args.strict_prospective:
        raise ValueError('Prospective scoring unavailable: target collection/return clocks and independent new data are missing')
    report = run(args.output)
    print(json.dumps({'output': str(args.output), 'summary': report['summary'], 'new_blind_validation': False}, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
