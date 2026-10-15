import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(path):
    with Path(path).open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream))


def close(actual, expected, label):
    if not math.isclose(float(actual), float(expected), rel_tol=1e-9, abs_tol=1e-10):
        raise ValueError('Independent numerical check failed: ' + label)


def verify(output):
    output = Path(output)
    protocol_path = ROOT / 'research/ale_lineage_screen_protocol.json'
    protocol = json.loads(protocol_path.read_text(encoding='utf-8'))
    report = json.loads((output / 'report.json').read_text(encoding='utf-8'))
    if report['protocol_sha256'] != digest(protocol_path):
        raise ValueError('Protocol changed')
    for name, expected in protocol['input_sha256'].items():
        if digest(ROOT / name) != expected:
            raise ValueError('Frozen input changed')
    for name, expected in protocol['implementation_sha256'].items():
        if digest(ROOT / name) != expected:
            raise ValueError('Frozen implementation changed')
    for name, expected in report['artifact_sha256'].items():
        if digest(output / name) != expected:
            raise ValueError('Output fingerprint changed')
    if any(report[key] is not False for key in ['new_blind_validation', 'prospective_selection_validated',
        'productivity_scored', 'omics_used_as_isolate_features', 'strain_optimization_gain_validated']):
        raise ValueError('Scientific scope changed')
    source = ROOT / 'research/results/ale_chassis_audit'
    endpoints = [row for row in rows(source / 'isolate_endpoints.csv') if row['design'] == 'F']
    observations = rows(source / 'isolate_observations.csv')
    labels = {row['isolate_id']: float(row['titer_g_L']) for row in endpoints}
    feature_rows = {}
    for endpoint in endpoints:
        identity = endpoint['isolate_id']
        points = sorted([row for row in observations if row['isolate_id'] == identity
            and float(row['recorded_time_h']) <= protocol['cutoff_h']], key=lambda row: float(row['recorded_time_h']))
        initial, latest = float(points[0]['od600']), float(points[-1]['od600'])
        feature_rows[identity] = {'initial_od600': initial, 'early_od600': latest,
            'early_od_fold': latest / initial, 'log_early_od_fold': float(np.log(latest / initial)),
            'early_recorded_time_h': float(points[-1]['recorded_time_h']), 'early_points': len(points),
            'lineage_group': endpoint['lineage_group']}
    saved_features = rows(output / 'features.csv')
    if len(saved_features) != 9 or {row['isolate_id'] for row in saved_features} != set(labels):
        raise ValueError('Feature identity coverage changed')
    for row in saved_features:
        expected = feature_rows[row['isolate_id']]
        if row['lineage_group'] != expected['lineage_group']:
            raise ValueError('Lineage identity changed')
        for key in ['initial_od600', 'early_od600', 'early_od_fold', 'log_early_od_fold', 'early_recorded_time_h', 'early_points']:
            close(row[key], expected[key], key)
        close(row['titer_g_L'], labels[row['isolate_id']], 'source titer')
    decisions = rows(output / 'candidate_decisions.csv')
    metric_rows = rows(output / 'lineage_metrics.csv')
    if len(decisions) != 36 or len(metric_rows) != 12:
        raise ValueError('Fold result coverage changed')
    models = json.loads((output / 'fitted_models.json').read_text(encoding='utf-8'))
    metrics, all_predictions = {}, {}
    expected_decision_keys = set()
    for group in sorted({row['lineage_group'] for row in feature_rows.values()}):
        train = sorted(identity for identity, row in feature_rows.items() if row['lineage_group'] != group)
        test = sorted(set(labels) - set(train))
        xtrain = np.array([[feature_rows[identity][key] for key in protocol['features']] for identity in train])
        xtest = np.array([[feature_rows[identity][key] for key in protocol['features']] for identity in test])
        center = xtrain.mean(axis=0)
        scales = xtrain.std(axis=0, ddof=0)
        scales[scales <= 1e-12] = 1.0
        matrix = np.column_stack([np.ones(6), (xtrain - center) / scales])
        penalty = np.diag([0.0, *([math.sqrt(protocol['ridge_penalty'])] * len(protocol['features']))])
        beta = np.linalg.lstsq(np.vstack([matrix, penalty]), np.r_[[labels[identity] for identity in train], np.zeros(3)], rcond=None)[0]
        prediction = np.column_stack([np.ones(3), (xtest - center) / scales]) @ beta
        model = models[group]
        if model['training_isolate_ids'] != train or set(model['training_lineage_groups']) != {feature_rows[identity]['lineage_group'] for identity in train}:
            raise ValueError('Training group isolation failed')
        for key, values in [('means', center), ('scales', scales), ('coefficients', beta)]:
            for actual, expected in zip(model[key], values):
                close(actual, expected, key)
        scores = {'uniform_random_expectation': np.zeros(3), 'early_density': np.array([feature_rows[i]['early_od600'] for i in test]),
            'early_density_fold': np.array([feature_rows[i]['early_od_fold'] for i in test]), 'fixed_ridge': prediction}
        all_predictions.update(dict(zip(test, prediction.tolist())))
        targets = np.array([labels[identity] for identity in test])
        for method, values in scores.items():
            tie = np.abs(values - values.max()) <= protocol['score_tie_absolute_tolerance']
            probability = tie / tie.sum()
            for identity, score, weight in zip(test, values, probability):
                expected_decision_keys.add((group, method, identity))
                matches = [row for row in decisions if (row['heldout_lineage'], row['method'], row['isolate_id']) == (group, method, identity)]
                if len(matches) != 1 or matches[0]['test_titer_used_to_score'] != 'False':
                    raise ValueError('Decision coverage or target access flag changed')
                close(matches[0]['candidate_score'], score, 'candidate score')
                close(matches[0]['selection_probability'], weight, 'selection probability')
                close(matches[0]['target_assay_budget'], 1, 'same assay budget')
                close(matches[0]['common_historical_titer_labels'], 6, 'same training history')
            expected_titer = float(probability @ targets)
            best = float(targets.max())
            value = {'expected_selected_titer_g_L': expected_titer, 'fraction_of_heldout_best': expected_titer / best,
                'absolute_regret_g_L': best - expected_titer,
                'top1_hit_probability': float(probability @ (np.abs(targets - best) <= 1e-12))}
            matches = [row for row in metric_rows if row['heldout_lineage'] == group and row['method'] == method]
            if len(matches) != 1:
                raise ValueError('Metric coverage changed')
            for key, expected in value.items():
                close(matches[0][key], expected, key)
            metrics[(group, method)] = value
    if {(row['heldout_lineage'], row['method'], row['isolate_id']) for row in decisions} != expected_decision_keys:
        raise ValueError('Unexpected decisions')
    summary = json.loads((output / 'summary.json').read_text(encoding='utf-8'))
    if summary != report['summary']:
        raise ValueError('Summary differs from report')
    for method in protocol['methods']:
        for key in ['expected_selected_titer_g_L', 'fraction_of_heldout_best', 'absolute_regret_g_L', 'top1_hit_probability']:
            close(summary[method][key], np.mean([value[key] for (_, m), value in metrics.items() if m == method]), 'lineage equal mean')
    close(summary['fixed_ridge']['prediction_mae_g_L'], np.mean([abs(all_predictions[i] - labels[i]) for i in labels]), 'prediction MAE')
    receipt = {'passed': True, 'independent_fit': 'NumPy augmented least squares; production uses scalar normal equations',
        'source_feature_rows_recomputed': 9, 'lineage_folds_refitted': 3,
        'candidate_decisions_recomputed': 36, 'fold_metrics_recomputed': 12,
        'report_sha256': digest(output / 'report.json'), 'new_blind_validation': False,
        'prospective_selection_validated': False, 'productivity_scored': False}
    (output / 'verification.json').write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    return receipt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT / 'research/results/ale_lineage_screen_v1')
    args = parser.parse_args()
    print(json.dumps(verify(args.output), ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
