import copy
import unittest

from ale_lineage_screen import build_features, evaluate_decisions, fold_decisions, predict, selection_weights


PROTOCOL = {'design': 'F', 'product': 'fumarate', 'cutoff_h': 18.0,
    'features': ['initial_od600', 'log_early_od_fold'], 'ridge_penalty': 1.0,
    'score_tie_absolute_tolerance': 1e-12,
    'methods': ['uniform_random_expectation', 'early_density', 'early_density_fold', 'fixed_ridge']}


def fixture():
    observations, endpoints = [], []
    for group in [1, 2, 3]:
        for index, clone in enumerate(['S', 'M', 'L']):
            identity = f'Chassis_F-{group}-{clone}'
            initial = 0.1 + group * 0.01 + index * 0.001
            for row_index, (time, value) in enumerate([(0, initial), (18, initial * (2 + group * 0.2 + index * 0.3)), (30, initial * 20)], 3):
                observations.append({'isolate_id': identity, 'recorded_time_h': time, 'od600': value, 'source_excel_row': row_index})
            endpoints.append({'design': 'F', 'isolate_id': identity, 'lineage_group': f'Chassis_F-{group}', 'product': 'fumarate',
                'titer_g_L': 0.2 + group * 0.01 + index * 0.04, 'source_sheet': 'fixture', 'source_excel_row': group * 3 + index,
                'assay_collection_time_h': '', 'assay_return_time_verified': 'False'})
    return observations, endpoints


class ALELineageScreenTests(unittest.TestCase):
    def test_full_population_is_held_out(self):
        records = build_features(*fixture(), PROTOCOL)
        model, decisions = fold_decisions(records, 'Chassis_F-2', PROTOCOL)
        self.assertEqual(len(model['training_isolate_ids']), 6)
        self.assertEqual(model['training_lineage_groups'], ['Chassis_F-1', 'Chassis_F-3'])
        self.assertEqual({row['lineage_group'] for row in decisions}, {'Chassis_F-2'})
        self.assertFalse(set(model['training_isolate_ids']) & {row['isolate_id'] for row in decisions})

    def test_heldout_targets_do_not_change_prediction_or_selection(self):
        records = build_features(*fixture(), PROTOCOL)
        changed = copy.deepcopy(records)
        for row in changed:
            if row['lineage_group'] == 'Chassis_F-2':
                row['titer_g_L'] *= 1000
        self.assertEqual(fold_decisions(records, 'Chassis_F-2', PROTOCOL), fold_decisions(changed, 'Chassis_F-2', PROTOCOL))

    def test_heldout_features_do_not_change_training_preprocessing(self):
        records = build_features(*fixture(), PROTOCOL)
        changed = copy.deepcopy(records)
        for row in changed:
            if row['lineage_group'] == 'Chassis_F-2':
                row['initial_od600'] *= 1000
                row['log_early_od_fold'] += 100
        before, original = fold_decisions(records, 'Chassis_F-2', PROTOCOL)
        after, perturbed = fold_decisions(changed, 'Chassis_F-2', PROTOCOL)
        self.assertEqual(before, after)
        self.assertNotEqual([row['candidate_score'] for row in original if row['method'] == 'fixed_ridge'],
            [row['candidate_score'] for row in perturbed if row['method'] == 'fixed_ridge'])

    def test_future_od_is_excluded(self):
        observations, endpoints = fixture()
        changed = copy.deepcopy(observations)
        for row in changed:
            if row['recorded_time_h'] > 18:
                row['od600'] *= 1000
        self.assertEqual(build_features(observations, endpoints, PROTOCOL), build_features(changed, endpoints, PROTOCOL))

    def test_uniform_ties_and_random_baseline_keep_budget_one(self):
        self.assertEqual(selection_weights([7, 7, 1], 1e-12), [0.5, 0.5, 0])
        records = build_features(*fixture(), PROTOCOL)
        _, decisions = fold_decisions(records, 'Chassis_F-2', PROTOCOL)
        metrics = evaluate_decisions(decisions, records, PROTOCOL)
        random = next(row for row in metrics if row['method'] == 'uniform_random_expectation')
        candidates = [row for row in records if row['lineage_group'] == 'Chassis_F-2']
        self.assertAlmostEqual(random['expected_selected_titer_g_L'], sum(row['titer_g_L'] for row in candidates) / 3)
        self.assertAlmostEqual(random['top1_hit_probability'], 1 / 3)
        self.assertEqual(random['target_assay_budget'], 1)

    def test_negative_predictions_are_preserved_without_titer_claim(self):
        model = {'features': ['initial_od600'], 'means': [0], 'scales': [1], 'coefficients': [-1, 0]}
        self.assertEqual(predict(model, {'initial_od600': 1}), -1)

    def test_missing_initial_and_duplicate_isolates_are_rejected(self):
        observations, endpoints = fixture()
        changed = copy.deepcopy(observations)
        changed[0]['recorded_time_h'] = 1
        with self.assertRaisesRegex(ValueError, 'Initial OD'):
            build_features(changed, endpoints, PROTOCOL)
        with self.assertRaisesRegex(ValueError, 'distinct fumarate isolates'):
            build_features(observations, endpoints + [endpoints[0]], PROTOCOL)

    def test_productivity_clock_is_not_invented(self):
        records = build_features(*fixture(), PROTOCOL)
        self.assertTrue(all(row['assay_collection_time_h'] == '' for row in records))
        self.assertTrue(all(row['assay_return_time_verified'] == 'False' for row in records))
        self.assertTrue(all('productivity_g_L_h' not in row for row in records))


if __name__ == '__main__':
    unittest.main()
