import csv
import json
import os
import unittest
from pathlib import Path

import ale_marge_demo as demo


def result_root():
    return Path(os.environ.get('MARGE_RESULTS_ROOT', demo.ROOT / 'research/results'))


class MargeDemoTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.protocol = json.loads(demo.PROTOCOL.read_text(encoding='utf-8'))
        cls.models = {design: demo.prepared_model(cls.protocol, design) for design in cls.protocol['designs']}

    def synthetic_protocol(self, rows, name):
        folder = demo.ROOT / 'tmp' / 'test_ale_marge_demo'
        folder.mkdir(parents=True, exist_ok=True)
        path = folder / name
        with path.open('w', encoding='utf-8', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=['omics', 'design', 'feature_id', 'systematic_id',
                'released_log2_ratio_evolved_over_parental', 'parental_n_positive', 'evolved_n_positive'])
            writer.writeheader(); writer.writerows(rows)
        return {**self.protocol, 'abundance_path': str(path)}

    def test_strict_paper_entry_stays_closed(self):
        status = demo.strict_status(self.protocol)
        self.assertFalse(status['strict_entry_allowed'])
        self.assertFalse(status['paper_exact_reproduction_completed'])
        self.assertFalse(status['extracellular_rate_constraints_reproduced'])
        self.assertEqual(status['missing_required_inputs'], self.protocol['missing_paper_inputs'])
        self.assertGreaterEqual(len(self.protocol['missing_paper_inputs']), 5)

    def test_all_limitations_are_declared(self):
        self.assertGreaterEqual(len(self.protocol['limitations']), 10)
        self.assertTrue(any('q<0.1' in item for item in self.protocol['limitations']))
        self.assertTrue(any('fitted separately' in item for item in self.protocol['limitations']))

    def test_proteome_balance_override_is_declared(self):
        override = self.protocol['proteome_balance_override']
        self.assertTrue(override['expected'])
        self.assertEqual(override['row_ids'], ['proteome_a', 'proteome_b'])
        self.assertEqual(override['expected_terms_per_row'], 1)
        self.assertTrue(any('capacity row under the same name' in item for item in self.protocol['limitations']))
        self.assertTrue(any('proteome budget does not constrain fluxes' in item for item in self.protocol['limitations']))

    def test_proteome_balance_row_collides_with_capacity_row_id(self):
        from reframed.core.transformation import gpr_transform
        for design, model in self.models.items():
            with self.subTest(design=design):
                transformed = gpr_transform(model, inplace=False, add_proteome=True)
                balance = transformed.metabolite_reaction_lookup()['proteome']
                self.assertIn('proteome_synth', balance)
                self.assertTrue(any(key.startswith('u_') for key in balance))
                # upstream names the joint balance row '<metabolite>_<condition>', i.e. 'proteome_a',
                # and later adds the capacity row under that same id
                self.assertEqual('proteome' + '_a', 'proteome_a')
                self.assertEqual('proteome' + '_b', 'proteome_b')

    def test_archived_cases_record_the_degeneracy(self):
        output = result_root() / 'ale_marge_demo_v1'
        if not (output / 'report.json').exists():
            self.skipTest('archived MARGE output not built')
        report = json.loads((output / 'report.json').read_text(encoding='utf-8'))
        self.assertFalse(report['proteome_budget_effective'])
        self.assertFalse(report['stage2_carries_information'])
        self.assertEqual(report['proteome_balance_override_confirmed_cases'], len(report['cases']))
        for summary in report['cases']:
            case = json.loads((output / summary['case_id'] / 'case.json').read_text(encoding='utf-8'))
            with self.subTest(case=summary['case_id']):
                self.assertTrue(case['proteome_mass_balance_absent'])
                self.assertFalse(case['proteome_budget_constrains_fluxes'])
                self.assertTrue(case['stage2_degenerate'])
                self.assertEqual(case['proteome_row_terms_a'], 1)
                self.assertEqual(case['proteome_row_terms_b'], 1)
                self.assertEqual(case['stage2_proteome_optimum'], 0.0)
                for condition in ['a', 'b']:
                    self.assertGreater(case['sum_gene_usage_' + condition], 1e-6)
                    self.assertAlmostEqual(case['proteome_synth_' + condition], 0.0)

    def test_balance_preserving_variant_restores_the_accounting(self):
        output = result_root() / 'ale_marge_balancefix_v1'
        if not (output / 'report.json').exists():
            self.skipTest('balance-preserving variant not built')
        report = json.loads((output / 'report.json').read_text(encoding='utf-8'))
        self.assertTrue(report['proteome_balance_kept'])
        self.assertTrue(report['proteome_budget_effective'])
        self.assertTrue(report['stage2_carries_information'])
        self.assertEqual(report['proteome_balance_override_confirmed_cases'], 0)
        self.assertEqual(report['proteome_capacity_binding_cases'], 0)
        for summary in report['cases']:
            case = json.loads((output / summary['case_id'] / 'case.json').read_text(encoding='utf-8'))
            with self.subTest(case=summary['case_id']):
                self.assertFalse(case['proteome_mass_balance_absent'])
                self.assertTrue(case['proteome_budget_constrains_fluxes'])
                self.assertFalse(case['stage2_degenerate'])
                self.assertGreater(case['stage2_proteome_optimum'], 1e-6)
                for condition in ['a', 'b']:
                    usage = case['sum_gene_usage_' + condition]
                    self.assertAlmostEqual(case['proteome_synth_' + condition], usage,
                        delta=1e-6 * max(1.0, usage))
                    self.assertGreater(case['proteome_capacity_slack_' + condition], 0.0)

    def test_fix_leaves_released_ranges_but_not_neutral_ranges_alone(self):
        baseline = result_root() / 'ale_marge_demo_v1'
        variant = result_root() / 'ale_marge_balancefix_v1'
        if not (baseline / 'report.json').exists() or not (variant / 'report.json').exists():
            self.skipTest('both MARGE variants must be built')

        def ranges(root_dir, name):
            with (root_dir / name / 'flux_ranges.csv').open(encoding='utf-8', newline='') as handle:
                rows = list(csv.DictReader(handle))
            return {(r['condition'], r['reaction_id']): (float(r['minimum']), float(r['maximum'])) for r in rows}

        for mode, expect_change in [('released_ratio', False), ('neutral_ratio', True)]:
            with self.subTest(mode=mode):
                moved = False
                for name in ['F_protein_' + mode, 'Sz_protein_' + mode, 'Mz_protein_' + mode]:
                    before, after = ranges(baseline, name), ranges(variant, name)
                    for key in before:
                        delta = max(abs(before[key][0] - after[key][0]), abs(before[key][1] - after[key][1]))
                        if delta > 1e-5:
                            moved = True
                self.assertEqual(moved, expect_change, f'range behaviour changed for {mode}')

    def test_source_bounds_and_maintenance_are_preserved(self):
        for design, model in self.models.items():
            with self.subTest(design=design):
                self.assertEqual(model.reactions['R_ATPM'].lb, 1)
                self.assertEqual(model.reactions['R_ATPM'].ub, 1)
                self.assertAlmostEqual(model.reactions['R_EX_glc_e_'].lb, -1.15)
                self.assertEqual((model.reactions['R_EX_gly_e_'].lb, model.reactions['R_EX_gly_e_'].ub), (0.0, 1000.0))

    def test_glycine_uptake_is_opened_only_under_a_hypothetical_cap(self):
        self.assertTrue(any('hypothetical' in item for item in self.protocol['limitations']))
        self.assertEqual(self.protocol['parental_glycine_uptake_cap'], 1.0)
        self.assertEqual(self.protocol['evolved_glycine_uptake_cap'], 0.0)
        for design, model in self.models.items():
            with self.subTest(design=design):
                self.assertEqual(model.reactions['R_EX_gly_e_'].lb, 0.0)

    def test_isoenzyme_knockout_follows_gpr_logic(self):
        model = self.models['F']
        self.assertEqual((model.reactions['R_PGCD'].lb, model.reactions['R_PGCD'].ub), (0, 0))
        self.assertEqual((model.reactions['R_FUM'].lb, model.reactions['R_FUM'].ub), (0, 0))
        self.assertEqual((model.reactions['R_FUMm'].lb, model.reactions['R_FUMm'].ub), (0, 0))

    def test_knockout_sets_are_design_specific(self):
        sz, mz = self.models['Sz'], self.models['Mz']
        self.assertNotEqual((sz.reactions['R_FUMm'].lb, sz.reactions['R_FUMm'].ub), (0, 0))
        self.assertEqual((sz.reactions['R_G6PDH2'].lb, sz.reactions['R_G6PDH2'].ub), (0, 0))
        self.assertEqual((mz.reactions['R_MDHm'].lb, mz.reactions['R_MDHm'].ub), (0, 0))
        self.assertEqual((mz.reactions['R_MDH'].lb, mz.reactions['R_MDH'].ub), (0, 0))
        self.assertEqual((mz.reactions['R_ME2m'].lb, mz.reactions['R_ME2m'].ub), (0, 0))

    def test_ambiguous_genotype_gene_is_rejected(self):
        folder = demo.ROOT / 'tmp' / 'test_ale_marge_demo'
        folder.mkdir(parents=True, exist_ok=True)
        table = folder / 'ambiguous_gene_table.tsv'
        table.write_text('YER081W\tSER3\nYIL074C\tSER3\n', encoding='utf-8')
        protocol = {**self.protocol, 'gene_name_table': str(table)}
        with self.assertRaises(ValueError):
            demo.prepared_model(protocol, 'F')

    def test_ratios_are_positive_two_to_the_log2(self):
        protocol = self.protocol
        accepted, excluded = demo.input_ratios(protocol, self.models['F'], 'protein', 'F')
        self.assertEqual(len(accepted), 579)
        self.assertTrue(excluded)
        for row in accepted:
            self.assertGreater(row['relative_expression_b_over_a'], 0)
            self.assertAlmostEqual(row['relative_expression_b_over_a'], 2 ** row['released_log2_ratio'])
            self.assertFalse(row['q_value_available'])
            self.assertFalse(row['paper_q_filter_applied'])

    def test_duplicate_annotation_is_excluded_not_averaged(self):
        protocol = self.synthetic_protocol([
            {'omics': 'protein', 'design': 'F', 'feature_id': 'f1', 'systematic_id': 'YER081W',
             'released_log2_ratio_evolved_over_parental': '1', 'parental_n_positive': '3', 'evolved_n_positive': '3'},
            {'omics': 'protein', 'design': 'F', 'feature_id': 'f2', 'systematic_id': 'YER081W',
             'released_log2_ratio_evolved_over_parental': '2', 'parental_n_positive': '3', 'evolved_n_positive': '3'},
        ], 'duplicate_annotation.csv')
        accepted, excluded = demo.input_ratios(protocol, self.models['F'], 'protein', 'F')
        self.assertEqual(accepted, [])
        self.assertEqual({row['reason'] for row in excluded}, {'multiple_released_features_map_to_same_gene'})

    def test_missing_ratio_is_not_imputed(self):
        protocol = self.synthetic_protocol([
            {'omics': 'protein', 'design': 'F', 'feature_id': 'f1', 'systematic_id': 'YER081W',
             'released_log2_ratio_evolved_over_parental': '', 'parental_n_positive': '1', 'evolved_n_positive': '2'},
        ], 'missing_ratio.csv')
        accepted, excluded = demo.input_ratios(protocol, self.models['F'], 'protein', 'F')
        self.assertEqual(accepted, [])
        self.assertEqual([row['reason'] for row in excluded], ['fewer_than_two_positive_observations_in_a_group'])

    def test_pseudo_and_unmapped_features_are_excluded(self):
        protocol = self.synthetic_protocol([
            {'omics': 'protein', 'design': 'F', 'feature_id': 'f1', 'systematic_id': 'S0001',
             'released_log2_ratio_evolved_over_parental': '0.5', 'parental_n_positive': '3', 'evolved_n_positive': '3'},
            {'omics': 'protein', 'design': 'F', 'feature_id': 'f2', 'systematic_id': '',
             'released_log2_ratio_evolved_over_parental': '0.5', 'parental_n_positive': '3', 'evolved_n_positive': '3'},
            {'omics': 'protein', 'design': 'F', 'feature_id': 'f3', 'systematic_id': 'YZZ999W',
             'released_log2_ratio_evolved_over_parental': '0.5', 'parental_n_positive': '3', 'evolved_n_positive': '3'},
        ], 'pseudo_features.csv')
        accepted, excluded = demo.input_ratios(protocol, self.models['F'], 'protein', 'F')
        self.assertEqual(accepted, [])
        self.assertEqual({row['reason'] for row in excluded}, {'no_model_gene_or_pseudo_gene'})

    def test_nonfinite_ratio_is_rejected(self):
        protocol = self.synthetic_protocol([
            {'omics': 'protein', 'design': 'F', 'feature_id': 'f1', 'systematic_id': 'YER081W',
             'released_log2_ratio_evolved_over_parental': '5000', 'parental_n_positive': '3', 'evolved_n_positive': '3'},
        ], 'nonfinite_ratio.csv')
        accepted, excluded = demo.input_ratios(protocol, self.models['F'], 'protein', 'F')
        self.assertEqual(accepted, [])
        self.assertEqual([row['reason'] for row in excluded], ['invalid_positive_ratio'])

    def test_underflowing_ratio_is_rejected_as_zero(self):
        protocol = self.synthetic_protocol([
            {'omics': 'protein', 'design': 'F', 'feature_id': 'f1', 'systematic_id': 'YER081W',
             'released_log2_ratio_evolved_over_parental': '-5000', 'parental_n_positive': '3', 'evolved_n_positive': '3'},
        ], 'underflowing_ratio.csv')
        accepted, excluded = demo.input_ratios(protocol, self.models['F'], 'protein', 'F')
        self.assertEqual(accepted, [])
        self.assertEqual([row['reason'] for row in excluded], ['invalid_positive_ratio'])

    def test_growth_fractions_follow_mu_ratio(self):
        for design in self.protocol['designs']:
            rates = self.protocol['growth_rates_h_inverse']
            with self.subTest(design=design):
                self.assertAlmostEqual(rates[design]['parental'] / rates['WT'], rates[design]['parental'] / 0.332)
                self.assertGreater(rates[design]['evolved'], rates[design]['parental'])

    def test_modes_are_released_and_neutral_control(self):
        self.assertEqual(sorted(self.protocol['modes']), ['neutral_ratio', 'released_ratio'])
        self.assertEqual(self.protocol['step2_relax'], 0.1)

    def test_tracked_reactions_track_product_exporter(self):
        expected = {'F': 'R_EX_fum_e_', 'Sz': 'R_EX_succ_e_', 'Mz': 'R_EX_mal_L_e_'}
        for design, exporter in expected.items():
            with self.subTest(design=design):
                tracked = self.protocol['tracked_reactions'][design]
                self.assertEqual(len(tracked), 5)
                self.assertIn(exporter, tracked)
                self.assertIn('R_MDHp', tracked)

    def test_upstream_module_fingerprints_are_pinned(self):
        self.assertTrue(self.protocol['upstream_module_sha256'])
        self.assertTrue(all(len(value) == 64 for value in self.protocol['upstream_module_sha256'].values()))
        self.assertEqual(self.protocol['reframed_version'], '1.2.1')
        self.assertTrue(self.protocol['input_sha256'])


if __name__ == '__main__':
    unittest.main()
