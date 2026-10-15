import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import project


def fixture(root):
    root = Path(root)
    folder = root / 'result'
    folder.mkdir()
    protocol = root / 'protocol.json'
    protocol.write_text('{}\n', encoding='utf-8')
    (folder / 'values.csv').write_text('value\n1\n', encoding='utf-8')
    report = {'protocol_sha256': project.sha(protocol),
              'artifact_sha256': {'values.csv': project.sha(folder / 'values.csv')}}
    (folder / 'report.json').write_text(json.dumps(report), encoding='utf-8')
    receipt = {'passed': True, 'report_sha256': project.sha(folder / 'report.json')}
    (folder / 'verification.json').write_text(json.dumps(receipt), encoding='utf-8')
    return folder, protocol


class ProjectEntryTests(unittest.TestCase):
    def test_stale_verification_cannot_pass(self):
        with tempfile.TemporaryDirectory() as directory:
            folder, protocol = fixture(directory)
            report = project.read_json(folder / 'report.json')
            report['modified'] = True
            (folder / 'report.json').write_text(json.dumps(report), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, '回执与报告不匹配'):
                project.checked_result(folder, protocol)

    def test_changed_artifact_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            folder, protocol = fixture(directory)
            (folder / 'values.csv').write_text('value\n999\n', encoding='utf-8')
            with self.assertRaisesRegex(ValueError, '结果文件摘要不匹配'):
                project.checked_result(folder, protocol)

    def test_failed_receipt_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            folder, protocol = fixture(directory)
            (folder / 'verification.json').write_text(json.dumps({'passed': False}), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, '未通过'):
                project.checked_result(folder, protocol)

    def test_manifest_cannot_escape_result_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            for name in ['../outside.csv', '..\\outside.csv', 'C:\\outside.csv']:
                with self.subTest(name=name), self.assertRaises(ValueError):
                    project.relative_file(directory, name)

    def test_nonfinite_metrics_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'bad.json'
            path.write_text('{"score": NaN}', encoding='utf-8')
            with self.assertRaisesRegex(ValueError, '非有限'):
                project.read_json(path)
        with self.assertRaises(ValueError):
            project.relative_change('inf', 1)

    def test_output_protects_data_and_existing_scientific_archives(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'data/raw').mkdir(parents=True)
            (root / 'research/results/frozen').mkdir(parents=True)
            for path in [root / 'data/raw/new', root / 'research/results/frozen/new']:
                with self.subTest(path=path), self.assertRaisesRegex(ValueError, '历史结果内'):
                    project.check_output(path, root)
            self.assertEqual(project.check_output(root / 'research/results/new', root),
                             (root / 'research/results/new').resolve())

    def test_actual_junction_raw_target_is_protected(self):
        output = (project.ROOT / 'data/raw').resolve() / 'project_entry_guard_fixture'
        with self.assertRaisesRegex(ValueError, '历史结果内'):
            project.check_output(output)

    def test_existing_empty_output_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(ValueError, '新目录'):
                project.check_output(directory)

    def test_positive_replay_does_not_create_prospective_gain_claim(self):
        report = {'isolates': 9, 'population_groups': 3, 'summary': {}}
        for name, value in [('uniform_random_expectation', .2), ('fixed_ridge', .3)]:
            report['summary'][name] = {'expected_selected_titer_g_L': value,
                                      'fraction_of_heldout_best': .8, 'top1_hit_probability': .3}
        summary = project.screening_summary(report)
        self.assertFalse(summary['prospective_selection_validated'])
        self.assertFalse(summary['strain_optimization_gain_validated'])
        self.assertAlmostEqual(summary['methods'][1]['relative_to_random_percent'], 50)

    def test_zero_reference_has_no_invented_relative_gain(self):
        self.assertIsNone(project.relative_change(.3, 0))

    def test_failed_command_is_logged_and_raised(self):
        with tempfile.TemporaryDirectory() as directory:
            records = []
            with self.assertRaisesRegex(RuntimeError, '步骤失败'):
                project.execute_step(sys.executable, ['-c', 'print("failed fixture"); raise SystemExit(7)'],
                                     Path(directory), 'failure_fixture', records)
            self.assertEqual(records[0]['exit_code'], 7)
            self.assertIn('failed fixture', Path(records[0]['log']).read_text(encoding='utf-8'))

    def test_pipeline_failure_stops_following_modules_and_records_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / 'new_run'
            with patch.object(project, 'preflight', return_value={'fixture': True}), \
                 patch.object(project, 'execute_step', side_effect=RuntimeError('producer failed')) as execute:
                with self.assertRaisesRegex(RuntimeError, 'producer failed'):
                    project.run_project(output, project.MODULES, sys.executable)
            self.assertEqual(execute.call_count, 1)
            self.assertEqual(project.read_json(output / 'run_receipt.json')['status'], 'failed')
            self.assertFalse((output / 'index.html').exists())

    def test_bad_module_remains_failed_in_combined_status(self):
        with patch.object(project, 'load_module', side_effect=ValueError('bad input')):
            data = project.collect(selected=('screening',))
        self.assertFalse(data['overall_integrity_passed'])
        self.assertFalse(data['modules']['screening']['integrity_passed'])
        self.assertFalse(data['same_strain_end_to_end_validated'])

    def test_incomplete_run_cannot_be_reused_as_passed_results(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            (path / 'run_receipt.json').write_text(json.dumps({'status': 'running'}), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, '已完成'):
                project.collect_run(path)

    def test_run_with_missing_independent_verification_step_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            for filename in ['project_summary.json', 'index.html']:
                (path / filename).write_text('fixture', encoding='utf-8')
            receipt = {'status': 'passed', 'modules_requested': ['screening'], 'steps': [],
                       'project_summary_sha256': project.sha(path / 'project_summary.json'),
                       'overview_sha256': project.sha(path / 'index.html')}
            (path / 'run_receipt.json').write_text(json.dumps(receipt), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, '缺少计算、核验或展示'):
                project.collect_run(path)


if __name__ == '__main__':
    unittest.main()
