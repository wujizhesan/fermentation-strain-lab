import copy
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
import unittest

from build_marge_viewer import VARIANTS, load_viewer_data, render_viewer


ROOT = Path(__file__).resolve().parents[1]
RESULTS = Path(os.environ.get('MARGE_RESULTS_ROOT', ROOT / 'research/results'))


def copy_display_inputs(destination):
    for name in VARIANTS.values():
        source = RESULTS / name
        folder = destination / name
        folder.mkdir()
        report = json.loads((source / 'report.json').read_text(encoding='utf-8'))
        files = ['report.json', 'verification.json', 'strict_reproduction_status.json', 'source_audit.json', 'input_gene_ratios.csv']
        files.extend(case['case_id'] + '/flux_ranges.csv' for case in report['cases'])
        for name in files:
            path = folder / name
            path.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source / name, path)


class MargeViewerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = load_viewer_data(RESULTS)

    def test_verified_data_preserves_the_widening_counterexample(self):
        self.assertEqual(self.data['summary'], {'cases': 24, 'joint_lps': 48, 'range_endpoints': 480, 'input_gene_ratios': 4159})
        for variant in self.data['variants'].values():
            self.assertEqual(variant['width_counts'], {'narrower': 59, 'same': 0, 'wider': 1})
            case = variant['wider_cases'][0]
            self.assertEqual((case['design'], case['layer'], case['condition'], case['reaction_id']), ('Sz', 'protein', 'a', 'R_EX_succ_e_'))
            self.assertAlmostEqual(case['neutral_width'], 1.0956104988966, places=7)
            self.assertAlmostEqual(case['released_width'], 2.0313152709923, places=7)

    def test_changed_display_csv_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            copy_display_inputs(base)
            path = base / 'ale_marge_demo_v1/F_protein_released_ratio/flux_ranges.csv'
            path.write_bytes(path.read_bytes().replace(b'1.694500', b'9.694500', 1))
            with self.assertRaisesRegex(ValueError, 'artifact fingerprint changed'):
                load_viewer_data(base)

    def test_stale_verification_cannot_be_displayed_as_passed(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            copy_display_inputs(base)
            path = base / 'ale_marge_demo_v1/verification.json'
            value = json.loads(path.read_text(encoding='utf-8'))
            value['report_sha256'] = 'stale'
            path.write_text(json.dumps(value), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'verified report required'):
                load_viewer_data(base)

    def test_resigned_duplicate_case_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            copy_display_inputs(base)
            folder = base / 'ale_marge_demo_v1'
            report = json.loads((folder / 'report.json').read_text(encoding='utf-8'))
            report['cases'].append(copy.deepcopy(report['cases'][0]))
            body = json.dumps(report).encode('utf-8')
            (folder / 'report.json').write_bytes(body)
            verification = json.loads((folder / 'verification.json').read_text(encoding='utf-8'))
            verification['report_sha256'] = hashlib.sha256(body).hexdigest()
            (folder / 'verification.json').write_text(json.dumps(verification), encoding='utf-8')
            with self.assertRaisesRegex(ValueError, 'case coverage'):
                load_viewer_data(base)

    def test_embedded_data_round_trips_special_feature_names(self):
        data = copy.deepcopy(self.data)
        data['genes'][0]['feature_id'] = '<MDH & </script> α>'
        page = render_viewer(data)
        encoded = page.split('<script type="application/json" id="marge-data">', 1)[1].split('</script>', 1)[0]
        self.assertEqual(json.loads(encoded), data)
        self.assertNotIn('</script>', encoded)


if __name__ == '__main__':
    unittest.main()
