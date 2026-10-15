from __future__ import annotations

import copy
import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import openpyxl
import external_screen_validation as producer
import verify_external_screen_validation as validator

ROOT = Path(__file__).resolve().parent.parent
RESULT = ROOT / "research/results/ale_external_screen_v1/report.json"


class ExternalScreenValidationTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.protocol = json.loads(producer.PROTOCOL.read_text(encoding="utf-8"))
        cls.report = json.loads(RESULT.read_text(encoding="utf-8"))

    def test_current_sources_and_replay_verify(self):
        checks = validator.verify(RESULT)
        self.assertTrue(all(checks.values()), checks)

    def test_endpoint_changes_cannot_change_early_decisions(self):
        for spec in self.protocol["panels"]:
            rows = [r for r in self.report["candidates"] if r["panel_id"] == spec["id"]]
            changed = copy.deepcopy(rows)
            for i, row in enumerate(changed):
                row["target_mean"] = 100 * (i + 1)
            before, _ = producer.score_panel(rows, spec, 1e-12)
            after, _ = producer.score_panel(changed, spec, 1e-12)
            self.assertEqual(before, after)

    def test_ties_use_uniform_expectation(self):
        self.assertEqual(producer.select_weights([1, 1, 0], 1e-12), [0.5, 0.5, 0])
        self.assertEqual(producer.select_weights([-2, -1, -1], 1e-12), [0, 0.5, 0.5])
        with self.assertRaises(ValueError):
            producer.select_weights([float("nan")], 1e-12)

    def test_identity_missing_values_and_future_time_rejected(self):
        spec = copy.deepcopy(self.protocol["panels"][0])
        for case in ("identity", "missing", "future", "duplicate_rep", "duplicate_strain"):
            with self.subTest(case=case):
                book = openpyxl.load_workbook(ROOT / spec["source"], data_only=False)
                altered = copy.deepcopy(spec)
                sheet = book[spec["early_sheet"]]
                if case == "identity":
                    sheet["E1"] = "unrelated_strain"
                elif case == "missing":
                    sheet["E3"] = None
                elif case == "future":
                    sheet["A3"] = 48
                elif case == "duplicate_rep":
                    altered["candidates"][0]["columns"] = ["E", "E", "G"]
                else:
                    altered["candidates"][1]["isolate_id"] = altered["candidates"][0]["isolate_id"]
                with self.assertRaises(ValueError):
                    producer.extract_panel(book, altered, 2)
                book.close()

    def test_raw_parser_keeps_missing_values_and_negative_signals(self):
        path = ROOT / self.protocol["panels"][1]["source"]
        cells = validator.sheet_cells(path, "Sheet1")
        self.assertEqual(cells["K17"], -0.00373)
        source = ROOT / "data/raw/ale_external_sources_20261015/PMC12663439/41467_2025_65793_MOESM14_ESM.xlsx"
        cells = validator.sheet_cells(source, "Fig. S5")
        self.assertIsNone(cells.get("C8"))

    def test_re_signed_wrong_reports_rejected(self):
        for case in ("identity", "early", "future", "decision", "metric", "scope", "unit", "replicates"):
            with self.subTest(case=case), tempfile.TemporaryDirectory() as directory:
                folder = Path(directory)
                for file in RESULT.parent.glob("*.csv"):
                    (folder / file.name).write_bytes(file.read_bytes())
                doc = copy.deepcopy(self.report)
                if case == "identity":
                    doc["candidates"][0]["isolate_id"] = "anonymous_1"
                elif case == "early":
                    doc["candidates"][0]["early_mean"] = doc["candidates"][0]["target_mean"]
                elif case == "future":
                    doc["candidates"][0]["cutoff_h"] = 96
                elif case == "decision":
                    doc["decisions"][0]["selection_probability"] = 1
                elif case == "metric":
                    doc["metrics"][1]["fraction_of_panel_best"] = 1
                elif case == "scope":
                    doc["scope"]["original_18h_fumarate_ridge_external_validation"] = True
                elif case == "unit":
                    doc["metrics"][-1]["target_unit"] = "g/L"
                else:
                    doc["candidates"][0]["replicate_count"] = 30
                path = folder / "report.json"
                path.write_text(json.dumps(doc, ensure_ascii=False), encoding="utf-8")
                checks = validator.verify(path)
                self.assertFalse(all(checks.values()), checks)
                self.assertTrue(checks["协议 SHA"])
                self.assertTrue(checks["原始输入 SHA"])

    def test_changed_csv_rejected_after_re_signing(self):
        with tempfile.TemporaryDirectory() as directory:
            folder = Path(directory)
            for file in RESULT.parent.glob("*.csv"):
                (folder / file.name).write_bytes(file.read_bytes())
            path = folder / "decisions.csv"
            with path.open(encoding="utf-8", newline="") as stream:
                rows = list(csv.DictReader(stream))
            rows[0]["score"] = "999"
            with path.open("w", encoding="utf-8", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            doc = copy.deepcopy(self.report)
            doc["artifact_sha256"]["decisions.csv"] = hashlib.sha256(path.read_bytes()).hexdigest()
            result = folder / "report.json"
            result.write_text(json.dumps(doc, ensure_ascii=False), encoding="utf-8")
            checks = validator.verify(result)
            self.assertTrue(checks["CSV SHA"])
            self.assertFalse(checks["CSV 内容绑定原始单元格"])

    def test_empty_output_reexecution_is_identical(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "fresh"
            proc = subprocess.run([sys.executable, "-X", "utf8", "-B", "external_screen_validation.py", "--output", str(path)], cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(proc.returncode, 0, proc.stdout + proc.stderr)
            fresh = json.loads((path / "report.json").read_text(encoding="utf-8"))
            original = copy.deepcopy(self.report)
            fresh.pop("generated_at_utc")
            original.pop("generated_at_utc")
            self.assertEqual(fresh, original)
            self.assertTrue(all(validator.verify(path / "report.json").values()))


if __name__ == "__main__":
    unittest.main()
