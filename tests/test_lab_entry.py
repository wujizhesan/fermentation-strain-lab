import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import lab
import project


def fixture(root):
    (root / "core").mkdir()
    (root / "core/run_receipt.json").write_text("{}", encoding="utf-8")
    module = {"title": "早期数据选株", "summary": {"endpoint": "g/L", "conclusion": "未验证收益", "methods": []},
              "integrity_passed": True, "evidence": [{"report_sha256": "fixture"}], "attachments": []}
    data = {"mode": "fresh_replay", "modules": {"screening": copy.deepcopy(module)},
            "same_strain_end_to_end_validated": False, "strain_optimization_gain_validated": False,
            "original_model_external_validation": False}
    receipt = {"status": "passed", "modules_requested": ["screening"], "external_steps": [],
               "core_run_receipt_sha256": lab.digest(root / "core/run_receipt.json")}
    save(root, data, receipt)
    return data, receipt, {"overall_integrity_passed": True, "modules": {"screening": module}}


def save(root, data, receipt):
    (root / "summary.json").write_text(json.dumps(data), encoding="utf-8")
    (root / "index.html").write_text(lab.render(data), encoding="utf-8")
    receipt["artifacts"] = {name: lab.digest(root / name) for name in ("summary.json", "index.html")}
    (root / "lab_receipt.json").write_text(json.dumps(receipt), encoding="utf-8")


class LabEntryTests(unittest.TestCase):
    def test_reference_preserves_negative_results_and_snapshot_status(self):
        data = lab.reference_status()
        self.assertEqual(set(data["modules"]), set(lab.MODULES))
        self.assertTrue(all(row["numeric_recomputed_now"] is False for row in data["modules"].values()))
        self.assertIn("未高于", data["modules"]["screening"]["summary"]["conclusion"])
        self.assertFalse(data["original_model_external_validation"])

    def test_summary_cannot_change_metrics_by_resigning_overview(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data, receipt, core = fixture(root)
            data["modules"]["screening"]["summary"]["conclusion"] = "已验证增产"
            save(root, data, receipt)
            with patch.object(project, "collect_run", return_value=core), self.assertRaisesRegex(ValueError, "总览与模块结果不一致"):
                lab.collect_run(root)

    def test_failed_core_cannot_be_wrapped_as_passed(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _, _, core = fixture(root)
            core["overall_integrity_passed"] = False
            with patch.object(project, "collect_run", return_value=core), self.assertRaisesRegex(ValueError, "核心模块证据检查失败"):
                lab.collect_run(root)

    def test_html_cannot_disagree_with_summary_after_resigning_hash(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            _, receipt, core = fixture(root)
            (root / "index.html").write_text("fake positive result", encoding="utf-8")
            receipt["artifacts"]["index.html"] = lab.digest(root / "index.html")
            (root / "lab_receipt.json").write_text(json.dumps(receipt), encoding="utf-8")
            with patch.object(project, "collect_run", return_value=core), self.assertRaisesRegex(ValueError, "HTML总览"):
                lab.collect_run(root)

    def test_missing_external_verification_step_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            data, receipt, core = fixture(root)
            receipt["modules_requested"].append("external")
            save(root, data, receipt)
            with patch.object(project, "collect_run", return_value=core), self.assertRaisesRegex(ValueError, "核验步骤缺失"):
                lab.collect_run(root)

    def test_failed_external_command_retains_failure_log(self):
        def fail(interpreter, args, output, label, records):
            log = output / "logs/failure.log"
            log.parent.mkdir()
            log.write_text("failed fixture", encoding="utf-8")
            records.append({"step": label, "exit_code": 7, "log": str(log), "log_sha256": lab.digest(log)})
            raise RuntimeError("failed fixture")

        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "new"
            with patch.object(project, "preflight"), patch.object(project, "execute_step", side_effect=fail):
                with self.assertRaisesRegex(RuntimeError, "failed fixture"):
                    lab.run(output, ("external",), sys.executable)
            receipt = project.read_json(output / "lab_receipt.json")
            self.assertEqual(receipt["status"], "failed")
            self.assertEqual(receipt["external_steps"][0]["exit_code"], 7)
            self.assertTrue((output / receipt["external_steps"][0]["log"]).is_file())
            self.assertFalse((output / "index.html").exists())

    def test_existing_output_is_rejected_before_execution(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch.object(project, "preflight") as preflight, self.assertRaisesRegex(ValueError, "新目录"):
                lab.run(Path(directory), ("screening",), sys.executable)
            preflight.assert_not_called()

    def test_failed_run_cannot_be_reused(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "lab_receipt.json").write_text('{"status":"failed"}', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "运行未完成"):
                lab.collect_run(root)


if __name__ == "__main__":
    unittest.main()
