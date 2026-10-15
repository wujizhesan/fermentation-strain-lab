from __future__ import annotations

import csv
from datetime import datetime, timedelta, timezone
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import ale_prospective_intake as intake


NOW = datetime(2026, 10, 15, 12, 0, tzinfo=timezone.utc)


def write_table(path: Path, name: str, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=intake.HEADERS[name])
        writer.writeheader()
        writer.writerows(rows)


def fixture(folder: Path, endpoints: bool = False) -> dict[str, Path]:
    candidates, cultures, measurements, decisions = [], [], [], []
    for population in ("new_pop_a", "new_pop_b"):
        for index, isolate in enumerate((f"{population}_clone_1", f"{population}_clone_2")):
            blind = f"blind_{population}_{index + 1}"
            candidates.append({
                "study_id": "ale_fumarate_validation_1", "cohort_role": "external_validation",
                "parent_strain_id": "Producer_F", "parent_genotype": "F parent",
                "population_id": population, "ale_condition_id": "ale_md20_v1", "isolate_id": isolate,
                "isolate_stock_id": f"stock_{isolate}", "genotype_qc_id": f"qc_{isolate}",
                "product": "fumarate", "species": "Saccharomyces cerevisiae", "blind_code": blind,
            })
            culture_id = f"culture_{isolate}_r1"
            cultures.append({
                "culture_id": culture_id, "study_id": "ale_fumarate_validation_1",
                "population_id": population, "isolate_id": isolate, "biological_replicate": "1",
                "production_condition_id": "md20_57h_v1", "medium": "MD20",
                "medium_batch_id": "medium_batch_1", "initial_glucose_g_L": "20",
                "initial_pH": "5", "temperature_C": "30", "agitation_rpm": "180",
                "vessel_type": "shake_flask", "vessel_volume_mL": "250", "working_volume_mL": "50",
                "inoculum_type": "washed_preculture", "inoculum_value": "0.1", "inoculum_unit": "OD600",
                "start_time_utc": "2026-10-13T12:00:00Z", "position": f"rack_{population}_{index}",
                "blind_code": blind, "protocol_id": "ale_fumarate_screen_v1",
            })
            for time_h, value, suffix in ((0, 0.1, "t0"), (18, 0.5 + index * 0.1, "t18")):
                sampled = datetime(2026, 10, 13, 12, tzinfo=timezone.utc) + timedelta(hours=time_h)
                measured = sampled + timedelta(minutes=5)
                available = measured + timedelta(minutes=5)
                measurements.append({
                    "measurement_id": f"m_{culture_id}_{suffix}", "assay_id": f"a_{culture_id}_{suffix}",
                    "study_id": "ale_fumarate_validation_1", "population_id": population,
                    "isolate_id": isolate, "culture_id": culture_id, "biological_replicate": "1",
                    "variable": "od600", "value": str(value), "unit": "OD600", "recorded_time_h": str(time_h),
                    "sampled_at_utc": sampled.isoformat().replace("+00:00", "Z"),
                    "measured_at_utc": measured.isoformat().replace("+00:00", "Z"),
                    "result_available_at_utc": available.isoformat().replace("+00:00", "Z"),
                    "assay_method": "OD600_reader", "assay_batch_id": "od_batch_1",
                    "sample_id": f"sample_{culture_id}_{suffix}", "raw_file_path": "", "raw_file_sha256": "",
                    "blind_code": blind, "technical_replicate": "1",
                })
            if endpoints:
                sampled = datetime(2026, 10, 13, 12, tzinfo=timezone.utc) + timedelta(hours=57)
                measured = sampled + timedelta(minutes=5)
                available = datetime(2026, 10, 16, 1, tzinfo=timezone.utc)
                measurements.append({
                    "measurement_id": f"m_{culture_id}_t57", "assay_id": f"a_{culture_id}_t57",
                    "study_id": "ale_fumarate_validation_1", "population_id": population,
                    "isolate_id": isolate, "culture_id": culture_id, "biological_replicate": "1",
                    "variable": "fumarate_titer", "value": str(0.2 + index * 0.1), "unit": "g/L",
                    "recorded_time_h": "57", "sampled_at_utc": sampled.isoformat().replace("+00:00", "Z"),
                    "measured_at_utc": measured.isoformat().replace("+00:00", "Z"),
                    "result_available_at_utc": available.isoformat().replace("+00:00", "Z"),
                    "assay_method": "UPLC_PDA_260nm_v1", "assay_batch_id": "fumarate_batch_1",
                    "sample_id": f"sample_{culture_id}_t57", "raw_file_path": "", "raw_file_sha256": "",
                    "blind_code": blind, "technical_replicate": "1",
                })
    decision_time = "2026-10-14T08:00:00Z"
    for candidate in candidates:
        population, isolate, blind = candidate["population_id"], candidate["isolate_id"], candidate["blind_code"]
        culture_id = f"culture_{isolate}_r1"
        ids = json.dumps([f"m_{culture_id}_t0", f"m_{culture_id}_t18"])
        for method in sorted(intake.METHODS):
            prediction_name = f"prediction_{method}_{population}.json"
            prediction_path = folder / prediction_name
            prediction_path.write_text(f"{method}:{population}", encoding="utf-8")
            top = isolate.endswith("clone_2")
            score = 2.0 if top else 1.0
            if method == "random_baseline":
                score = 0.9 if top else 0.1
            decisions.append({
                "decision_id": f"d_{method}_{isolate}", "study_id": "ale_fumarate_validation_1",
                "population_id": population, "isolate_id": isolate, "selection_method": method,
                "model_version": intake.MODEL_VERSION if method == "fixed_ridge" else f"{method}_v1",
                "score": str(score), "rank": "1" if top else "2", "selected": str(top).lower(),
                "assay_budget": "1", "decision_time_utc": decision_time, "input_cutoff_h": "18",
                "score_input_measurement_ids": "[]" if method == "random_baseline" else ids,
                "random_seed": "42" if method == "random_baseline" else "",
                "tie_break_rule": intake.TIE_BREAK_RULE,
                "endpoint_used_to_score": "false", "blind_code": blind,
                "prediction_file_path": prediction_name,
                "prediction_file_sha256": hashlib.sha256(f"{method}:{population}".encode()).hexdigest(),
            })
    paths = {name: folder / f"{name}.csv" for name in ("candidates", "cultures", "measurements", "decisions")}
    for name, rows in (("candidates", candidates), ("cultures", cultures),
                       ("measurements", measurements), ("decisions", decisions)):
        write_table(paths[name], name, rows)
    return paths


class AleProspectiveIntakeTests(unittest.TestCase):
    def test_freeze_then_evaluate_preserves_blind_and_timing_receipts(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            paths = fixture(folder)
            freeze = folder / "freeze.json"
            with patch.object(intake, "utc_now", return_value=NOW):
                frozen = intake.run("freeze", paths["candidates"], paths["cultures"], paths["measurements"], paths["decisions"], freeze)
            self.assertTrue(frozen["integrity_passed"])
            self.assertFalse(frozen["endpoint_labels_present"])
            self.assertFalse(frozen["strict_external_validation_allowed"])
            receipt = json.loads(freeze.read_text(encoding="utf-8"))
            self.assertEqual(receipt["record_type"], "pre_endpoint_selection_freeze")
            endpoint_paths = fixture(folder, endpoints=True)
            with patch.object(intake, "utc_now", return_value=datetime(2026, 10, 16, 2, tzinfo=timezone.utc)):
                report = intake.run("evaluate", endpoint_paths["candidates"], endpoint_paths["cultures"],
                                    endpoint_paths["measurements"], endpoint_paths["decisions"],
                                    folder / "evaluation.json", freeze)
            self.assertTrue(report["recorded_temporal_order_passed"])
            self.assertFalse(report["independent_blinding_proven"])
            self.assertEqual(report["candidate_count"], 4)
            self.assertEqual(report["population_count"], 2)
            self.assertEqual(report["endpoint_measurement_count"], 4)
            self.assertFalse(report["strict_external_validation_allowed"])
            self.assertFalse(report["wet_lab_gain_validated"])

    def test_freeze_rejects_endpoint_values_and_late_od(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            paths = fixture(folder, endpoints=True)
            with self.assertRaisesRegex(ValueError, "禁止导入任何富马酸终点"):
                intake.run("freeze", paths["candidates"], paths["cultures"], paths["measurements"], paths["decisions"], folder / "freeze.json")
            paths = fixture(folder)
            rows = intake.read_table(paths["measurements"], "measurements")
            late = dict(next(row for row in rows if row["recorded_time_h"] == "18"))
            late.update(measurement_id=late["measurement_id"] + "_late", assay_id=late["assay_id"] + "_late",
                        recorded_time_h="19", technical_replicate="2",
                        sampled_at_utc="2026-10-14T07:00:00Z", measured_at_utc="2026-10-14T07:05:00Z",
                        result_available_at_utc="2026-10-14T07:10:00Z", sample_id=late["sample_id"] + "_late")
            rows.append(late)
            write_table(paths["measurements"], "measurements", rows)
            with self.assertRaisesRegex(ValueError, "18 h 后"):
                intake.run("freeze", paths["candidates"], paths["cultures"], paths["measurements"], paths["decisions"], folder / "freeze2.json")

    def test_mismatched_culture_identity_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            paths = fixture(folder)
            rows = intake.read_table(paths["measurements"], "measurements")
            rows[0]["isolate_id"] = "some_other_clone"
            write_table(paths["measurements"], "measurements", rows)
            with self.assertRaisesRegex(ValueError, "身份字段冲突"):
                intake.run("freeze", paths["candidates"], paths["cultures"], paths["measurements"], paths["decisions"], folder / "freeze.json")

    def test_prediction_artifact_hash_is_rechecked(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            paths = fixture(folder)
            decisions = intake.read_table(paths["decisions"], "decisions")
            prediction = folder / decisions[0]["prediction_file_path"]
            prediction.write_text("changed after score export", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "预测文件 SHA256 不匹配"):
                intake.run("freeze", paths["candidates"], paths["cultures"], paths["measurements"], paths["decisions"], folder / "freeze.json")

    def test_nonmatching_production_recipe_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            paths = fixture(folder)
            cultures = intake.read_table(paths["cultures"], "cultures")
            cultures[0]["initial_glucose_g_L"] = "15"
            write_table(paths["cultures"], "cultures", cultures)
            with self.assertRaisesRegex(ValueError, "培养条件偏离已恢复"):
                intake.run("freeze", paths["candidates"], paths["cultures"], paths["measurements"], paths["decisions"], folder / "freeze.json")

    def test_decision_must_follow_all_feature_results(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            paths = fixture(folder)
            rows = intake.read_table(paths["decisions"], "decisions")
            rows[0]["decision_time_utc"] = "2026-10-14T06:02:00Z"
            write_table(paths["decisions"], "decisions", rows)
            with self.assertRaisesRegex(ValueError, "早于输入 OD 结果可用时刻"):
                intake.run("freeze", paths["candidates"], paths["cultures"], paths["measurements"], paths["decisions"], folder / "freeze.json")

    def test_target_must_follow_freeze_and_be_available_at_evaluation(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            paths = fixture(folder)
            freeze = folder / "freeze.json"
            with patch.object(intake, "utc_now", return_value=datetime(2026, 10, 16, 2, tzinfo=timezone.utc)):
                intake.run("freeze", paths["candidates"], paths["cultures"], paths["measurements"], paths["decisions"], freeze)
            endpoints = fixture(folder, endpoints=True)
            rows = intake.read_table(endpoints["measurements"], "measurements")
            for row in rows:
                if row["variable"] == "fumarate_titer":
                    row["result_available_at_utc"] = "2026-10-16T01:00:00Z"
            write_table(endpoints["measurements"], "measurements", rows)
            with patch.object(intake, "utc_now", return_value=datetime(2026, 10, 16, 2, tzinfo=timezone.utc)):
                with self.assertRaisesRegex(ValueError, "冻结回执生成前已可用"):
                    intake.run("evaluate", endpoints["candidates"], endpoints["cultures"], endpoints["measurements"],
                               endpoints["decisions"], folder / "evaluation.json", freeze)

    def test_target_endpoint_time_is_fixed_at_57_hours(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            paths = fixture(folder)
            freeze = folder / "freeze.json"
            with patch.object(intake, "utc_now", return_value=NOW):
                intake.run("freeze", paths["candidates"], paths["cultures"], paths["measurements"], paths["decisions"], freeze)
            endpoints = fixture(folder, endpoints=True)
            rows = intake.read_table(endpoints["measurements"], "measurements")
            for row in rows:
                if row["variable"] == "fumarate_titer":
                    row["recorded_time_h"] = "56"
                    row["sampled_at_utc"] = "2026-10-15T20:00:00Z"
                    row["measured_at_utc"] = "2026-10-15T20:05:00Z"
            write_table(endpoints["measurements"], "measurements", rows)
            with patch.object(intake, "utc_now", return_value=datetime(2026, 10, 16, 2, tzinfo=timezone.utc)):
                with self.assertRaisesRegex(ValueError, "必须有来源定义的 57 h"):
                    intake.run("evaluate", endpoints["candidates"], endpoints["cultures"], endpoints["measurements"],
                               endpoints["decisions"], folder / "evaluation.json", freeze)

    def test_decision_file_and_early_data_are_bound_to_freeze(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            paths = fixture(folder)
            freeze = folder / "freeze.json"
            with patch.object(intake, "utc_now", return_value=NOW):
                intake.run("freeze", paths["candidates"], paths["cultures"], paths["measurements"], paths["decisions"], freeze)
            endpoints = fixture(folder, endpoints=True)
            rows = intake.read_table(endpoints["measurements"], "measurements")
            next(row for row in rows if row["variable"] == "od600" and row["recorded_time_h"] == "18")["value"] = "0.55"
            write_table(endpoints["measurements"], "measurements", rows)
            with patch.object(intake, "utc_now", return_value=datetime(2026, 10, 16, 2, tzinfo=timezone.utc)):
                with self.assertRaisesRegex(ValueError, "OD 输入与冻结回执不一致"):
                    intake.run("evaluate", endpoints["candidates"], endpoints["cultures"], endpoints["measurements"],
                               endpoints["decisions"], folder / "evaluation.json", freeze)

    def test_budget_and_selected_flags_must_agree(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            paths = fixture(folder)
            rows = intake.read_table(paths["decisions"], "decisions")
            rows[0]["selected"] = "false" if rows[0]["selected"] == "true" else "true"
            write_table(paths["decisions"], "decisions", rows)
            with self.assertRaisesRegex(ValueError, "selected 数量与 assay_budget"):
                intake.run("freeze", paths["candidates"], paths["cultures"], paths["measurements"], paths["decisions"], folder / "freeze.json")

    def test_existing_receipt_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as temp:
            folder = Path(temp)
            paths = fixture(folder)
            output = folder / "freeze.json"
            output.write_text("protected", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "必须是新文件"):
                intake.run("freeze", paths["candidates"], paths["cultures"], paths["measurements"], paths["decisions"], output)
            self.assertEqual(output.read_text(encoding="utf-8"), "protected")


if __name__ == "__main__":
    unittest.main()
