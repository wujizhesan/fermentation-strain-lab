from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import posixpath
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path
from urllib.parse import parse_qs, urlsplit

ROOT = Path(__file__).resolve().parent
PROTOCOL = ROOT / "research/ale_external_screen_protocol_20261015.json"
DEFAULT_RESULT = ROOT / "research/results/ale_external_screen_v1/report.json"
NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
RID = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sheet_cells(path: Path, name: str) -> dict:
    with zipfile.ZipFile(path) as archive:
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relations = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {r.get("Id"): r.get("Target") for r in relations}
        reference = next(s for s in workbook.findall("m:sheets/m:sheet", NS) if s.get("name") == name)
        target = targets[reference.get(RID)]
        location = target.lstrip("/") if target.startswith("/") else posixpath.normpath("xl/" + target)
        shared = []
        if "xl/sharedStrings.xml" in archive.namelist():
            shared = ["".join(node.itertext()) for node in ET.fromstring(archive.read("xl/sharedStrings.xml")).findall("m:si", NS)]
        result = {}
        for cell in ET.fromstring(archive.read(location)).findall(".//m:sheetData/m:row/m:c", NS):
            coordinate = cell.get("r")
            if cell.find("m:f", NS) is not None:
                result[coordinate] = "formula"
                continue
            value = cell.findtext("m:v", default=None, namespaces=NS)
            kind = cell.get("t")
            if kind == "s":
                value = shared[int(value)]
            elif kind == "inlineStr":
                value = "".join(node.text or "" for node in cell.findall(".//m:t", NS))
            elif kind == "b":
                value = value == "1"
            elif value is not None:
                value = float(value)
            result[coordinate] = value
        return result


def equivalent(actual, expected) -> bool:
    if isinstance(expected, dict):
        return isinstance(actual, dict) and actual.keys() == expected.keys() and all(equivalent(actual[k], expected[k]) for k in expected)
    if isinstance(expected, list):
        return isinstance(actual, list) and len(actual) == len(expected) and all(equivalent(a, b) for a, b in zip(actual, expected))
    if isinstance(expected, bool) or expected is None:
        return actual is expected
    if isinstance(expected, (int, float)):
        return type(actual) in (int, float) and math.isfinite(actual) and math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-12)
    return actual == expected


def _numeric(cells: dict, coordinate: str) -> float:
    value = cells.get(coordinate)
    if type(value) not in (float, int) or not math.isfinite(value):
        raise ValueError(f"原始数值缺失或非法：{coordinate}")
    return float(value)


def verify(result: Path) -> dict[str, bool]:
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    report = json.loads(result.read_text(encoding="utf-8"))
    manifest = json.loads((ROOT / protocol["source_manifest"]).read_text(encoding="utf-8"))
    source_ok = True
    for row in manifest["files"]:
        metadata = json.loads((ROOT / manifest["metadata"][row["pmcid"]]["path"]).read_text(encoding="utf-8"))
        blob = (ROOT / row["path"]).read_bytes()
        resource = next(url for url in metadata["media_urls"] if urlsplit(url).path.rsplit("/", 1)[-1] == Path(row["path"]).name)
        source_ok &= row["downloaded"] is True and row["md5_verified"] is True
        source_ok &= hashlib.md5(blob).hexdigest() == row["provider_md5"] == parse_qs(urlsplit(resource).query)["md5"][0]
        source_ok &= sha(ROOT / row["path"]) == row["sha256"] and len(blob) == row["bytes"]
        source_ok &= row["doi"] == metadata["doi"] and row["license"] == metadata["license_code"]
    records, decisions, metrics = [], [], []
    identity_ok, time_ok = True, True
    for spec in protocol["panels"]:
        early, target = sheet_cells(ROOT / spec["source"], spec["early_sheet"]), sheet_cells(ROOT / spec["source"], spec["target_sheet"])
        identity_ok &= len(spec["candidates"]) >= 2 and len({r["isolate_id"] for r in spec["candidates"]}) == len(spec["candidates"])
        time_ok &= spec["cutoff_h"] < spec["endpoint_h"]
        for role, value in (("early", spec["cutoff_h"]), ("target", spec["endpoint_h"])):
            cells = early if role == "early" else target
            observed = cells[spec["time_cells"][role]]
            time_ok &= float(str(observed).split()[0]) == value
        if spec["initial_row"] is not None:
            time_ok &= early[spec["time_cells"]["initial"]] == spec["initial_h"] == 0
        panel_records = []
        for candidate in spec["candidates"]:
            columns = candidate["columns"]
            identity_ok &= len(set(columns)) == len(columns) and len(columns) >= protocol["minimum_replicates_per_candidate"]
            for cells in (early, target):
                identity_ok &= cells.get(columns[0] + str(spec.get("header_row", 1))) == candidate["isolate_id"]
            coordinates = {"early": [c + str(spec["early_row"]) for c in columns], "target": [c + str(spec["target_row"]) for c in columns], "initial": [c + str(spec["initial_row"]) for c in columns] if spec["initial_row"] is not None else []}
            values = {role: [_numeric(target if role == "target" else early, coordinate) for coordinate in refs] for role, refs in coordinates.items()}
            mean = lambda values: math.fsum(values) / len(values)
            if values["initial"] and any(value <= 0 for value in values["initial"] + values["early"]):
                raise ValueError("OD 必须为正")
            if any(value < 0 for value in values["target"]):
                raise ValueError("终点为负，禁止静默截断")
            row = {"panel_id": spec["id"], "isolate_id": candidate["isolate_id"], "source_path": spec["source"], "early_sheet": spec["early_sheet"], "target_sheet": spec["target_sheet"], "cutoff_h": spec["cutoff_h"], "endpoint_h": spec["endpoint_h"], "early_cells": coordinates["early"], "target_cells": coordinates["target"], "initial_cells": coordinates["initial"], "early_replicates": values["early"], "target_replicates": values["target"], "initial_replicates": values["initial"], "early_mean": mean(values["early"]), "target_mean": mean(values["target"]), "initial_mean": mean(values["initial"]) if values["initial"] else None, "replicate_count": len(columns), "target_unit": spec["target_unit"]}
            panel_records.append(row)
        records.extend(panel_records)
        best = max(r["target_mean"] for r in panel_records)
        if best <= 0:
            raise ValueError("全零终点不能计算最优比例")
        for method in spec["methods"]:
            scores = [0.0 if method == "uniform_random_expectation" else r["early_mean"] / r["initial_mean"] if method == "early_density_fold" else r["early_mean"] for r in panel_records]
            if method not in ("uniform_random_expectation", "early_density", "early_density_fold", "early_product"):
                raise ValueError("未声明的筛选规则")
            chosen = [i for i, score in enumerate(scores) if math.isclose(score, max(scores), rel_tol=0, abs_tol=protocol["tie_absolute_tolerance"])]
            probabilities = [1 / len(chosen) if i in chosen else 0.0 for i in range(len(scores))]
            selected = math.fsum(r["target_mean"] * probabilities[i] for i, r in enumerate(panel_records))
            for row, score, probability in zip(panel_records, scores, probabilities):
                decisions.append({"panel_id": spec["id"], "method": method, "isolate_id": row["isolate_id"], "score": score, "selection_probability": probability, "test_endpoint_used_to_score": False})
            metrics.append({"panel_id": spec["id"], "method": method, "candidate_strain_count": len(panel_records), "expected_endpoint_value": selected, "target_unit": spec["target_unit"], "fraction_of_panel_best": selected / best, "absolute_regret": best - selected, "top1_hit_probability": math.fsum(probabilities[i] for i, r in enumerate(panel_records) if math.isclose(r["target_mean"], best, rel_tol=0, abs_tol=protocol["tie_absolute_tolerance"])), "early_assay": "product assay on all candidates" if method == "early_product" else "OD on all candidates" if method.startswith("early_density") else "none", "selection_budget_strains": 1, "actual_endpoint_assays_saved_verified": False})
    scope_keys = ("source_values_and_targets_seen_before_lock", "new_blind_validation", "original_18h_fumarate_ridge_external_validation", "publication_superiority_claim", "same_strain_omics_process_closed_loop", "model_training", "prospective_assay_savings_verified")
    scope = {key: protocol[key] for key in scope_keys}
    file_hashes = {name: sha(result.parent / name) for name in ("candidates.csv", "decisions.csv", "metrics.csv")}
    csv_ok = True
    for name, expected in (("candidates.csv", records), ("decisions.csv", decisions), ("metrics.csv", metrics)):
        with (result.parent / name).open(encoding="utf-8", newline="") as stream:
            actual = list(csv.DictReader(stream))
        projection = [{key: json.dumps(value, ensure_ascii=False) if isinstance(value, list) else "" if value is None else str(value) for key, value in row.items()} for row in expected]
        csv_ok &= actual == projection
    ogataea = sheet_cells(ROOT / "data/raw/ale_external_sources_20261015/PMC12663439/41467_2025_65793_MOESM14_ESM.xlsx", "Fig. S5")
    yarrowia = sheet_cells(ROOT / protocol["panels"][0]["source"], "Supplementary Fig. 3")
    return {
        "协议 SHA": report.get("protocol_sha256") == sha(PROTOCOL),
        "原始输入 SHA": report.get("input_sha256") == protocol["input_sha256"] and all(sha(ROOT / path) == digest for path, digest in protocol["input_sha256"].items()),
        "公开元数据 MD5 与来源身份": source_ok,
        "候选身份与重复数": identity_ok,
        "早期与终点时钟": time_ok,
        "逐候选原始单元格独立复算": equivalent(report.get("candidates"), records),
        "全部选择分数与概率独立复算": equivalent(report.get("decisions"), decisions),
        "全部面板指标独立复算": equivalent(report.get("metrics"), metrics),
        "CSV 内容绑定原始单元格": csv_ok,
        "CSV SHA": report.get("artifact_sha256") == file_hashes,
        "实现 SHA": report.get("implementation_sha256") == {name: sha(ROOT / name) for name in ("external_screen_validation.py", "verify_external_screen_validation.py")},
        "面板与排除项绑定协议": report.get("panel_specs") == protocol["panels"] and report.get("excluded_panels") == protocol["excluded_panels"],
        "保留复现范围": report.get("scope") == scope and scope["source_values_and_targets_seen_before_lock"] is True and all(scope[key] is False for key in scope_keys[1:]),
        "未伪造化验节省与时间可用性": protocol["assay_return_times_verified"] is False and protocol["observation_available_at_sampling_time_is_assumed"] is True and protocol["methods_can_differ_in_early_assay_cost"] is True,
        "甲醇酵母缺失重复有原始依据": all(ogataea.get(column + str(row)) is None for row in (8, 9, 24, 25) for column in ("B", "C", "D", "E", "F")),
        "油脂酵母匿名终点库未拼接": [yarrowia.get("A" + str(i)) for i in range(2, 22)] == list(range(1, 21)) and yarrowia.get("C1") == "SA titer (g/L)",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result", default=str(DEFAULT_RESULT))
    parser.add_argument("--receipt")
    args = parser.parse_args()
    result = Path(args.result)
    try:
        checks = verify(result)
    except (OSError, ValueError, KeyError, TypeError, StopIteration, zipfile.BadZipFile) as exc:
        checks = {"输入与证据有效": False}
        print(f"[FAIL] {exc}")
    passed = all(checks.values())
    receipt = {"passed": passed, "report_sha256": sha(result), "protocol_sha256": sha(PROTOCOL), "checks": checks}
    path = Path(args.receipt) if args.receipt else result.with_name("verification.json")
    path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    for name, ok in checks.items():
        print(f"[{'PASS' if ok else 'FAIL'}] {name}")
    print(f"{sum(checks.values())}/{len(checks)} PASS")
    return int(not passed)


if __name__ == "__main__":
    raise SystemExit(main())
