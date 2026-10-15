from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from datetime import datetime, timezone
from pathlib import Path

import openpyxl

ROOT = Path(__file__).resolve().parent
PROTOCOL = ROOT / "research/ale_external_screen_protocol_20261015.json"
DEFAULT_OUTPUT = ROOT / "research/results/ale_external_screen_v1"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def number(value) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"要求有限原始数值，不能补齐缺失或使用未计算公式：{value!r}")
    return float(value)


def select_weights(scores: list[float], tolerance: float) -> list[float]:
    if not scores or any(not math.isfinite(score) for score in scores):
        raise ValueError("选株分数缺失或非有限")
    maximum = max(scores)
    selected = [i for i, score in enumerate(scores) if abs(score - maximum) <= tolerance]
    return [1 / len(selected) if i in selected else 0.0 for i in range(len(scores))]


def extract_panel(book, spec: dict, minimum_replicates: int) -> list[dict]:
    early_sheet, target_sheet = book[spec["early_sheet"]], book[spec["target_sheet"]]
    expected_time = {"early": spec["cutoff_h"], "target": spec["endpoint_h"]}
    if spec["initial_row"] is not None:
        expected_time["initial"] = spec["initial_h"]
    for role, value in expected_time.items():
        sheet = target_sheet if role == "target" else early_sheet
        raw = sheet[spec["time_cells"][role]].value
        observed = number(raw) if isinstance(raw, (int, float)) else float(str(raw).split()[0])
        if observed != value:
            raise ValueError("原始时间与固定窗口不符")
    if spec["cutoff_h"] >= spec["endpoint_h"]:
        raise ValueError("早期时钟必须严格早于终点")
    ids = [candidate["isolate_id"] for candidate in spec["candidates"]]
    if len(ids) != len(set(ids)) or len(ids) < 2:
        raise ValueError("选株需要至少两个不同的候选菌株")
    records = []
    for candidate in spec["candidates"]:
        columns = candidate["columns"]
        if len(columns) != len(set(columns)) or len(columns) < minimum_replicates:
            raise ValueError("重复列无效或重复数不足")
        header = spec.get("header_row", 1)
        if any(sheet[f"{columns[0]}{header}"].value != candidate["isolate_id"] for sheet in (early_sheet, target_sheet)):
            raise ValueError("早期与终点菌株身份不一致")
        early_cells = [f"{column}{spec['early_row']}" for column in columns]
        target_cells = [f"{column}{spec['target_row']}" for column in columns]
        early = [number(early_sheet[cell].value) for cell in early_cells]
        target = [number(target_sheet[cell].value) for cell in target_cells]
        initial_cells, initial = [], []
        if spec["initial_row"] is not None:
            initial_cells = [f"{column}{spec['initial_row']}" for column in columns]
            initial = [number(early_sheet[cell].value) for cell in initial_cells]
            if any(value <= 0 for value in initial + early):
                raise ValueError("OD 初值和早期值必须为正")
        if any(value < 0 for value in target):
            raise ValueError("终点数值为负；需核查来源，不截断修分")
        records.append({
            "panel_id": spec["id"], "isolate_id": candidate["isolate_id"],
            "source_path": spec["source"], "early_sheet": spec["early_sheet"], "target_sheet": spec["target_sheet"],
            "cutoff_h": spec["cutoff_h"], "endpoint_h": spec["endpoint_h"],
            "early_cells": early_cells, "target_cells": target_cells, "initial_cells": initial_cells,
            "early_replicates": early, "target_replicates": target, "initial_replicates": initial,
            "early_mean": statistics.fmean(early), "target_mean": statistics.fmean(target),
            "initial_mean": statistics.fmean(initial) if initial else None,
            "replicate_count": len(columns), "target_unit": spec["target_unit"],
        })
    return records


def score_panel(records: list[dict], spec: dict, tolerance: float) -> tuple[list, list]:
    decisions, metrics = [], []
    targets = [row["target_mean"] for row in records]
    best = max(targets)
    if best <= 0:
        raise ValueError("全部终点为零，无法计算归一化选株指标")
    for method in spec["methods"]:
        if method == "uniform_random_expectation":
            scores = [0.0] * len(records)
        elif method in ("early_density", "early_product"):
            scores = [row["early_mean"] for row in records]
        elif method == "early_density_fold":
            scores = [row["early_mean"] / row["initial_mean"] for row in records]
        else:
            raise ValueError("不支持的固定筛选规则")
        weights = select_weights(scores, tolerance)
        selected = math.fsum(weight * target for weight, target in zip(weights, targets))
        for row, value, weight in zip(records, scores, weights):
            decisions.append({"panel_id": spec["id"], "method": method, "isolate_id": row["isolate_id"], "score": value, "selection_probability": weight, "test_endpoint_used_to_score": False})
        metrics.append({
            "panel_id": spec["id"], "method": method, "candidate_strain_count": len(records),
            "expected_endpoint_value": selected, "target_unit": spec["target_unit"],
            "fraction_of_panel_best": selected / best, "absolute_regret": best - selected,
            "top1_hit_probability": math.fsum(weight for weight, target in zip(weights, targets) if abs(target - best) <= tolerance),
            "early_assay": "product assay on all candidates" if method == "early_product" else "OD on all candidates" if method.startswith("early_density") else "none",
            "selection_budget_strains": 1, "actual_endpoint_assays_saved_verified": False,
        })
    return decisions, metrics


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(value, ensure_ascii=False) if isinstance(value, list) else value for key, value in row.items()})


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", default=str(DEFAULT_OUTPUT))
    args = parser.parse_args()
    protocol = json.loads(PROTOCOL.read_text(encoding="utf-8"))
    for rel, digest in protocol["input_sha256"].items():
        if sha(ROOT / rel) != digest:
            raise ValueError(f"原始输入变化，拒绝覆盖旧结果：{rel}")
    out = Path(args.output)
    if out.exists() and any(out.iterdir()):
        raise ValueError("结果目录非空；使用新目录保留既有回放")
    out.mkdir(parents=True, exist_ok=True)
    all_records, all_decisions, all_metrics = [], [], []
    for spec in protocol["panels"]:
        book = openpyxl.load_workbook(ROOT / spec["source"], data_only=False, read_only=True)
        try:
            records = extract_panel(book, spec, protocol["minimum_replicates_per_candidate"])
        finally:
            book.close()
        decisions, metrics = score_panel(records, spec, protocol["tie_absolute_tolerance"])
        all_records.extend(records)
        all_decisions.extend(decisions)
        all_metrics.extend(metrics)
    write_csv(out / "candidates.csv", all_records)
    write_csv(out / "decisions.csv", all_decisions)
    write_csv(out / "metrics.csv", all_metrics)
    report = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(), "protocol_sha256": sha(PROTOCOL),
        "input_sha256": protocol["input_sha256"],
        "scope": {key: protocol[key] for key in ("source_values_and_targets_seen_before_lock", "new_blind_validation", "original_18h_fumarate_ridge_external_validation", "publication_superiority_claim", "same_strain_omics_process_closed_loop", "model_training", "prospective_assay_savings_verified")},
        "panel_specs": protocol["panels"], "excluded_panels": protocol["excluded_panels"],
        "candidates": all_records, "decisions": all_decisions, "metrics": all_metrics,
        "artifact_sha256": {name: sha(out / name) for name in ("candidates.csv", "decisions.csv", "metrics.csv")},
        "implementation_sha256": {name: sha(ROOT / name) for name in ("external_screen_validation.py", "verify_external_screen_validation.py")},
    }
    (out / "report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    for metric in all_metrics:
        print(f"{metric['panel_id']} / {metric['method']}: 最优比例 {metric['fraction_of_panel_best']:.4f}; top1 {metric['top1_hit_probability']:.4f}")
    print("两个面板独立计分；原18h富马酸岭回归未获得外部验证；无盲测或实测节省收益")
    print(f"输出：{out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
