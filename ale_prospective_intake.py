from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parent
MODEL_PATH = ROOT / "research/results/ale_external_model_v1/fitted_model.json"
STATUS_PATH = ROOT / "research/ale_prospective_validation_status_20261015.json"
SOURCE_ALIGNMENT_PATH = ROOT / "research/ale_source_process_alignment_20261015.json"
EARLY_CUTOFF_H = 18.0
ENDPOINT_H = 57.0
TIME_TOLERANCE_H = 0.25
TARGET_VARIABLE = "fumarate_titer"
MODEL_VERSION = "locked_fixed_ridge_v1"
TIE_BREAK_RULE = "score_descending_then_isolate_id_ascending"
METHODS = {"fixed_ridge", "early_od_mean", "early_od_fold", "random_baseline"}

HEADERS = {
    "candidates": (
        "study_id", "cohort_role", "parent_strain_id", "parent_genotype", "population_id",
        "ale_condition_id", "isolate_id", "isolate_stock_id", "genotype_qc_id", "product",
        "species", "blind_code",
    ),
    "cultures": (
        "culture_id", "study_id", "population_id", "isolate_id", "biological_replicate",
        "production_condition_id", "medium", "medium_batch_id", "initial_glucose_g_L",
        "initial_pH", "temperature_C", "agitation_rpm", "vessel_type", "vessel_volume_mL",
        "working_volume_mL", "inoculum_type", "inoculum_value", "inoculum_unit",
        "start_time_utc", "position", "blind_code", "protocol_id",
    ),
    "measurements": (
        "measurement_id", "assay_id", "study_id", "population_id", "isolate_id", "culture_id",
        "biological_replicate", "variable", "value", "unit", "recorded_time_h", "sampled_at_utc",
        "measured_at_utc", "result_available_at_utc", "assay_method", "assay_batch_id", "sample_id",
        "raw_file_path", "raw_file_sha256", "blind_code", "technical_replicate",
    ),
    "decisions": (
        "decision_id", "study_id", "population_id", "isolate_id", "selection_method",
        "model_version", "score", "rank", "selected", "assay_budget", "decision_time_utc",
        "input_cutoff_h", "score_input_measurement_ids", "random_seed", "tie_break_rule", "endpoint_used_to_score",
        "blind_code", "prediction_file_path", "prediction_file_sha256",
    ),
    "omics": (
        "omics_sample_id", "study_id", "population_id", "isolate_id", "culture_id",
        "biological_replicate", "modality", "collection_time_h", "collection_time_utc",
        "condition_id", "aliquot_id", "raw_file_path", "raw_file_sha256", "processing_protocol_id",
    ),
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows_sha256(rows: list[dict]) -> str:
    payload = json.dumps(sorted(rows, key=lambda row: row["measurement_id"]), ensure_ascii=False,
                         sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def utc_now() -> datetime:
    return datetime.now(timezone.utc)


def parse_time(value: str, label: str) -> datetime:
    try:
        result = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as exc:
        raise ValueError(f"{label} 必须是带时区的 ISO 8601 时间：{value!r}") from exc
    if result.tzinfo is None or result.utcoffset() is None:
        raise ValueError(f"{label} 必须包含时区：{value!r}")
    return result.astimezone(timezone.utc)


def number(value: str, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} 必须是数值：{value!r}") from exc
    if not math.isfinite(result):
        raise ValueError(f"{label} 必须是有限数值：{value!r}")
    return result


def positive_int(value: str, label: str) -> int:
    try:
        result = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} 必须是正整数：{value!r}") from exc
    if result < 1 or str(result) != value.strip():
        raise ValueError(f"{label} 必须是正整数：{value!r}")
    return result


def boolean(value: str, label: str) -> bool:
    normalized = value.strip().lower()
    if normalized not in {"true", "false"}:
        raise ValueError(f"{label} 只能是 true 或 false：{value!r}")
    return normalized == "true"


def read_table(path: Path, name: str) -> list[dict]:
    if not path.is_file():
        raise ValueError(f"找不到{name} CSV：{path}")
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != HEADERS[name]:
            raise ValueError(f"{name} 表头不匹配；请使用 research/templates 中对应模板")
        result = []
        for line, row in enumerate(reader, start=2):
            if row is None or all(not (value or "").strip() for value in row.values()):
                continue
            if None in row:
                raise ValueError(f"{name} 第 {line} 行列数与表头不一致")
            result.append({key: (value or "").strip() for key, value in row.items()})
    if not result:
        raise ValueError(f"{name} 表没有数据行")
    for line, row in enumerate(result, start=2):
        missing = [key for key, value in row.items() if not value and key not in {"raw_file_path", "raw_file_sha256", "random_seed"}]
        if missing:
            raise ValueError(f"{name} 第 {line} 行缺少必填字段：{', '.join(missing)}")
    return result


def validate_model() -> tuple[dict, str]:
    if not MODEL_PATH.is_file() or not STATUS_PATH.is_file():
        raise ValueError("锁定模型或验证状态文件缺失")
    model_sha = sha256(MODEL_PATH)
    model = json.loads(MODEL_PATH.read_text(encoding="utf-8"))
    status = json.loads(STATUS_PATH.read_text(encoding="utf-8"))
    if status.get("external_model_freeze", {}).get("model_sha256") != model_sha:
        raise ValueError("冻结模型 SHA 与前瞻验证状态不一致")
    if (model.get("endpoint", {}).get("product") != "fumarate"
            or model.get("endpoint", {}).get("unit") != "g/L"
            or model.get("endpoint", {}).get("time_h") != ENDPOINT_H
            or model.get("model", {}).get("features") != ["initial_od600", "log_early_od_fold"]):
        raise ValueError("冻结模型的目标或特征定义与入口协议不一致")
    if status.get("development_task", {}).get("early_cutoff_h") != EARLY_CUTOFF_H:
        raise ValueError("验证状态中的早期截止时间与入口协议不一致")
    return model, model_sha


def validate_omics(rows: list[dict], cultures: dict, candidates: dict) -> None:
    seen = set()
    for row in rows:
        sample_id = row["omics_sample_id"]
        aliquot_id = row["aliquot_id"]
        if sample_id in seen or aliquot_id in seen:
            raise ValueError("omics_sample_id 和 aliquot_id 必须唯一")
        seen.update((sample_id, aliquot_id))
        key = (row["study_id"], row["population_id"], row["isolate_id"])
        if key not in candidates:
            raise ValueError(f"组学样本找不到候选菌株：{key}")
        culture = cultures.get(row["culture_id"])
        if culture is None or (culture["study_id"], culture["population_id"], culture["isolate_id"]) != key:
            raise ValueError(f"组学样本 culture_id 关系冲突：{sample_id}")
        if positive_int(row["biological_replicate"], "omics biological_replicate") != int(culture["biological_replicate"]):
            raise ValueError(f"组学样本重复号与培养物不一致：{sample_id}")
        time_h = number(row["collection_time_h"], "collection_time_h")
        if time_h < 0:
            raise ValueError("collection_time_h 不能为负")
        collected = parse_time(row["collection_time_utc"], "collection_time_utc")
        started = parse_time(culture["start_time_utc"], "start_time_utc")
        elapsed = (collected - started).total_seconds() / 3600
        if elapsed < 0 or abs(elapsed - time_h) > TIME_TOLERANCE_H:
            raise ValueError(f"组学采样时刻与培养起始时间不匹配：{sample_id}")
        if row["condition_id"] != culture["production_condition_id"]:
            raise ValueError(f"组学 condition_id 与培养条件不一致：{sample_id}")
        validate_optional_sha(row, sample_id)


def validate_optional_sha(row: dict, label: str) -> None:
    digest, path = row.get("raw_file_sha256", ""), row.get("raw_file_path", "")
    if bool(digest) != bool(path):
        raise ValueError(f"{label} 的 raw_file_path 与 SHA256 必须同时填写或同时留空")
    if digest and not re.fullmatch(r"[0-9a-fA-F]{64}", digest):
        raise ValueError(f"{label} 的 raw_file_sha256 格式错误")


def validate_raw_files(rows: list[dict], base: Path, label: str) -> tuple[int, int]:
    verified = 0
    missing = 0
    for row in rows:
        path_value = row.get("raw_file_path", "")
        if not path_value:
            missing += 1
            continue
        path = Path(path_value)
        if not path.is_absolute():
            path = base / path
        if not path.is_file():
            raise ValueError(f"{label} 原始文件不存在：{path_value}")
        if sha256(path).lower() != row["raw_file_sha256"].lower():
            raise ValueError(f"{label} 原始文件 SHA256 不匹配：{path_value}")
        verified += 1
    return verified, missing


def validate_tables(candidate_rows: list[dict], culture_rows: list[dict], measurement_rows: list[dict],
                    decision_rows: list[dict], stage: str, freeze_receipt: dict | None = None,
                    omics_rows: list[dict] | None = None, now: datetime | None = None) -> dict:
    model, model_sha = validate_model()
    now = now or utc_now()
    if now.tzinfo is None or now.utcoffset() is None:
        raise ValueError("验证时钟必须带时区")
    now = now.astimezone(timezone.utc)

    candidates = {}
    by_population = {}
    group_parent = {}
    blind_codes = set()
    for row in candidate_rows:
        key = (row["study_id"], row["population_id"], row["isolate_id"])
        if key in candidates:
            raise ValueError(f"候选菌株键重复：{key}")
        if row["cohort_role"] != "external_validation":
            raise ValueError(f"本入口只接收 external_validation 队列：{key}")
        if row["product"].lower() != "fumarate" or row["species"] != "Saccharomyces cerevisiae":
            raise ValueError(f"候选菌株与冻结的酿酒酵母富马酸任务不匹配：{key}")
        if row["blind_code"] in blind_codes:
            raise ValueError(f"blind_code 必须唯一：{row['blind_code']}")
        blind_codes.add(row["blind_code"])
        if row["isolate_id"] in model["model"]["training_isolate_ids"] or row["population_id"] in model["model"]["training_lineage_groups"]:
            raise ValueError(f"候选菌株或群体与模型开发数据重叠：{key}")
        candidates[key] = row
        by_population.setdefault((row["study_id"], row["population_id"]), set()).add(row["isolate_id"])
        pop_key = (row["study_id"], row["population_id"])
        parent = (row["parent_strain_id"], row["parent_genotype"], row["ale_condition_id"])
        if pop_key in group_parent and group_parent[pop_key] != parent:
            raise ValueError(f"同一 ALE 群体的亲本或演化条件冲突：{pop_key}")
        group_parent[pop_key] = parent
    if len({row["study_id"] for row in candidate_rows}) != 1:
        raise ValueError("一次验证入口只允许一个 study_id")
    for population, isolates in by_population.items():
        if len(isolates) < 2:
            raise ValueError(f"每个独立 ALE 群体至少要有两个候选菌株：{population}")

    cultures = {}
    cultures_by_candidate = {}
    condition_signatures = {}
    replicate_keys = set()
    for row in culture_rows:
        culture_id = row["culture_id"]
        if culture_id in cultures:
            raise ValueError(f"culture_id 重复：{culture_id}")
        key = (row["study_id"], row["population_id"], row["isolate_id"])
        if key not in candidates:
            raise ValueError(f"培养物找不到候选菌株：{culture_id}")
        candidate = candidates[key]
        if row["blind_code"] != candidate["blind_code"]:
            raise ValueError(f"培养物盲码与候选菌株不一致：{culture_id}")
        rep = positive_int(row["biological_replicate"], "biological_replicate")
        condition_rep_key = (*key, row["production_condition_id"], rep)
        if condition_rep_key in replicate_keys:
            raise ValueError(f"同一菌株生产条件下生物学重复号重复：{condition_rep_key}")
        replicate_keys.add(condition_rep_key)
        start = parse_time(row["start_time_utc"], "start_time_utc")
        inoculum = number(row["inoculum_value"], "inoculum_value")
        if inoculum <= 0 or row["inoculum_unit"] != "OD600":
            raise ValueError(f"model 特征要求以 OD600 记录正值接种量：{culture_id}")
        numeric = {
            "initial_glucose_g_L": (number(row["initial_glucose_g_L"], "initial_glucose_g_L"), 0, None),
            "initial_pH": (number(row["initial_pH"], "initial_pH"), 0, 14),
            "temperature_C": (number(row["temperature_C"], "temperature_C"), 0, None),
            "agitation_rpm": (number(row["agitation_rpm"], "agitation_rpm"), 0, None),
            "vessel_volume_mL": (number(row["vessel_volume_mL"], "vessel_volume_mL"), 0, None),
            "working_volume_mL": (number(row["working_volume_mL"], "working_volume_mL"), 0, None),
        }
        if any(value <= lower or (upper is not None and value > upper) for value, lower, upper in numeric.values()):
            raise ValueError(f"培养条件数值超出允许范围：{culture_id}")
        if numeric["working_volume_mL"][0] >= numeric["vessel_volume_mL"][0]:
            raise ValueError(f"working_volume_mL 必须小于 vessel_volume_mL：{culture_id}")
        if (row["medium"] != "MD20"
                or not math.isclose(numeric["initial_glucose_g_L"][0], 20.0, abs_tol=1e-9)
                or not math.isclose(numeric["initial_pH"][0], 5.0, abs_tol=1e-9)
                or not math.isclose(numeric["temperature_C"][0], 30.0, abs_tol=1e-9)
                or not math.isclose(numeric["agitation_rpm"][0], 180.0, abs_tol=1e-9)
                or not math.isclose(inoculum, 0.1, rel_tol=0.05, abs_tol=1e-8)):
            raise ValueError(f"培养条件偏离已恢复的 MD20/20 g/L 葡萄糖/pH 5/30°C/180 rpm/初始 OD600 0.1 任务：{culture_id}")
        signature = (row["production_condition_id"], row["medium"], row["initial_glucose_g_L"],
                     row["initial_pH"], row["temperature_C"], row["agitation_rpm"], row["vessel_type"],
                     row["vessel_volume_mL"], row["working_volume_mL"], row["inoculum_unit"])
        study = row["study_id"]
        if study in condition_signatures and condition_signatures[study] != signature:
            raise ValueError(f"同一 study_id 的生产培养条件不一致：{culture_id}")
        condition_signatures[study] = signature
        culture = dict(row)
        culture["_start"] = start
        culture["_replicate"] = rep
        cultures[culture_id] = culture
        cultures_by_candidate.setdefault(key, []).append(culture_id)
    if set(cultures_by_candidate) != set(candidates):
        raise ValueError("每个候选菌株必须至少对应一个独立生产培养物")

    measurements = {}
    measurements_by_culture = {}
    assays = set()
    measurement_keys = set()
    for row in measurement_rows:
        measurement_id, assay_id = row["measurement_id"], row["assay_id"]
        if measurement_id in measurements or assay_id in assays:
            raise ValueError(f"measurement_id 和 assay_id 必须唯一：{measurement_id}")
        assays.add(assay_id)
        culture = cultures.get(row["culture_id"])
        if culture is None:
            raise ValueError(f"measurement 找不到 culture_id：{measurement_id}")
        expected = (culture["study_id"], culture["population_id"], culture["isolate_id"], culture["biological_replicate"], culture["blind_code"])
        observed = (row["study_id"], row["population_id"], row["isolate_id"], row["biological_replicate"], row["blind_code"])
        if observed != expected:
            raise ValueError(f"measurement 与 culture/candidate 身份字段冲突：{measurement_id}")
        value = number(row["value"], "measurement value")
        time_h = number(row["recorded_time_h"], "recorded_time_h")
        if time_h < 0:
            raise ValueError(f"recorded_time_h 不能为负：{measurement_id}")
        if row["variable"] not in {"od600", TARGET_VARIABLE}:
            raise ValueError(f"本入口只接收 od600 和 fumarate_titer：{measurement_id}")
        if row["variable"] == "od600" and (row["unit"] != "OD600" or value <= 0):
            raise ValueError(f"OD600 必须使用 OD600 单位并为正值：{measurement_id}")
        if row["variable"] == TARGET_VARIABLE and (row["unit"] != "g/L" or value < 0):
            raise ValueError(f"fumarate_titer 必须使用 g/L 单位并为非负数：{measurement_id}")
        technical = positive_int(row["technical_replicate"], "technical_replicate")
        row_key = (row["culture_id"], row["variable"], time_h, technical)
        if row_key in measurement_keys:
            raise ValueError(f"同一培养物时间点技术重复号重复：{row_key}")
        measurement_keys.add(row_key)
        sampled = parse_time(row["sampled_at_utc"], "sampled_at_utc")
        measured = parse_time(row["measured_at_utc"], "measured_at_utc")
        available = parse_time(row["result_available_at_utc"], "result_available_at_utc")
        if sampled < culture["_start"] or measured < sampled or available < measured:
            raise ValueError(f"样品、测量和结果可用时刻顺序错误：{measurement_id}")
        elapsed = (sampled - culture["_start"]).total_seconds() / 3600
        if abs(elapsed - time_h) > TIME_TOLERANCE_H:
            raise ValueError(f"记录小时数与样品采集时间不匹配：{measurement_id}")
        validate_optional_sha(row, measurement_id)
        parsed = dict(row)
        parsed.update(_value=value, _time_h=time_h, _available=available, _sampled=sampled)
        measurements[measurement_id] = parsed
        measurements_by_culture.setdefault(row["culture_id"], []).append(parsed)

    if stage == "freeze" and any(row["variable"] == TARGET_VARIABLE for row in measurement_rows):
        raise ValueError("freeze 阶段禁止导入任何富马酸终点结果")
    if stage == "evaluate" and not any(row["variable"] == TARGET_VARIABLE for row in measurement_rows):
        raise ValueError("evaluate 阶段缺少 57 h 富马酸终点结果")

    early_by_candidate = {}
    for candidate_key, culture_ids in cultures_by_candidate.items():
        required_ids = set()
        for culture_id in culture_ids:
            od_rows = [row for row in measurements_by_culture.get(culture_id, []) if row["variable"] == "od600"]
            initial = [row for row in od_rows if row["_time_h"] == 0]
            early = [row for row in od_rows if 0 < row["_time_h"] <= EARLY_CUTOFF_H]
            if not initial or not early:
                raise ValueError(f"每个培养物必须有 t=0 和 18 h 内实测 OD：{culture_id}")
            inoculum_value = number(cultures[culture_id]["inoculum_value"], "inoculum_value")
            initial_mean = math.fsum(row["_value"] for row in initial) / len(initial)
            if not math.isclose(initial_mean, inoculum_value, rel_tol=0.05, abs_tol=1e-8):
                raise ValueError(f"培养物接种 OD 与 t=0 实测 OD 相差超过 5%：{culture_id}")
            last_time = max(row["_time_h"] for row in early)
            latest = [row for row in early if row["_time_h"] == last_time]
            if math.fsum(row["_value"] for row in latest) <= 0:
                raise ValueError(f"18 h 内早期 OD 无效：{culture_id}")
            required_ids.update(row["measurement_id"] for row in initial + latest)
            if stage == "evaluate":
                endpoints = [row for row in measurements_by_culture.get(culture_id, [])
                             if row["variable"] == TARGET_VARIABLE and row["_time_h"] == ENDPOINT_H]
                if not endpoints:
                    raise ValueError(f"每个候选培养物都必须有来源定义的 57 h 富马酸滴度：{culture_id}")
                if any(row["_available"] > now for row in endpoints):
                    raise ValueError(f"终点结果尚未到达 result_available_at_utc：{culture_id}")
        early_by_candidate[candidate_key] = required_ids

    decision_ids = set()
    decision_groups = {}
    for row in decision_rows:
        decision_id = row["decision_id"]
        if decision_id in decision_ids:
            raise ValueError(f"decision_id 重复：{decision_id}")
        decision_ids.add(decision_id)
        key = (row["study_id"], row["population_id"], row["isolate_id"])
        if key not in candidates:
            raise ValueError(f"decision 找不到候选菌株：{decision_id}")
        candidate = candidates[key]
        if row["blind_code"] != candidate["blind_code"]:
            raise ValueError(f"decision 盲码与候选菌株不一致：{decision_id}")
        method = row["selection_method"]
        if method not in METHODS:
            raise ValueError(f"不支持的 selection_method：{method}")
        score = number(row["score"], "score")
        rank = positive_int(row["rank"], "rank")
        budget = positive_int(row["assay_budget"], "assay_budget")
        selected = boolean(row["selected"], "selected")
        cutoff = number(row["input_cutoff_h"], "input_cutoff_h")
        if cutoff != EARLY_CUTOFF_H:
            raise ValueError(f"input_cutoff_h 必须为冻结的 {EARLY_CUTOFF_H:g} h")
        if row["tie_break_rule"] != TIE_BREAK_RULE:
            raise ValueError(f"tie_break_rule 必须为 {TIE_BREAK_RULE}：{decision_id}")
        decision_time = parse_time(row["decision_time_utc"], "decision_time_utc")
        if decision_time > now:
            raise ValueError(f"decision_time_utc 晚于当前验证时钟：{decision_id}")
        target_used = boolean(row["endpoint_used_to_score"], "endpoint_used_to_score")
        if target_used:
            raise ValueError(f"选株分数明确使用了终点标签：{decision_id}")
        try:
            input_ids = json.loads(row["score_input_measurement_ids"])
        except json.JSONDecodeError as exc:
            raise ValueError(f"score_input_measurement_ids 必须是 JSON 字符串数组：{decision_id}") from exc
        if not isinstance(input_ids, list) or any(not isinstance(item, str) for item in input_ids) or len(input_ids) != len(set(input_ids)):
            raise ValueError(f"score_input_measurement_ids 必须是无重复 ID 数组：{decision_id}")
        if method == "random_baseline":
            if input_ids or not row["random_seed"]:
                raise ValueError(f"random_baseline 必须记录随机种子且不引用测量值：{decision_id}")
        else:
            if row["random_seed"]:
                raise ValueError(f"非随机规则不应填写 random_seed：{decision_id}")
            expected_ids = early_by_candidate[key]
            if set(input_ids) != expected_ids:
                raise ValueError(f"选株输入必须且只能引用该菌株各培养物的初始 OD 和最后一个 18 h 内实测 OD：{decision_id}")
            for measurement_id in input_ids:
                measurement = measurements.get(measurement_id)
                if measurement is None or measurement["variable"] != "od600" or measurement["_time_h"] > EARLY_CUTOFF_H:
                    raise ValueError(f"模型输入引用了不存在或超出早期窗口的记录：{decision_id}/{measurement_id}")
                if measurement["_available"] > decision_time:
                    raise ValueError(f"选株决定早于输入 OD 结果可用时刻：{decision_id}")
            if method == "fixed_ridge" and row["model_version"] != MODEL_VERSION:
                raise ValueError(f"fixed_ridge 必须使用 {MODEL_VERSION}：{decision_id}")
        if not re.fullmatch(r"[0-9a-fA-F]{64}", row["prediction_file_sha256"]):
            raise ValueError(f"prediction_file_sha256 格式错误：{decision_id}")
        group_key = (row["study_id"], row["population_id"], method)
        decision_groups.setdefault(group_key, []).append({**row, "_score": score, "_rank": rank,
            "_budget": budget, "_selected": selected, "_time": decision_time, "_input_ids": input_ids})

    for group_key, group_rows in decision_groups.items():
        study, population, method = group_key
        expected_isolates = by_population[(study, population)]
        observed_isolates = {row["isolate_id"] for row in group_rows}
        if observed_isolates != expected_isolates or len(group_rows) != len(expected_isolates):
            raise ValueError(f"每个 selection_method 必须恰好覆盖该群体全部候选：{group_key}")
        if {row["_rank"] for row in group_rows} != set(range(1, len(group_rows) + 1)):
            raise ValueError(f"rank 必须为 1 到候选数的唯一序号：{group_key}")
        budgets = {row["_budget"] for row in group_rows}
        times = {row["_time"] for row in group_rows}
        hashes = {row["prediction_file_sha256"] for row in group_rows}
        prediction_paths = {row["prediction_file_path"] for row in group_rows}
        versions = {row["model_version"] for row in group_rows}
        if len(budgets) != 1 or next(iter(budgets)) > len(group_rows):
            raise ValueError(f"assay_budget 在组内必须相同且不超过候选数：{group_key}")
        if len(times) != 1 or len(hashes) != 1 or len(prediction_paths) != 1 or len(versions) != 1:
            raise ValueError(f"同一筛选组必须共享决定时刻、预测文件、文件摘要和模型版本：{group_key}")
        selected_rows = [row for row in group_rows if row["_selected"]]
        if len(selected_rows) != next(iter(budgets)):
            raise ValueError(f"selected 数量与 assay_budget 不相等：{group_key}")
        order = sorted(group_rows, key=lambda row: (-row["_score"], row["isolate_id"]))
        if [row["isolate_id"] for row in sorted(group_rows, key=lambda row: row["_rank"])] != [row["isolate_id"] for row in order]:
            raise ValueError(f"rank 必须按 score 降序排列，同分按 isolate_id 排序：{group_key}")
        if {row["isolate_id"] for row in sorted(group_rows, key=lambda row: row["_rank"])[:next(iter(budgets))]} != {row["isolate_id"] for row in selected_rows}:
            raise ValueError(f"selected 必须对应 rank 前 assay_budget 名：{group_key}")
        seeds = {row["random_seed"] for row in group_rows}
        if method == "random_baseline" and len(seeds) != 1:
            raise ValueError(f"同一随机基线组必须使用同一个随机种子：{group_key}")
    for population in by_population:
        found = {method for study, pop, method in decision_groups if (study, pop) == population}
        if found != {"fixed_ridge", "early_od_mean", "early_od_fold", "random_baseline"}:
            raise ValueError(f"每个群体必须同时记录固定模型、两种 OD 基线和随机基线：{population}")
    if not decision_groups:
        raise ValueError("缺少选株决策记录")

    freeze_doc = None
    if stage == "freeze":
        for culture_id, records in measurements_by_culture.items():
            if any(row["variable"] != "od600" for row in records):
                raise ValueError("freeze 阶段只允许导入 18 h 截止内的 OD 记录")
            if any(row["_time_h"] > EARLY_CUTOFF_H for row in records):
                raise ValueError(f"freeze 阶段不允许导入 18 h 后的测量：{culture_id}")
        decisions_latest = max(row["_time"] for rows in decision_groups.values() for row in rows)
        if decisions_latest > now:
            raise ValueError("尚未到选株决定时间")
        receipt_inputs = {
            "candidates_csv_sha256": None,
            "cultures_csv_sha256": None,
            "early_measurements_sha256": rows_sha256([row for row in measurement_rows if row["variable"] == "od600" and number(row["recorded_time_h"], "recorded_time_h") <= EARLY_CUTOFF_H]),
            "decisions_csv_sha256": None,
            "frozen_model_sha256": model_sha,
        }
        freeze_doc = {"schema_version": 1, "record_type": "pre_endpoint_selection_freeze",
                      "frozen_at_utc": now.isoformat().replace("+00:00", "Z"),
                      "timestamp_trust": "local_system_clock; not a trusted timestamp or digital signature",
                      "study_id": candidate_rows[0]["study_id"], "candidate_count": len(candidate_rows),
                      "population_count": len(by_population), "decision_row_count": len(decision_rows),
                      "inputs": receipt_inputs}
        return {"stage": stage, "freeze_receipt": freeze_doc,
                "integrity_passed": True, "endpoint_labels_present": False,
                "independent_blinding_proven": False,
                "strict_external_validation_allowed": False,
                "strict_validation_blockers": strict_blockers(),
                "model_sha256": model_sha}

    if stage != "evaluate":
        raise ValueError(f"未知 stage：{stage}")
    if not freeze_receipt or freeze_receipt.get("record_type") != "pre_endpoint_selection_freeze":
        raise ValueError("evaluate 阶段必须提供 freeze 生成的决策冻结回执")
    receipt_time = parse_time(freeze_receipt["frozen_at_utc"], "freeze_receipt.frozen_at_utc")
    if freeze_receipt.get("inputs", {}).get("frozen_model_sha256") != model_sha:
        raise ValueError("冻结回执绑定的模型 SHA 与当前锁定模型不一致")
    decisions_latest = max(row["_time"] for rows in decision_groups.values() for row in rows)
    if decisions_latest > receipt_time:
        raise ValueError("冻结回执时间早于其中记录的选株决定")
    early_hash = rows_sha256([row for row in measurement_rows if row["variable"] == "od600" and number(row["recorded_time_h"], "recorded_time_h") <= EARLY_CUTOFF_H])
    if freeze_receipt["inputs"].get("early_measurements_sha256") != early_hash:
        raise ValueError("evaluate 阶段 18 h 内 OD 输入与冻结回执不一致")
    endpoint_rows = [row for row in measurements.values() if row["variable"] == TARGET_VARIABLE]
    for row in endpoint_rows:
        if row["_time_h"] != ENDPOINT_H:
            raise ValueError(f"只接受固定 57 h 富马酸终点：{row['measurement_id']}")
        if row["_available"] <= receipt_time:
            raise ValueError(f"终点结果在冻结回执生成前已可用，不能作为前瞻盲测：{row['measurement_id']}")
    for group_key, rows in decision_groups.items():
        population_endpoints = [row for row in endpoint_rows if (row["study_id"], row["population_id"]) == group_key[:2]]
        if not population_endpoints:
            raise ValueError(f"群体缺少富马酸终点结果：{group_key[:2]}")
        if any(max(row["_time"] for row in rows) >= measurement["_available"] for measurement in population_endpoints):
            raise ValueError(f"选株决定必须早于该群体所有终点结果可用时间：{group_key}")
    if len({row["assay_method"] for row in endpoint_rows}) != 1:
        raise ValueError("同一验证队列的富马酸终点必须使用单一、预先声明的 assay_method")
    if omics_rows:
        validate_omics(omics_rows, cultures, candidates)
    return {
        "schema_version": 1,
        "stage": stage,
        "generated_at_utc": now.isoformat().replace("+00:00", "Z"),
        "study_id": candidate_rows[0]["study_id"],
        "candidate_count": len(candidate_rows),
        "population_count": len(by_population),
        "culture_count": len(cultures),
        "measurement_count": len(measurement_rows),
        "early_od_measurement_count": sum(row["variable"] == "od600" for row in measurement_rows),
        "endpoint_measurement_count": len(endpoint_rows),
        "endpoint_time_h": ENDPOINT_H,
        "endpoint_unit": "g/L",
        "recovered_process_conditions_match": True,
        "decision_methods": sorted(METHODS),
        "freeze_receipt_sha256": hashlib.sha256(json.dumps(freeze_receipt, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest(),
        "frozen_model_sha256": model_sha,
        "data_integrity_passed": True,
        "recorded_temporal_order_passed": True,
        "independent_blinding_proven": False,
        "strict_external_validation_allowed": False,
        "strict_validation_blockers": strict_blockers(),
        "wet_lab_gain_validated": False,
        "score_recomputed_by_intake": False,
        "timestamp_trust": "local system clock and recorded LIMS timestamps; not independently attested",
        "interpretation": "输入关系、18 h 特征窗口、决策冻结与记录的 57 h 结果时序通过结构核验；回执不能证明分析人员此前未接触终点，也不计算预测效果或解除来源工艺阻塞。",
    }


def strict_blockers() -> list[str]:
    alignment = json.loads(SOURCE_ALIGNMENT_PATH.read_text(encoding="utf-8"))
    return list(alignment["unresolved_for_strict_process_reproduction"])


def freeze_inputs(candidate_path: Path, culture_path: Path, decision_path: Path) -> dict:
    return {
        "candidates_csv_sha256": sha256(candidate_path),
        "cultures_csv_sha256": sha256(culture_path),
        "decisions_csv_sha256": sha256(decision_path),
    }


def run(stage: str, candidate_path: Path, culture_path: Path, measurement_path: Path,
        decision_path: Path, output_path: Path, freeze_receipt_path: Path | None = None,
        omics_path: Path | None = None) -> dict:
    candidates = read_table(candidate_path, "candidates")
    cultures = read_table(culture_path, "cultures")
    measurements = read_table(measurement_path, "measurements")
    decisions = read_table(decision_path, "decisions")
    for row in decisions:
        prediction_path = Path(row["prediction_file_path"])
        if not prediction_path.is_absolute():
            prediction_path = decision_path.parent / prediction_path
        if not prediction_path.is_file():
            raise ValueError(f"找不到预测文件：{row['prediction_file_path']}")
        if sha256(prediction_path).lower() != row["prediction_file_sha256"].lower():
            raise ValueError(f"预测文件 SHA256 不匹配：{row['decision_id']}")
    omics = read_table(omics_path, "omics") if omics_path else None
    raw_verified, raw_missing = validate_raw_files(measurements, measurement_path.parent, "measurement")
    if omics:
        omics_verified, omics_missing = validate_raw_files(omics, omics_path.parent, "omics")
    else:
        omics_verified, omics_missing = 0, 0
    receipt = json.loads(freeze_receipt_path.read_text(encoding="utf-8")) if freeze_receipt_path else None
    result = validate_tables(candidates, cultures, measurements, decisions, stage, receipt, omics)
    if stage == "evaluate":
        result["raw_file_traceability"] = {
            "measurement_raw_files_sha_verified": raw_verified,
            "measurement_rows_without_raw_file": raw_missing,
            "omics_raw_files_sha_verified": omics_verified,
            "omics_rows_without_raw_file": omics_missing,
        }
    if output_path.exists():
        raise ValueError(f"回执输出路径必须是新文件：{output_path}")
    if stage == "freeze":
        result["freeze_receipt"]["inputs"].update(freeze_inputs(candidate_path, culture_path, decision_path))
        payload = result["freeze_receipt"]
    else:
        for key, path in (("candidates_csv_sha256", candidate_path), ("cultures_csv_sha256", culture_path),
                          ("decisions_csv_sha256", decision_path)):
            if receipt["inputs"].get(key) != sha256(path):
                raise ValueError(f"evaluate 阶段的 {key} 与冻结回执不一致")
        payload = result
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return result


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="ALE 外部队列前瞻数据入口和盲测时序核验")
    parser.add_argument("--stage", choices=("freeze", "evaluate"), required=True)
    parser.add_argument("--candidates", type=Path, default=ROOT / "research/templates/ale_prospective_candidates.csv")
    parser.add_argument("--cultures", type=Path, default=ROOT / "research/templates/ale_prospective_cultures.csv")
    parser.add_argument("--measurements", type=Path, default=ROOT / "research/templates/ale_prospective_measurements.csv")
    parser.add_argument("--decisions", type=Path, default=ROOT / "research/templates/ale_prospective_decisions.csv")
    parser.add_argument("--omics", type=Path)
    parser.add_argument("--freeze-receipt", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        if args.stage == "evaluate" and args.freeze_receipt is None:
            raise ValueError("evaluate 阶段必须通过 --freeze-receipt 指定冻结回执")
        result = run(args.stage, args.candidates, args.cultures, args.measurements,
                     args.decisions, args.output, args.freeze_receipt, args.omics)
        print(json.dumps({key: result[key] for key in result if key != "freeze_receipt"}, ensure_ascii=False, indent=2))
        print(f"回执：{args.output}")
        return 0
    except Exception as exc:
        print(f"失败：{exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
