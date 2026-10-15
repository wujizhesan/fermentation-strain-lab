from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import html
import json
from pathlib import Path
import sys
import urllib.request

import project

ROOT = Path(__file__).resolve().parent
MODULES = (*project.MODULES, "external")


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path: Path):
    return project.read_json(path)


def package_integrity(root: Path = ROOT) -> dict:
    path = root / "PACKAGE_MANIFEST.json"
    if not path.exists():
        return {"available": False, "passed": None, "files_checked": 0}
    manifest = read_json(path)
    for rel, record in manifest["files"].items():
        source = project.relative_file(root, rel)
        if source.stat().st_size != record["bytes"] or digest(source) != record["sha256"]:
            raise ValueError("发行文件发生变化：" + rel)
    return {"available": True, "passed": True, "files_checked": len(manifest["files"]), "version": manifest["version"]}


def reference_status(root: Path = ROOT) -> dict:
    package = package_integrity(root)
    modules = {}
    for name in ("screening", "softsensor"):
        folder = root / project.ARCHIVES[name]
        report, proof = project.checked_result(folder, root / project.PROTOCOLS[name])
        summary = project.screening_summary(report) if name == "screening" else project.softsensor_summary(report)
        modules[name] = {"title": project.TITLES[name], "summary": summary, "evidence": [proof], "snapshot_checked": True, "numeric_recomputed_now": False}
    from build_marge_viewer import load_viewer_data
    viewer = load_viewer_data(root / "research/results")
    modules["metabolism"] = {"title": project.TITLES["metabolism"], "snapshot_checked": True, "numeric_recomputed_now": False,
        "summary": {"conclusion": "比较上游原样与保留蛋白组平衡版本；未验证实测通量或菌株增产。", "endpoint": "模型可行通量区间", "cases": viewer["summary"]["cases"], "strict_paper_reproduction": False}}
    folder = root / "research/results/ale_external_screen_v1"
    report, proof = project.checked_result(folder, root / "research/ale_external_screen_protocol_20261015.json")
    if any(report["scope"][key] is not False for key in ("new_blind_validation", "original_18h_fumarate_ridge_external_validation", "publication_superiority_claim")):
        raise ValueError("外部案例证据范围发生变化")
    modules["external"] = {"title": "外部公开案例回放", "snapshot_checked": True, "numeric_recomputed_now": False,
        "summary": {"endpoint": "各面板内部归一化终点", "metrics": report["metrics"], "conclusion": "两项24h筛选规则的回顾性回放；原18h富马酸模型未获外部验证。"}, "evidence": [proof]}
    input_status = {"ready": None, "missing": [], "changed": []}
    source_path = root / "RELEASE_SOURCES.json"
    if source_path.exists():
        for row in read_json(source_path)["files"]:
            path = project.relative_file(root, row["path"])
            if not path.is_file():
                input_status["missing"].append(row["path"])
            elif row["validation"] == "archive_members":
                from marge_sources import archive_manifest, validate_archive_members
                validate_archive_members(path.read_bytes(), archive_manifest(root))
            elif digest(path) != row["sha256"]:
                input_status["changed"].append(row["path"])
        input_status["ready"] = not input_status["missing"] and not input_status["changed"]
    return {"project": "Fermentation Strain Lab", "schema_version": 1, "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "mode": "reference_snapshot_integrity", "package": package, "raw_inputs": input_status, "modules": modules,
        "same_strain_end_to_end_validated": False, "strain_optimization_gain_validated": False, "original_model_external_validation": False}


def external_summary(folder: Path) -> dict:
    from verify_external_screen_validation import verify
    checks = verify(folder / "report.json")
    if not all(checks.values()):
        raise ValueError("外部案例复算不通过：" + ",".join(k for k, v in checks.items() if not v))
    report, proof = project.checked_result(folder, ROOT / "research/ale_external_screen_protocol_20261015.json")
    return {"title": "外部公开案例回放", "numeric_recomputed_now": True, "summary": {
        "endpoint": "各面板内部归一化终点", "metrics": report["metrics"], "conclusion": "回顾性24h规则回放，不验证原18h富马酸模型。"}, "evidence": [proof]}


def render(data: dict) -> str:
    cards = []
    for name, module in data["modules"].items():
        summary = module["summary"]
        table = ""
        if "methods" in summary:
            rows = []
            for row in summary["methods"]:
                value = row.get("expected_selected_titer_g_L", row.get("batch_macro_RMSE_g_L"))
                rows.append(f"<tr><td>{html.escape(row['label'])}</td><td>{value:.6f}</td></tr>")
            table = "<table><thead><tr><th>方法</th><th>g/L</th></tr></thead><tbody>" + "".join(rows) + "</tbody></table>"
        elif "metrics" in summary:
            rows = [f"<tr><td>{html.escape(row['panel_id'])}</td><td>{html.escape(row['method'])}</td><td>{100*row['fraction_of_panel_best']:.2f}%</td></tr>" for row in summary["metrics"]]
            table = "<div class='scroll'><table><thead><tr><th>面板</th><th>规则</th><th>选中终点／面板最优</th></tr></thead><tbody>" + "".join(rows) + "</tbody></table></div>"
        else:
            table = f"<p>情景数：{summary.get('cases', 24)}；严格论文复现未完成。</p>"
        cards.append(f"<section><h2>{html.escape(module['title'])}</h2><p>{html.escape(summary['conclusion'])}</p>{table}</section>")
    mode = "本次计算并独立核验" if data["mode"] == "fresh_replay" else "参考快照摘要检查；本次未重新计算数值"
    return "<!doctype html><html lang='zh-CN'><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>Fermentation Strain Lab</title><style>body{font:16px/1.6 system-ui,sans-serif;margin:0;background:#f3f6f7;color:#17323b}main{max-width:980px;margin:auto;padding:28px 18px}h1{font-size:30px}section{background:white;padding:22px;margin:20px 0;border-radius:12px}table{border-collapse:collapse;width:100%;font-size:14px}th,td{text-align:left;padding:9px;border-bottom:1px solid #dae2e5}.scroll{overflow:auto}.note{padding:16px;background:#e7eef2;border-radius:10px}code{overflow-wrap:anywhere}</style><main><h1>菌株筛选与发酵建模研究原型</h1><p>主线：早期观测 → 菌株排序 → 终点评价。</p><p class='note'>" + mode + "。三个研究模块使用不同数据；尚无同株实验闭环、原模型外部验证或菌株增产证据。</p>" + "".join(cards) + "<p>完整数值与来源记录见同目录 summary.json。</p></main></html>"


def save_overview(data: dict, output: Path) -> None:
    project.write_json(output / "summary.json", data)
    (output / "index.html").write_text(render(data), encoding="utf-8")


def prepare(download: bool = False) -> dict:
    package_integrity()
    source_file = ROOT / "RELEASE_SOURCES.json"
    if not source_file.exists():
        raise ValueError("此命令用于精简发行目录；原工作目录直接使用既有来源入口")
    sources = read_json(source_file)
    missing = []
    for row in sources["files"]:
        path = project.relative_file(ROOT, row["path"])
        if path.is_file():
            if row["validation"] == "sha256" and digest(path) != row["sha256"]:
                raise ValueError("冻结输入发生变化，拒绝覆盖：" + row["path"])
            continue
        if row["validation"] == "archive_member" or row["validation"] == "archive_members":
            continue
        if not download:
            missing.append(row["path"])
            continue
        if not row.get("url"):
            raise ValueError("必须由离线输入包恢复冻结快照：" + row["path"])
        print("获取公开输入：" + row["path"], flush=True)
        chunks = []
        with urllib.request.urlopen(row["url"], timeout=30) as response:
            while chunk := response.read(65536):
                chunks.append(chunk)
        body = b"".join(chunks)
        if hashlib.sha256(body).hexdigest() != row["sha256"]:
            raise ValueError("公开服务内容与冻结SHA不符；未写入：" + row["path"])
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("xb") as stream:
            stream.write(body)
    from prepare_marge_demo import prepare as prepare_marge
    if not missing:
        prepare_marge(check_only=not download)
    else:
        raise FileNotFoundError("缺少冻结输入；使用 prepare --download 或离线输入包：" + ";".join(missing))
    status = reference_status()["raw_inputs"]
    if status["ready"] is not True:
        raise ValueError("输入未完整准备：" + json.dumps(status, ensure_ascii=False))
    return status


def collect_run(folder: Path) -> dict:
    receipt = read_json(folder / "lab_receipt.json")
    if receipt["status"] != "passed":
        raise ValueError("运行未完成或未通过")
    for name in ("summary.json", "index.html"):
        if digest(folder / name) != receipt["artifacts"][name]:
            raise ValueError("运行总览发生变化：" + name)
    selected = receipt["modules_requested"]
    if not selected or len(set(selected)) != len(selected) or any(name not in MODULES for name in selected):
        raise ValueError("运行模块集合不合法")
    modules = {}
    if receipt.get("core_run_receipt_sha256"):
        if digest(folder / "core/run_receipt.json") != receipt["core_run_receipt_sha256"]:
            raise ValueError("核心运行回执变化")
        core = project.collect_run(folder / "core")
        if core["overall_integrity_passed"] is not True:
            raise ValueError("核心模块证据检查失败")
        modules.update(core["modules"])
    for step in receipt["external_steps"]:
        log = project.relative_file(folder, step["log"])
        if step["exit_code"] != 0 or digest(log) != step["log_sha256"]:
            raise ValueError("外部案例执行日志变化")
    if "external" in selected:
        if [step["step"] for step in receipt["external_steps"]] != ["external_calculate", "external_verify"]:
            raise ValueError("外部案例计算／核验步骤缺失")
        modules["external"] = external_summary(folder / "external")
    elif receipt["external_steps"]:
        raise ValueError("未申请外部案例却包含外部执行步骤")
    if set(modules) != set(selected):
        raise ValueError("实际运行模块与申请模块不一致")
    data = read_json(folder / "summary.json")
    if data["mode"] != "fresh_replay" or any(data.get(key) is not False for key in ("same_strain_end_to_end_validated", "strain_optimization_gain_validated", "original_model_external_validation")):
        raise ValueError("运行报告证据范围变化")
    if set(data["modules"]) != set(modules):
        raise ValueError("总览模块与实际结果不一致")
    for name, module in modules.items():
        saved = data["modules"][name]
        for field in ("title", "summary", "evidence"):
            if saved.get(field) != module[field]:
                raise ValueError("总览与模块结果不一致：" + name + "/" + field)
        if name != "external":
            if saved.get("integrity_passed") is not True or saved.get("attachments") != module["attachments"]:
                raise ValueError("核心模块展示与核验状态不一致：" + name)
        elif saved.get("numeric_recomputed_now") is not True:
            raise ValueError("外部案例运行状态不一致")
    if (folder / "index.html").read_text(encoding="utf-8") != render(data):
        raise ValueError("HTML总览与核验汇总不一致")
    return data


def run(output: Path, selected: tuple, interpreter: str) -> dict:
    package_integrity()
    output = project.check_output(output)
    core_selected = tuple(name for name in selected if name in project.MODULES)
    project.preflight(interpreter, tuple(set(core_selected) | ({"softsensor"} if "external" in selected else set())))
    output.mkdir(parents=True)
    receipt = {"status": "running", "modules_requested": list(selected), "external_steps": [], "core_run_receipt_sha256": None}
    data = {"project": "Fermentation Strain Lab", "schema_version": 1, "mode": "fresh_replay", "modules": {},
        "same_strain_end_to_end_validated": False, "strain_optimization_gain_validated": False, "original_model_external_validation": False}
    try:
        if core_selected:
            core = project.run_project(output / "core", core_selected, interpreter)
            data["modules"].update(core["modules"])
            receipt["core_run_receipt_sha256"] = digest(output / "core/run_receipt.json")
        if "external" in selected:
            steps = []
            try:
                for label, arguments in (("external_calculate", ["external_screen_validation.py", "--output", str(output / "external")]), ("external_verify", ["verify_external_screen_validation.py", "--result", str(output / "external/report.json")])):
                    project.execute_step(interpreter, arguments, output, label, steps)
            finally:
                receipt["external_steps"] = [{**step, "log": Path(step["log"]).relative_to(output).as_posix()} for step in steps]
            data["modules"]["external"] = external_summary(output / "external")
        save_overview(data, output)
        receipt.update(status="passed", artifacts={name: digest(output / name) for name in ("summary.json", "index.html")})
    except Exception as exc:
        receipt.update(status="failed", error=str(exc))
        project.write_json(output / "lab_receipt.json", receipt)
        raise
    project.write_json(output / "lab_receipt.json", receipt)
    collect_run(output)
    return data


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description="公开数据驱动的菌株筛选与发酵建模研究原型")
    commands = parser.add_subparsers(dest="command", required=True)
    status = commands.add_parser("status")
    status.add_argument("--json", action="store_true")
    status.add_argument("--run-dir", type=Path)
    preview = commands.add_parser("report")
    preview.add_argument("--output", type=Path, required=True)
    preview.add_argument("--run-dir", type=Path)
    preparation = commands.add_parser("prepare")
    preparation.add_argument("--download", action="store_true")
    execution = commands.add_parser("run")
    execution.add_argument("--module", choices=("all", *MODULES), default="screening")
    execution.add_argument("--output", type=Path, required=True)
    execution.add_argument("--python")
    verification = commands.add_parser("verify")
    verification.add_argument("--run-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == "prepare":
            print(json.dumps(prepare(args.download), ensure_ascii=False, indent=2))
            return 0
        if args.command in ("status", "report"):
            data = collect_run(args.run_dir) if args.run_dir else reference_status()
            if args.command == "report":
                output = project.check_output(args.output)
                output.mkdir(parents=True)
                save_overview(data, output)
                print("总览：" + str(output / "index.html"))
                return 0
        elif args.command == "verify":
            collect_run(args.run_dir)
            print("运行回执、日志、数值结果与证据范围检查通过")
            return 0
        else:
            selected = MODULES if args.module == "all" else (args.module,)
            data = run(args.output, selected, project.resolve_interpreter(args.python))
        if getattr(args, "json", False):
            print(json.dumps(data, ensure_ascii=False, indent=2))
        else:
            print("科研原型；" + ("参考快照检查，本次未重算" if data["mode"] == "reference_snapshot_integrity" else "本次计算并独立核验"))
            for module in data["modules"].values():
                print(module["title"] + "：" + module["summary"]["conclusion"])
            if data.get("raw_inputs", {}).get("ready") is False:
                print("重算所需原始输入尚未完整取得；运行 prepare --download 或恢复离线输入包")
        return 0
    except Exception as exc:
        print("失败：" + str(exc), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
