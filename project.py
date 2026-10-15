import argparse
from datetime import datetime, timezone
import hashlib
import html
import json
import math
import os
from pathlib import Path, PureWindowsPath
import shutil
import subprocess
import sys
from urllib.parse import quote


ROOT = Path(__file__).resolve().parent
MODULES = ('screening', 'softsensor', 'metabolism')
TITLES = {'screening': '早期数据选株', 'softsensor': '发酵软测量', 'metabolism': '组学代谢建模'}
PROTOCOLS = {'screening': 'research/ale_lineage_screen_protocol.json',
             'softsensor': 'research/ambr_dynamic_protocol.json',
             'metabolism': 'research/ale_marge_protocol.json'}
ARCHIVES = {'screening': 'research/results/ale_lineage_screen_v1',
            'softsensor': 'research/results/ambr_dynamic', 'metabolism': 'research/results'}
METHOD_NAMES = {'uniform_random_expectation': '随机选择（期望）', 'early_density': '早期OD',
                'early_density_fold': '早期OD倍数', 'fixed_ridge': '固定岭回归',
                'training_mean': '训练标签均值', 'time_recipe_ridge': '时间＋计划配方',
                'instantaneous_ridge': '当前在线信号', 'history_ridge': '当前信号＋历史窗口'}


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def reject_constant(value):
    raise ValueError('JSON包含非有限数值：' + value)


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'), parse_constant=reject_constant)


def relative_file(base, name):
    if Path(name).is_absolute() or PureWindowsPath(name).drive or '..' in Path(name).parts or '..' in PureWindowsPath(name).parts:
        raise ValueError('文件清单必须使用目录内的相对路径：' + name)
    path = (Path(base) / name).resolve()
    if not path.is_relative_to(Path(base).resolve()):
        raise ValueError('结果文件超出所属目录：' + name)
    return path


def finite(value):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError('指标必须是有限数值')
    return value


def relative_change(candidate, reference):
    candidate, reference = finite(candidate), finite(reference)
    return 100 * (candidate / reference - 1) if reference != 0 else None


def checked_result(folder, protocol_path):
    folder = Path(folder)
    report = read_json(folder / 'report.json')
    receipt = read_json(folder / 'verification.json')
    if receipt.get('passed', receipt.get('verified')) is not True:
        raise ValueError('独立核验回执未通过：' + str(folder))
    if receipt.get('report_sha256') != sha(folder / 'report.json'):
        raise ValueError('核验回执与报告不匹配：' + str(folder))
    if report.get('protocol_sha256') != sha(protocol_path):
        raise ValueError('评价方案与报告不匹配：' + str(folder))
    if 'protocol_sha256' in receipt and receipt['protocol_sha256'] != report['protocol_sha256']:
        raise ValueError('核验回执与评价方案不匹配')
    artifacts = report.get('artifact_sha256')
    if not isinstance(artifacts, dict) or not artifacts:
        raise ValueError('报告缺少结果文件清单')
    for name, expected in artifacts.items():
        if sha(relative_file(folder, name)) != expected:
            raise ValueError('结果文件摘要不匹配：' + name)
    return report, {'report_path': str((folder / 'report.json').resolve()),
                    'report_sha256': sha(folder / 'report.json'),
                    'verification_sha256': sha(folder / 'verification.json'),
                    'protocol_sha256': sha(protocol_path), 'artifacts_checked': len(artifacts)}


def check_inputs(root, protocol):
    for name, expected in {**protocol.get('input_sha256', {}),
                           **protocol.get('implementation_sha256', {})}.items():
        if sha(Path(root) / name) != expected:
            raise ValueError('冻结输入或计算程序已变化：' + name)


def screening_summary(report):
    methods = report['summary']
    random = finite(methods['uniform_random_expectation']['expected_selected_titer_g_L'])
    rows = []
    for name in METHOD_NAMES:
        if name not in methods:
            continue
        row = methods[name]
        rows.append({'method': name, 'label': METHOD_NAMES[name],
                     'expected_selected_titer_g_L': finite(row['expected_selected_titer_g_L']),
                     'fraction_of_heldout_best': finite(row['fraction_of_heldout_best']),
                     'top1_hit_probability': finite(row['top1_hit_probability']),
                     'relative_to_random_percent': relative_change(row['expected_selected_titer_g_L'], random)})
    ridge = methods['fixed_ridge']['expected_selected_titer_g_L']
    conclusion = '当前岭回归的选株滴度未高于随机选择期望。' if finite(ridge) <= random else '回放观察到岭回归滴度高于随机期望，独立收益尚未验证。'
    return {'dataset': '公开酿酒酵母富马酸ALE', 'endpoint': '终点滴度，g/L',
            'isolates': report['isolates'], 'independent_population_groups': report['population_groups'],
            'comparison': '每折从留出谱系3株中选1株，三个谱系等权；随机值为期望。',
            'methods': rows, 'conclusion': conclusion, 'prospective_selection_validated': False,
            'strain_optimization_gain_validated': False,
            'limitations': ['数据此前已查看，只有三个群体，属于探索性回放。',
                            '化验采样与返回时钟未确认，尚未验证节省化验或生产率。',
                            '没有与这些分离株可靠对应的组学数据。']}


def softsensor_summary(report):
    rows = [row for row in report['summaries'] if row['direction'] == '1_to_2' and row['target'] == 'CDW_g_L']
    methods = [{'method': row['strategy'], 'label': METHOD_NAMES[row['strategy']],
                'batch_macro_RMSE_g_L': finite(row['batch_macro_RMSE'])} for row in rows]
    by_name = {row['method']: row['batch_macro_RMSE_g_L'] for row in methods}
    delta = relative_change(by_name['history_ridge'], by_name['instantaneous_ridge'])
    return {'dataset': 'Tamminen 2022公开真菌发酵', 'endpoint': '菌体干重，g/L',
            'cultures': report['cultures'], 'independent_rounds': report['independent_rounds'],
            'comparison': '第一轮训练、第二轮测试；十个测试培养的RMSE等权平均，越低越好。',
            'methods': methods, 'history_RMSE_relative_change_percent': delta,
            'conclusion': '历史窗口主目标误差未改善。' if delta >= 0 else '主方向观察到历史窗口误差降低，部署收益尚未验证。',
            'strain_optimization_gain_validated': False,
            'limitations': ['两轮实验共享条件，培养罐数不等于独立实验轮数。',
                            '实际在线信号到达时间未确认，属于条件性回放。',
                            '黏液酸标签经过终点体积校正，不能视为即时罐内浓度。']}


def load_module(name, root=ROOT, folder=None, fresh=False):
    root = Path(root)
    folder = Path(folder) if folder is not None else root / ARCHIVES[name]
    protocol_path = root / PROTOCOLS[name]
    protocol = read_json(protocol_path)
    if name == 'metabolism':
        from marge_sources import validate_marge_inputs
        validate_marge_inputs(root, protocol)
        variants, evidence = [], []
        for directory, kept in [('ale_marge_demo_v1', False), ('ale_marge_balancefix_v1', True)]:
            report, proof = checked_result(folder / directory, protocol_path)
            if report['input_sha256'] != protocol['input_sha256'] or report['proteome_balance_kept'] is not kept:
                raise ValueError('MARGE输入或版本身份不匹配')
            values = [finite(row['stage2_proteome_optimum']) for row in report['cases']]
            variants.append({'label': '保留平衡对照' if kept else '上游原样',
                             'cases': len(values), 'stage2_range': [min(values), max(values)],
                             'proteome_balance_kept': kept})
            evidence.append(proof)
        summary = {'dataset': '公开酵母亲本与进化株RNA／蛋白组', 'endpoint': '模型允许的通量范围',
                   'variants': variants, 'conclusion': '两套计算已独立复核，通量与菌株增产尚未实测验证。',
                   'strain_optimization_gain_validated': False, 'paper_exact_reproduction_completed': False,
                   'limitations': ['上游原样在当前OptLang／GLPK路径存在蛋白组平衡行覆盖。',
                                   '保留平衡版本是本地诊断对照，作者CPLEX路径未验证。',
                                   '严格论文复现缺少五类输入；未与选株分离株逐株对应。']}
    else:
        check_inputs(root, protocol)
        report, proof = checked_result(folder, protocol_path)
        evidence = [proof]
        if name == 'screening':
            if report['input_sha256'] != protocol['input_sha256']:
                raise ValueError('选株报告输入不匹配')
            summary = screening_summary(report)
        else:
            spec = protocol['source_file']
            if sha(root / spec['path']) != spec['sha256'] or report['source_sha256'] != spec['sha256']:
                raise ValueError('软测量公开工作簿摘要不匹配')
            if sha(root / 'ambr_dynamic.py') != report['implementation_sha256']:
                raise ValueError('软测量计算程序与报告不匹配')
            summary = softsensor_summary(report)
    attachments = []
    candidates = {'screening': [('选株对照图', folder / 'selection_comparison.png')],
                  'softsensor': [('逐罐误差对照图', folder / 'biomass_round_holdout.png')],
                  'metabolism': [('交互式通量区间报告', folder / 'marge-demo.html'),
                                 ('交互式通量区间报告', folder / 'ale_marge_viewer_v1/marge-demo.html')]}
    for label, path in candidates[name]:
        if path.is_file():
            attachments.append({'label': label, 'path': str(path.resolve()), 'sha256': sha(path)})
            if name == 'metabolism':
                break
    return {'title': TITLES[name], 'integrity_passed': True,
            'verification_mode': '本次运行并独立复算' if fresh else '检查既有核验回执与文件摘要',
            'summary': summary, 'evidence': evidence, 'attachments': attachments}


def collect(root=ROOT, folders=None, selected=MODULES, fresh=False):
    modules = {}
    for name in selected:
        try:
            modules[name] = load_module(name, root, (folders or {}).get(name), fresh)
        except Exception as exc:
            modules[name] = {'title': TITLES[name], 'integrity_passed': False, 'error': str(exc)}
    return {'project': '发酵分析与适应性进化菌株评价', 'schema_version': 1,
            'created_at_utc': datetime.now(timezone.utc).isoformat(),
            'overall_integrity_passed': all(module['integrity_passed'] for module in modules.values()),
            'same_strain_end_to_end_validated': False, 'new_blind_validation': False,
            'publication_superiority_claim': False, 'modules': modules}


def check_output(output, root=ROOT):
    output, root = Path(output).resolve(), Path(root).resolve()
    if output.exists():
        raise ValueError('输出目录必须是新目录：' + str(output))
    protected = [root / 'data', root / 'demo', root / 'dist']
    result_root = root / 'research/results'
    if result_root.exists():
        protected.extend(path for path in result_root.iterdir() if path.is_dir())
    for path in protected:
        if output.is_relative_to(path.resolve()):
            raise ValueError('输出目录位于数据、演示、发行包或历史结果内：' + str(output))
    return output


def write_json(path, data):
    with Path(path).open('x', encoding='utf-8') as stream:
        json.dump(data, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write('\n')


def report_html(data, output):
    cards, sections = [], []
    for name, module in data['modules'].items():
        summary = module['summary']
        cards.append('<a class="card" href="#' + name + '"><span>' + html.escape(module['title']) +
                     '</span><strong>' + html.escape(summary['conclusion']) + '</strong><small>' +
                     html.escape(module['verification_mode']) + '</small></a>')
        if name == 'screening':
            header = '<tr><th>方法</th><th>平均选中滴度 g/L</th><th>留出组最佳值比例</th><th>相对随机期望</th></tr>'
            body = ''.join('<tr><td>{}</td><td>{:.6f}</td><td>{:.2%}</td><td>{}</td></tr>'.format(
                html.escape(row['label']), row['expected_selected_titer_g_L'], row['fraction_of_heldout_best'],
                '不可计算' if row['relative_to_random_percent'] is None else '{:+.2f}%'.format(row['relative_to_random_percent']))
                for row in summary['methods'])
        elif name == 'softsensor':
            header = '<tr><th>方法</th><th>平均逐罐RMSE g/L</th></tr>'
            body = ''.join('<tr><td>{}</td><td>{:.6f}</td></tr>'.format(html.escape(row['label']), row['batch_macro_RMSE_g_L'])
                           for row in summary['methods'])
        else:
            header = '<tr><th>版本</th><th>情景数</th><th>第二阶段目标值范围</th></tr>'
            body = ''.join('<tr><td>{}</td><td>{}</td><td>{:.3f} — {:.3f}</td></tr>'.format(
                html.escape(row['label']), row['cases'], *row['stage2_range']) for row in summary['variants'])
        links = []
        for proof in module['evidence']:
            path = proof['report_path']
            url = quote(os.path.relpath(path, output).replace(os.sep, '/'), safe='/')
            links.append('<li><a href="' + html.escape(url, quote=True) + '">来源报告</a> · SHA256 <code>' +
                         html.escape(proof['report_sha256']) + '</code></li>')
        limits = ''.join('<li>' + html.escape(value) + '</li>' for value in summary['limitations'])
        attachments = []
        for artifact in module.get('attachments', []):
            url = quote(os.path.relpath(artifact['path'], output).replace(os.sep, '/'), safe='/')
            attachments.append('<a href="' + html.escape(url, quote=True) + '">' + html.escape(artifact['label']) + '</a>')
        sections.append('<section id="' + name + '"><h2>' + html.escape(module['title']) + '</h2><p>' +
                        html.escape(summary['dataset']) + ' · ' + html.escape(summary['endpoint']) + '</p><p>' +
                        html.escape(summary.get('comparison', '原样与保留平衡对照；第二阶段目标值为模型计算结果。')) +
                        '</p><div class="table"><table><thead>' + header + '</thead><tbody>' + body +
                        '</tbody></table></div><h3>结论与适用范围</h3><p>' + html.escape(summary['conclusion']) +
                        '</p><ul>' + limits + '</ul><p>' + ' · '.join(attachments) + '</p><details><summary>查看来源与计算核验记录</summary><ul>' +
                        ''.join(links) + '</ul></details></section>')
    return '''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>发酵与菌株评价 · 项目总览</title>
<style>:root{color-scheme:light}*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;background:#f3f6fa;color:#20334b;font:16px/1.7 system-ui,"Microsoft YaHei",sans-serif}main{max-width:1120px;margin:auto;padding:32px 24px}h1{font-size:clamp(26px,5vw,40px);line-height:1.25;margin:12px 0}h2{font-size:24px}h3{font-size:17px}.eyebrow{color:#267c72;font-weight:700}.intro{max-width:800px}.notice{padding:16px 20px;border-left:4px solid #ac720b;background:#fff7e6;border-radius:8px}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(240px,1fr));gap:16px;margin:24px 0}.card,section{background:white;border:1px solid #dce5ed;border-radius:14px;padding:24px}.card{text-decoration:none;color:inherit;display:flex;flex-direction:column;gap:12px}.card:hover{border-color:#267c72}.card span{color:#267c72}.card strong{font-size:18px;line-height:1.5}.card small{color:#67788d}section{margin:20px 0;scroll-margin-top:20px}section h2{margin-top:0}.table{overflow:auto}table{width:100%;border-collapse:collapse;text-align:left;white-space:nowrap;font-variant-numeric:tabular-nums}th,td{padding:12px;border-bottom:1px solid #e7edf3}th{background:#f3f7fa;font-size:14px}a{color:#176e65}a:focus-visible,summary:focus-visible{outline:3px solid #267c72;outline-offset:4px}summary{cursor:pointer}code{font-size:12px;overflow-wrap:anywhere}footer{color:#67788d;font-size:14px;padding:16px 0}@media(max-width:600px){main{padding:20px 12px}.card,section{padding:18px}.cards{grid-template-columns:1fr}th,td{padding:10px}}</style></head><body><main>
<div class="eyebrow">公开数据 · 三条研究线 · 一个入口</div><h1>发酵分析与菌株评价</h1>
<p class="intro">主问题：能否利用菌株的早期数据，选出最终产物表现更好的菌株？下方分别展示选株、软测量与代谢建模的计算证据。</p>
<div class="notice"><strong>当前证据范围：</strong>三个模块使用不同公开数据，尚未在同一套菌株实验中贯通。独立计算核验通过，不代表已验证选株收益、菌株增产或论文优势。</div>
<nav class="cards" aria-label="研究模块">''' + ''.join(cards) + '</nav>' + ''.join(sections) + '''
<footer>本页由统一入口生成，查看本页没有重新训练模型。<a href="project_summary.json">下载汇总JSON</a></footer>
</main></body></html>'''


def save_report(data, output):
    if not data['overall_integrity_passed']:
        raise ValueError('至少一个模块的证据检查未通过，不能生成通过状态的总览')
    write_json(output / 'project_summary.json', data)
    with (output / 'index.html').open('x', encoding='utf-8') as stream:
        stream.write(report_html(data, output))


def resolve_interpreter(value):
    value = value or sys.executable
    path = Path(value)
    resolved = str(path.resolve()) if path.is_file() else shutil.which(value)
    if not resolved:
        raise ValueError('找不到Python解释器；请用--python指定实际可执行文件路径')
    return resolved


def preflight(interpreter, selected):
    packages = {'screening': ['numpy', 'matplotlib'], 'softsensor': ['numpy', 'openpyxl', 'matplotlib'],
                'metabolism': ['numpy', 'scipy', 'cobra', 'optlang']}
    required = sorted({package for name in selected for package in packages[name]})
    script = 'import importlib,json,sys; names=' + repr(required) + '; [importlib.import_module(n) for n in names]; '
    if 'metabolism' in selected:
        script += "import importlib.util; assert importlib.util.find_spec('reframed'), 'Missing reframed'; "
    script += "import importlib.metadata as m; print(json.dumps({'python':sys.version,'executable':sys.executable,'package_versions':{n:m.version(n) for n in names}}))"
    completed = subprocess.run([interpreter, '-X', 'utf8', '-B', '-c', script], cwd=ROOT,
                               capture_output=True, text=True, encoding='utf-8', timeout=60)
    if completed.returncode:
        raise ValueError('计算环境检查失败，请指定对应环境的--python：\n' + completed.stderr[-1600:])
    return json.loads(completed.stdout)


def execute_step(interpreter, arguments, output, label, records):
    log = output / 'logs' / (label + '.log')
    log.parent.mkdir(exist_ok=True)
    command = [interpreter, '-X', 'utf8', '-B', *arguments]
    print('执行：' + label, file=sys.stderr, flush=True)
    environment = dict(os.environ)
    environment.update({'PYTHONUTF8': '1', 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1'})
    with log.open('x', encoding='utf-8') as stream:
        completed = subprocess.run(command, cwd=ROOT, env=environment, stdout=stream,
                                   stderr=subprocess.STDOUT, timeout=900)
    records.append({'step': label, 'command': command, 'exit_code': completed.returncode,
                    'log': str(log.resolve()), 'log_sha256': sha(log)})
    if completed.returncode:
        raise RuntimeError('步骤失败：' + label + '；请查看 ' + str(log))


def run_project(output, selected, interpreter, download=False):
    output = check_output(output)
    runtime = preflight(interpreter, selected)
    output.mkdir(parents=True, exist_ok=False)
    receipt = {'status': 'running', 'runtime': runtime, 'modules_requested': list(selected), 'steps': [],
               'same_strain_end_to_end_validated': False, 'strain_optimization_gain_validated': False}
    folders = {}
    try:
        for name in selected:
            folder = output / name
            folders[name] = folder
            if name == 'metabolism':
                execute_step(interpreter, ['run_marge_demo.py', '--output', str(folder)], output, 'metabolism_pipeline', receipt['steps'])
            else:
                producer, verifier, plotter = ('ale_lineage_screen.py', 'verify_ale_lineage_screen.py', 'plot_ale_lineage_screen.py') if name == 'screening' else ('ambr_dynamic.py', 'verify_ambr_dynamic.py', 'plot_ambr_dynamic.py')
                arguments = [producer, '--output', str(folder)]
                if name == 'softsensor' and download:
                    arguments.append('--download')
                execute_step(interpreter, arguments, output, name + '_calculate', receipt['steps'])
                execute_step(interpreter, [verifier, '--output', str(folder)], output, name + '_verify', receipt['steps'])
                execute_step(interpreter, [plotter, '--output', str(folder)], output, name + '_plot', receipt['steps'])
        data = collect(folders=folders, selected=selected, fresh=True)
        save_report(data, output)
        receipt.update({'status': 'passed', 'project_summary_sha256': sha(output / 'project_summary.json'),
                        'overview_sha256': sha(output / 'index.html')})
    except Exception as exc:
        receipt.update({'status': 'failed', 'error': str(exc)})
        write_json(output / 'run_receipt.json', receipt)
        raise
    write_json(output / 'run_receipt.json', receipt)
    return data


def collect_run(run_dir):
    run_dir = Path(run_dir).resolve()
    receipt = read_json(run_dir / 'run_receipt.json')
    if receipt.get('status') != 'passed':
        raise ValueError('只能汇总已完成并通过核验的运行目录')
    selected = receipt.get('modules_requested', [])
    if not selected or len(set(selected)) != len(selected) or any(name not in MODULES for name in selected):
        raise ValueError('运行回执中的模块列表无效')
    for field, filename in [('project_summary_sha256', 'project_summary.json'), ('overview_sha256', 'index.html')]:
        if receipt.get(field) != sha(run_dir / filename):
            raise ValueError('运行回执与汇总产物不匹配：' + filename)
    expected_steps = [step for name in selected for step in
                      ([name + '_pipeline'] if name == 'metabolism' else
                       [name + '_calculate', name + '_verify', name + '_plot'])]
    if [row['step'] for row in receipt.get('steps', [])] != expected_steps:
        raise ValueError('运行回执缺少计算、核验或展示步骤')
    for row in receipt['steps']:
        path = Path(row['log']).resolve()
        if row['exit_code'] != 0 or not path.is_relative_to(run_dir / 'logs') or sha(path) != row['log_sha256']:
            raise ValueError('运行步骤或执行日志检查失败：' + row['step'])
    data = collect(folders={name: run_dir / name for name in selected}, selected=selected)
    data['source_run_receipt_sha256'] = sha(run_dir / 'run_receipt.json')
    return data


def print_status(data):
    print(data['project'])
    for module in data['modules'].values():
        print('\n' + module['title'] + '：' + ('证据检查通过' if module['integrity_passed'] else '检查失败'))
        if not module['integrity_passed']:
            print(module['error'])
            continue
        print(module['verification_mode'])
        print(module['summary']['endpoint'])
        if 'comparison' in module['summary']:
            print(module['summary']['comparison'])
        print(module['summary']['conclusion'])
        for row in module['summary'].get('methods', []):
            value = row.get('expected_selected_titer_g_L', row.get('batch_macro_RMSE_g_L'))
            print('  {}：{:.6f} g/L'.format(row['label'], value))
    print('\n三个模块尚未在同一套菌株实验中贯通，真实菌株优化收益未验证。')


def main(argv=None):
    parser = argparse.ArgumentParser(description='统一查看、运行和汇总选株、软测量与MARGE；模型成绩与运行状态分别报告。')
    commands = parser.add_subparsers(dest='command', required=True)
    status = commands.add_parser('status', help='只读检查现有归档并显示科学结论')
    status.add_argument('--json', action='store_true')
    status.add_argument('--run-dir', type=Path, help='检查统一入口产生的已完成运行目录')
    report = commands.add_parser('report', help='从已核验归档生成统一总览，不重新拟合')
    report.add_argument('--output', type=Path, required=True, help='尚不存在的新目录')
    report.add_argument('--json', action='store_true')
    report.add_argument('--run-dir', type=Path, help='从指定的已完成运行目录生成总览')
    run = commands.add_parser('run', help='计算、独立核验并展示指定模块')
    run.add_argument('--module', choices=['all', *MODULES], default='all')
    run.add_argument('--output', type=Path, required=True, help='尚不存在的新目录')
    run.add_argument('--python', help='用于计算的Python可执行文件；默认当前解释器')
    run.add_argument('--download', action='store_true', help='允许软测量入口下载缺失的冻结工作簿')
    run.add_argument('--json', action='store_true')
    args = parser.parse_args(argv)
    try:
        if args.command == 'run':
            selected = MODULES if args.module == 'all' else (args.module,)
            data = run_project(args.output, selected, resolve_interpreter(args.python), args.download)
        else:
            data = collect_run(args.run_dir) if args.run_dir else collect()
            if args.command == 'report':
                output = check_output(args.output)
                if not data['overall_integrity_passed']:
                    raise ValueError('归档检查未通过：' + json.dumps(data, ensure_ascii=False))
                output.mkdir(parents=True, exist_ok=False)
                save_report(data, output)
        if args.json:
            print(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False))
        else:
            print_status(data)
            if args.command != 'status':
                print('\n项目总览：' + str(args.output.resolve() / 'index.html'))
        return 0 if data['overall_integrity_passed'] else 1
    except (OSError, ValueError, RuntimeError, subprocess.SubprocessError) as exc:
        print('未完成：' + str(exc), file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
