import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from datetime import datetime, timezone

from prepare_marge_demo import prepare
from build_marge_viewer import build_viewer


ROOT = Path(__file__).resolve().parent


def run(arguments, environment):
    completed = subprocess.run([sys.executable, '-X', 'utf8', '-B', *arguments], cwd=ROOT, env=environment)
    if completed.returncode:
        raise RuntimeError('步骤未通过：' + ' '.join(arguments))


def compare_variants(destinations):
    original, corrected = [json.loads((path / 'report.json').read_text(encoding='utf-8')) for path in destinations]
    before, after = [{r['case_id']: r for r in report['cases']} for report in [original, corrected]]
    result = {'cases': len(before), 'point_solutions_may_be_nonunique': True,
        'verified_backend': 'ReFramed 1.2.1 / OptLang / GLPK', 'paper_cplex_behavior_verified': False,
        'stage1_max_absolute_difference': max(abs(before[key]['stage1_error_optimum'] - after[key]['stage1_error_optimum']) for key in before),
        'stage2_original_range': [min(r['stage2_proteome_optimum'] for r in before.values()), max(r['stage2_proteome_optimum'] for r in before.values())],
        'stage2_corrected_range': [min(r['stage2_proteome_optimum'] for r in after.values()), max(r['stage2_proteome_optimum'] for r in after.values())],
        'interval_widths': {}}
    for path in destinations:
        counts = {'narrower': 0, 'same': 0, 'wider': 0}
        examples = []
        for design in ['F', 'Sz', 'Mz']:
            for layer in ['protein', 'rna']:
                modes = {}
                for mode in ['released_ratio', 'neutral_ratio']:
                    with (path / f'{design}_{layer}_{mode}' / 'flux_ranges.csv').open(encoding='utf-8', newline='') as handle:
                        modes[mode] = {(r['condition'], r['reaction_id']): float(r['maximum']) - float(r['minimum']) for r in csv.DictReader(handle)}
                for key, width in modes['released_ratio'].items():
                    neutral = modes['neutral_ratio'][key]
                    kind = 'narrower' if width < neutral - 1e-7 else ('wider' if width > neutral + 1e-7 else 'same')
                    counts[kind] += 1
                    if kind == 'wider':
                        examples.append({'design': design, 'layer': layer, 'condition': key[0],
                            'reaction_id': key[1], 'neutral_width': neutral, 'released_width': width})
        result['interval_widths'][path.name] = {'counts': counts, 'wider_cases': examples}
    return result


def main():
    parser = argparse.ArgumentParser(description='准备公开输入、运行两套MARGE实现、独立核验并测试。')
    parser.add_argument('--output', type=Path, default=ROOT / 'research/results')
    args = parser.parse_args()
    base = args.output.resolve()
    destinations = [base / name for name in ['ale_marge_demo_v1', 'ale_marge_balancefix_v1']]
    for path in destinations:
        if path.exists() and any(path.iterdir()):
            raise ValueError('结果目录已有内容；请用--output指定新目录：' + str(path))
    if (base / 'run_receipt.json').exists():
        raise ValueError('运行回执已存在；请指定新目录')
    checks = prepare()
    environment = dict(os.environ)
    environment.update({'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1', 'MKL_NUM_THREADS': '1',
        'MARGE_RESULTS_ROOT': str(base)})
    for index, path in enumerate(destinations):
        flag = ['--keep-proteome-balance'] if index else []
        run(['ale_marge_demo.py', '--output', str(path), *flag], environment)
        run(['verify_ale_marge_demo.py', '--output', str(path), *flag], environment)
    run(['-m', 'unittest', 'discover', '-s', 'gem_tests', '-v'], environment)
    receipt = {'passed': True, 'completed_at_utc': datetime.now(timezone.utc).isoformat(),
        'input_validation': checks, 'python_version': sys.version,
        'variants': {}, 'paper_exact_reproduction_completed': False}
    for path in destinations:
        verification = json.loads((path / 'verification.json').read_text(encoding='utf-8'))
        digest = hashlib.sha256((path / 'report.json').read_bytes()).hexdigest()
        if not verification['passed'] or verification['report_sha256'] != digest:
            raise ValueError('核验回执与当前结果不一致')
        receipt['variants'][path.name] = {'report_sha256': digest,
            'joint_lps_resolved': verification['joint_lps_resolved'],
            'range_endpoints_resolved': verification['range_endpoints_resolved']}
    comparison = compare_variants(destinations)
    (base / 'comparison.json').write_text(json.dumps(comparison, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    receipt['comparison_sha256'] = hashlib.sha256((base / 'comparison.json').read_bytes()).hexdigest()
    viewer = build_viewer(base, base / 'marge-demo.html', 'run_result')
    receipt['viewer'] = {key: value for key, value in viewer.items() if key != 'path'}
    receipt['viewer']['path'] = 'marge-demo.html'
    (base / 'run_receipt.json').write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
    print('运行完成，回执：' + str(base / 'run_receipt.json'), flush=True)
    print('交互式报告：' + str(base / 'marge-demo.html'), flush=True)


if __name__ == '__main__':
    main()
