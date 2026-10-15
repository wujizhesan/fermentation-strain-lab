import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent


def rows(path):
    with path.open(encoding='utf-8', newline='') as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=ROOT / 'research/results/ale_lineage_screen_v1')
    args = parser.parse_args()
    folder = args.output
    sha = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    report = json.loads((folder / 'report.json').read_text(encoding='utf-8'))
    verification = json.loads((folder / 'verification.json').read_text(encoding='utf-8'))
    if not verification['passed'] or verification['report_sha256'] != sha(folder / 'report.json'):
        raise ValueError('Current independently verified results required')
    for name in ['lineage_metrics.csv', 'candidate_decisions.csv', 'features.csv']:
        if sha(folder / name) != report['artifact_sha256'][name]:
            raise ValueError('Plot input fingerprint changed')
    for name in ['selection_comparison.png', 'selection_comparison.svg', 'visualization.json']:
        if (folder / name).exists():
            raise ValueError('Plot outputs already exist')
    metrics = rows(folder / 'lineage_metrics.csv')
    decisions = rows(folder / 'candidate_decisions.csv')
    features = {row['isolate_id']: row for row in rows(folder / 'features.csv')}
    groups = ['Chassis_F-1', 'Chassis_F-2', 'Chassis_F-3']
    methods = ['uniform_random_expectation', 'early_density', 'early_density_fold', 'fixed_ridge']
    labels = ['均匀随机期望', '早期OD', '早期OD倍数', '固定Ridge']
    colors = ['#8793a5', '#2456a4', '#b37d12', '#b24e53']
    plt.rcParams.update({'font.sans-serif': ['Microsoft YaHei', 'DejaVu Sans'], 'axes.unicode_minus': False,
        'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.fonttype': 'none'})
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.7), gridspec_kw={'width_ratios': [1.15, 1]})
    positions = np.arange(3)
    for index, (method, label, color) in enumerate(zip(methods, labels, colors)):
        values = [float(next(row for row in metrics if row['heldout_lineage'] == group and row['method'] == method)['fraction_of_heldout_best']) for group in groups]
        axes[0].bar(positions + (index - 1.5) * .19, values, width=.18, color=color, label=label)
    axes[0].set_xticks(positions, ['留出谱系1', '留出谱系2', '留出谱系3'])
    axes[0].set_ylim(0, 1.27)
    axes[0].set_ylabel('所选滴度 / 该留出谱系最高滴度')
    axes[0].set_title('A  同一候选池、同一选株预算', loc='left', fontsize=12)
    axes[0].legend(frameon=False, ncol=2, fontsize=9, loc='upper center')
    ridge = [row for row in decisions if row['method'] == 'fixed_ridge']
    maximum = max(float(row['candidate_score']) for row in ridge) * 1.08
    axes[1].plot([0, maximum], [0, maximum], '--', color='#8591a0', linewidth=1, label='预测 = 观测')
    for group, color, marker in zip(groups, ['#b24e53', '#2456a4', '#21766b'], ['o', 's', '^']):
        selected = [row for row in ridge if row['heldout_lineage'] == group]
        axes[1].scatter([float(features[row['isolate_id']]['titer_g_L']) for row in selected],
            [float(row['candidate_score']) for row in selected], s=48, color=color, marker=marker,
            edgecolors='white', linewidths=.5, label='谱系' + group[-1])
    axes[1].set_xlim(0, maximum)
    axes[1].set_ylim(0, maximum)
    axes[1].set_xlabel('公开终点滴度（g/L）')
    axes[1].set_ylabel('留出Ridge预测（g/L）')
    axes[1].set_title('B  9株的整谱系留出预测', loc='left', fontsize=12)
    axes[1].legend(frameon=False, fontsize=9)
    for axis in axes:
        axis.set_axisbelow(True)
        axis.grid(axis='y', color='#e5eaf0', linewidth=.7)
    fig.suptitle('公开酵母富马酸ALE：选株回放', x=.065, ha='left', fontsize=15, fontweight='bold')
    fig.text(.065, .015, '每折从留出谱系3株中选1株；随机值为期望。仅3个群体，数据已在此前查看；无新盲测、误差条或显著性主张。', fontsize=9, color='#586779')
    fig.tight_layout(rect=[0, .055, 1, .93])
    outputs = []
    for suffix in ['png', 'svg']:
        path = folder / ('selection_comparison.' + suffix)
        fig.savefig(path, dpi=180, facecolor='white')
        outputs.append({'path': path.name, 'sha256': sha(path)})
    plt.close(fig)
    (folder / 'visualization.json').write_text(json.dumps({'report_sha256': sha(folder / 'report.json'), 'outputs': outputs}, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(outputs, indent=2))


if __name__ == '__main__':
    main()
