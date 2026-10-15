import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rows(path):
    with Path(path).open(encoding='utf-8',newline='') as stream:
        return list(csv.DictReader(stream))


def plot(output):
    report = json.loads((output/'report.json').read_text(encoding='utf-8'))
    verified = json.loads((output/'verification.json').read_text(encoding='utf-8'))
    assert verified['verified'] and verified['report_sha256']==sha(output/'report.json')
    for name,expected in report['artifact_sha256'].items():
        assert sha(output/name)==expected
    pairs = [r for r in rows(output/'paired_batches.csv') if r['direction']=='1_to_2' and r['target']=='CDW_g_L']
    summary = rows(output/'summary.csv')
    names = ['training_mean','time_recipe_ridge','instantaneous_ridge','history_ridge']
    colors = ['#8b94a3','#376b9d','#229b8d','#cf783a']
    fig,axes = plt.subplots(1,2,figsize=(13.5,6),gridspec_kw={'width_ratios':[1,1.4]})
    ax = axes[0]
    for offset,direction,label in [(-.18,'1_to_2','Primary: round 1 to 2'),(.18,'2_to_1','Reverse sensitivity: round 2 to 1')]:
        values = [float(next(r for r in summary if r['direction']==direction and r['target']=='CDW_g_L' and r['strategy']==n)['batch_macro_RMSE']) for n in names]
        ax.barh(np.arange(4)+offset,values,height=.33,color=colors,alpha=1 if offset<0 else .45,label=label)
        for y,value in zip(np.arange(4)+offset,values):
            ax.text(value+.055,y,f'{value:.3f}',va='center',fontsize=9)
    ax.set_yticks(range(4),['Training mean','Time + recipe','Instantaneous signals','Signals + causal history'])
    ax.invert_yaxis()
    ax.set_xlim(0,6.8)
    ax.set_xlabel('Mean culture RMSE (CDW, g/L)')
    ax.set_title('Fixed models; whole-round holdout')
    handles,labels = ax.get_legend_handles_labels()
    ax = axes[1]
    deltas = [float(r['delta_RMSE']) for r in pairs]
    ax.barh(range(len(pairs)),deltas,color=['#229b8d' if d<0 else '#cf783a' for d in deltas])
    ax.set_yticks(range(len(pairs)),[r['batch_id'].removesuffix(' (2)') for r in pairs])
    ax.invert_yaxis()
    ax.axvline(0,color='#344357',lw=1)
    ax.set_xlabel('History minus instantaneous RMSE (g/L)')
    ax.set_title('Primary comparison: 2/10 cultures improve')
    ax.text(.98,.03,f'Mean delta: +{np.mean(deltas):.3f} g/L',ha='right',transform=ax.transAxes,fontsize=10)
    for ax in axes:
        ax.spines[['top','right']].set_visible(False)
        ax.grid(axis='x',alpha=.15)
        ax.set_axisbelow(True)
    fig.suptitle('Tamminen 2022 public fungal fermentation data: biomass replay',fontsize=14)
    fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.5,.93),ncol=2,fontsize=9,frameon=False)
    fig.text(.5,.02,'Ten matched conditions within two correlated rounds; descriptive results, no significance claim.',ha='center',fontsize=9,color='#555555')
    fig.tight_layout(rect=[0,.05,1,.88],w_pad=3)
    path = output/'biomass_round_holdout.png'
    fig.savefig(path,dpi=180)
    plt.close(fig)
    provenance = {'created_at_utc':datetime.now(timezone.utc).isoformat(),'plotter_sha256':sha(__file__),
                  'report_sha256':sha(output/'report.json'),'verification_sha256':sha(output/'verification.json'),
                  'input_sha256':{n:sha(output/n) for n in ['paired_batches.csv','summary.csv']},
                  'output_sha256':{path.name:sha(path)},'matplotlib_version':matplotlib.__version__,
                  'no_scientific_refitting':True}
    (output/'visualization.json').write_text(json.dumps(provenance,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'figure':str(path),'sha256':sha(path)}))


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=ROOT/'research/results/ambr_dynamic')
    args = parser.parse_args()
    plot(args.output)
