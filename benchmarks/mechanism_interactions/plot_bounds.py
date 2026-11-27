"""Plot the held-out paired errors beside the retrospective within-assay limits."""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.ticker import PercentFormatter

HERE = Path(__file__).resolve().parent


def rows(source):
    with source.open(newline='') as stream:
        return list(csv.DictReader(stream))


def plot(output):
    font_manager.findfont('Noto Sans CJK JP', fallback_to_default=False)
    plt.rcParams.update({'font.family': 'Noto Sans CJK JP', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig, (pose, rank) = plt.subplots(1, 2, figsize=(11.6, 5.2), layout='constrained',
                                    gridspec_kw={'width_ratios': [1, 1.3]})
    transfer = rows(HERE.parent / 'ch_pi_search/gla_transfer/results/pairs.csv')
    for i, point in enumerate(transfer):
        delta = float(point['both_minus_baseline'])
        color = '#156963' if delta < 0 else '#ac451c'
        pose.hlines(i, min(delta, 0), max(delta, 0), color=color, lw=3)
        pose.scatter(delta, i, color=color, zorder=3, s=35)
        label = f'{delta:+.3f}' if abs(delta) < .01 else f'{delta:+.2f}'
        pose.annotate(label, (delta, i), xytext=(0, 8),
                      textcoords='offset points', ha='center', fontsize=9)
    pose.axvline(0, color='#657079', lw=.8)
    pose.set(yticks=range(len(transfer)), yticklabels=[f"{r['receptor']} · {r['seed']}" for r in transfer],
             ylim=(len(transfer)-.55,-.65), xlim=(-14,11), xlabel='联合方案减原生首姿 RMSD（Å），负值表示改善')
    pose.set_title('GLA 同一配体的六个配对条件', loc='left', pad=16)

    assays = [r for r in rows(HERE / 'shared_population/results/summary.csv') if r['reference'] == 'pool']
    labels = []
    styles = [('equal_credit','两状态等权','#243e57','o'), ('shared_ceiling','共同占比上限','#13786f','s'),
              ('separate_ceiling','逐对占比上限','#c87928','|')]
    for i, assay in enumerate(assays):
        name = assay['target'] if assay['target'] != 'BACE1' else 'BACE1 · '+assay['assay'].removeprefix('CHEMBL')
        labels.append(f"{name}（{assay['pairs']} 对）")
        scores = [float(assay[key])/float(assay['pairs']) for key,_,_,_ in styles]
        rank.hlines(i, min(scores),max(scores), color='#b5bfc5', lw=1.4)
        for index, (x, (_, label, color, marker)) in enumerate(zip(scores, styles)):
            rank.plot(x,i+.14*(index-1), marker=marker, ls='none', ms=8 if marker=='|' else 5.7,
                      color=color, label=label if i==0 else None, mew=1.65 if marker=='|' else .6)
    rank.set(yticks=range(len(assays)), yticklabels=labels, ylim=(len(assays)-.5,-.6), xlim=(-.03,1.04),
             xlabel='同一测定内 ≥3 倍分离配体对的排序一致率')
    rank.xaxis.set_major_formatter(PercentFormatter(1, decimals=0))
    rank.set_title('两份受体结构可解释的排序范围', loc='left', pad=16)
    rank.legend(loc='upper center',bbox_to_anchor=(.5,-.15), ncol=3, frameon=False, fontsize=9,
                handletextpad=.2, columnspacing=.9)
    for ax in (pose,rank):
        ax.grid(axis='x',color='#23323e', alpha=.1)
        ax.set_axisbelow(True)
        ax.spines['left'].set_visible(False)
        ax.tick_params(axis='y',length=0)
    fig.savefig(output, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output', type=Path)
    plot(parser.parse_args().output)
