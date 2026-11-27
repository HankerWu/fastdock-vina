"""Read the first emitted AD4 score from each raw pool and independently recount Table II."""
import argparse
import csv
from itertools import permutations
from pathlib import Path

import numpy as np
from scipy.stats import kendalltau, spearmanr

# Holt et al., JACS (1993), doi:10.1021/ja00075a008, Table II, p. 9927.
KI_NM = {4:(186,186),5:(110,110),6:(12,12),7:(250,250),
         8:(10,10),9:(7,7),10:(300,600),11:(300,300)}


def first_score(path):
    with path.open() as f:
        return float(next(line for line in f if line.startswith('REMARK VINA RESULT:')).split()[3])


def check(work, main, out):
    with main.open() as f:
        summaries = list(csv.DictReader(f))
    pairs = [(i,j) for i,j in permutations(KI_NM,2) if KI_NM[j][0]>=3*KI_NM[i][1]]
    assert len(pairs)==15
    point_ids = [i for i,(lo,hi) in KI_NM.items() if lo==hi]
    actual = [KI_NM[i][0] for i in point_ids]
    by_arm = {}
    for row in summaries:
        arm, receiver = row['arm'],row['receiver']
        recs = ('1FKG','1FKH') if receiver=='mean' else (receiver,)
        values = {i:np.mean([first_score(work/'pools'/arm/f'{rec}_{i}.pdbqt') for rec in recs]) for i in KI_NM}
        predicted = [values[i] for i in point_ids]
        margin = np.array([values[j]-values[i] for i,j in pairs])
        credit = np.where(margin>0,1,np.where(margin==0,.5,0))
        metrics = [sum(credit),spearmanr(predicted,actual).statistic,kendalltau(predicted,actual).statistic,
                   values[10]-values[8],values[11]-values[9]]
        np.testing.assert_allclose(metrics,[float(row[k]) for k in ('correct','rho','tau_b','cliff8_10','cliff9_11')],
                                   atol=1e-10,rtol=0)
        by_arm[arm,receiver] = credit
    reports=[]
    for (arm,receiver),credit in by_arm.items():
        cov = by_arm['coverage',receiver]
        def labels(mask):
            return ' '.join(f'{i}/{j}' for (i,j),ok in zip(pairs,mask,strict=True) if ok)
        reports.append(dict(receiver=receiver,compared_to=arm,coverage_correct=sum(cov),other_correct=sum(credit),
            coverage_only=labels(cov>credit),other_only=labels(credit>cov),
            both_correct=labels((credit==1)&(cov==1)),both_wrong=labels((credit==0)&(cov==0))))
    with out.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(reports[0]))
        w.writeheader()
        w.writerows(reports)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for n in ('work','assay_csv','out'):
        p.add_argument(n,type=Path)
    a=p.parse_args()
    check(a.work,a.assay_csv,a.out)
