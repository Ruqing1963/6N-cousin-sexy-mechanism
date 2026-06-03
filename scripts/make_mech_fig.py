#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build the 2-panel mechanism figure from ../data/cs_mechanism_S10_data.csv
(produced by cs_mechanism.py with default MAXK=10). Left: two-centre CRT model
(lines) vs measured (markers) for all four constellations. Right: residuals
(<=3% for omega<=6; omega=7 small-sample).
"""
import csv, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
rows=list(csv.DictReader(open('../data/cs_mechanism_S10_data.csv')))
def series(pair):
    om=[];m=[];mo=[]
    for r in rows:
        if r['pair']==pair:
            om.append(int(r['omega'])); m.append(float(r['measured'])); mo.append(float(r['model']))
    return map(np.array,(om,m,mo))
fig,axes=plt.subplots(1,2,figsize=(14,5.4))
cfg=[('twin','#c0392b','o','twin (diff 2)'),('cousin','#185FA5','s','cousin (diff 4)'),
     ('sexyA','#3aa0d0','^','sexy-A (diff 6, L wings)'),('sexyB','#2ca25f','D','sexy-B (diff 6, R wings)')]
for pair,col,mk,lab in cfg:
    om,m,mo=series(pair)
    axes[0].plot(om,m,mk,color=col,ms=8,zorder=4)
    axes[0].plot(om,mo,'-',color=col,lw=1.8,alpha=.8,label=lab,zorder=3)
axes[0].axhline(1,color='gray',ls=':',lw=1)
axes[0].set_xlabel(r'$\omega_{>3}$ of left centre $N$',fontsize=11)
axes[0].set_ylabel(r'pair rate, normalised to $\omega=1$',fontsize=11)
axes[0].set_title('Two-centre CRT model (lines) vs measured (markers)',fontsize=12)
axes[0].legend(fontsize=9,loc='upper left'); axes[0].grid(alpha=.25); axes[0].set_xticks(range(1,8))
axes[1].axhline(0,color='gray',lw=1)
axes[1].axhspan(-3,3,color='#2ca25f',alpha=.10,label='$\\pm3\\%$ band')
for pair,col,mk,lab in cfg:
    om,m,mo=series(pair)
    err=100*(mo-m)/m
    axes[1].plot(om,err,mk+'-',color=col,lw=1.6,ms=7,label=pair)
axes[1].set_xlabel(r'$\omega_{>3}$ of left centre $N$',fontsize=11)
axes[1].set_ylabel('model $-$ measured  (%)',fontsize=11)
axes[1].set_title('residuals: $\\leq2.3\\%$ for $\\omega\\leq6$;\n$\\omega=7$ small-sample (sexy-A $0.9\\sigma$)',fontsize=12)
axes[1].legend(fontsize=9); axes[1].grid(alpha=.25); axes[1].set_xticks(range(1,8)); axes[1].set_ylim(-16,8)
plt.suptitle('Two-centre CRT mechanism reproduces all three $\\omega$-distortions in $S_{10}$ (no fitted parameters)',fontsize=12.5,y=1.02)
plt.tight_layout()
plt.savefig('fig_paper10_mechanism.pdf',bbox_inches='tight')
plt.savefig('fig_paper10_mechanism.png',dpi=160,bbox_inches='tight')
print("figure saved")
