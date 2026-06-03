#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Two-centre CRT mechanism for the cousin/sexy omega-distortions (Part IX open
problem). Extends the Part VII closed-form survival to pairs straddling two
consecutive centres (N, N+1). For each prime q>3 and each pair geometry, the CRT
factor is:
    q | N : deterministic -- each member 6(N+off)+s is q-safe unless off ≡ -s*6^{-1} (mod q)
    q ! N : averaged over admissible nonzero residues r of N
The pair-rate model is K * prod_q f_q, evaluated per centre by its true small-prime
divisibility, averaged within each left-centre omega stratum, normalised to omega=1.
No fitted parameters; the only inputs are each pair's geometry:
    twin   = ((0,-1),(0,+1))   single centre, both wings
    cousin = ((0,+1),(1,-1))   N right wing + (N+1) left wing
    sexyA  = ((0,-1),(1,-1))   N left wing + (N+1) left wing
    sexyB  = ((0,+1),(1,+1))   N right wing + (N+1) right wing
Memory-light streaming. Default S10. Requires: numpy.
"""
# CRT factorisation, extended to two centres (N, N+1), reproduces the measured
# omega-shapes of twin/cousin/sexyA/sexyB.
import numpy as np, math, os
def primes_upto(n):
    s=np.ones(n+1,bool); s[:2]=False
    for i in range(2,int(math.isqrt(n))+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0].astype(np.int64)
MAXK=int(os.environ.get("MAXK",10))
LO=10**(MAXK-1)//6+1; HI=10**MAXK//6; SEG=4_000_000
PB=int(math.isqrt(6*HI+250))+1; BP=primes_upto(PB)
POOL=[5,7,11,13,17,19,23,29,31,37,41,43,47]
poolbit={q:1<<i for i,q in enumerate(POOL)}
def dead_wing(q,s):
    inv=pow(6,q-2,q); return (-s*inv)%q
PAIRS={'twin':((0,-1),(0,+1)),'cousin':((0,+1),(1,-1)),'sexyA':((0,-1),(1,-1)),'sexyB':((0,+1),(1,+1))}
def pair_factor(q, offw, qdiv):
    if qdiv:
        for off,s in offw:
            if off%q==dead_wing(q,s)%q: return 0.0
        return 1.0
    cnt=tot=0
    for r in range(q):
        if r==0: continue
        tot+=1
        if all((r+off)%q!=dead_wing(q,s) for off,s in offw): cnt+=1
    return cnt/tot
fac={}
for name,offw in PAIRS.items():
    fac[name]=({q:pair_factor(q,offw,True) for q in POOL},{q:pair_factor(q,offw,False) for q in POOL})
OMAX=7
# accumulators: per omega -> [n_centers, sum_meas_pair, sum_model] for each pair
acc={name:{'tot':np.zeros(OMAX+2),'meas':np.zeros(OMAX+2),'mod':np.zeros(OMAX+2)} for name in PAIRS}
prev_pm=prev_pp=None; prev_om=None; prev_lastN=None
n=LO
while n<=HI:
    nh=min(n+SEG,HI+1); sz=nh-n
    rem=np.arange(n,nh,dtype=np.int64); ob=np.zeros(sz,np.int16); mk=np.zeros(sz,np.int32)
    for p in BP:
        if p*p>nh-1: break
        f=((n+p-1)//p)*p
        if f>=nh: continue
        idx=np.arange(f-n,sz,p)
        if idx.size==0: continue
        sub=rem[idx]; m=(sub%p)==0
        while m.any(): sub[m]//=p; m=(sub%p)==0
        rem[idx]=sub
        if p>3:
            ob[idx]+=1
            if p in poolbit: mk[idx]|=poolbit[p]
    ob[rem>1]+=1
    Narr=np.arange(n,nh,dtype=np.int64)
    vlo=6*n-1; vhi=6*(nh-1)+1; span=vhi-vlo+1
    comp=np.zeros(span,bool); sq=int(math.isqrt(vhi))+1
    for p in BP:
        if p>sq: break
        st=max(p*p,((vlo+p-1)//p)*p)
        if st>vhi: continue
        comp[st-vlo:span:p]=True
    pm=~comp[(6*Narr-1)-vlo]; pp=~comp[(6*Narr+1)-vlo]
    omc=np.clip(ob,0,OMAX+1)
    # model factor arrays per pair
    modarr={}
    for name in PAIRS:
        fqN,fqn=fac[name]; out=np.ones(sz)
        for q in POOL:
            has=(mk&poolbit[q])>0
            out*=np.where(has,fqN[q],fqn[q])
        modarr[name]=out
    # twin: single center
    tw=pm&pp
    for om in range(1,OMAX+1):
        s=(omc==om); ns=s.sum()
        if ns==0: continue
        acc['twin']['tot'][om]+=ns
        acc['twin']['meas'][om]+=(tw&s).sum()
        acc['twin']['mod'][om]+=modarr['twin'][s].sum()
    # straddling within segment
    if sz>=2:
        leftom=omc[:-1]
        pairmask={'cousin':pp[:-1]&pm[1:],'sexyA':pm[:-1]&pm[1:],'sexyB':pp[:-1]&pp[1:]}
        for name in ['cousin','sexyA','sexyB']:
            for om in range(1,OMAX+1):
                s=(leftom==om); ns=s.sum()
                if ns==0: continue
                acc[name]['tot'][om]+=ns
                acc[name]['meas'][om]+=pairmask[name][s].sum()
                acc[name]['mod'][om]+=modarr[name][:-1][s].sum()
    n=nh
print(f"S{MAXK} done")
def shape(name):
    a=acc[name]; oms=[om for om in range(1,OMAX+1) if a['tot'][om]>=20000]
    meas=np.array([a['meas'][om]/a['tot'][om] for om in oms])
    mod=np.array([a['mod'][om]/a['tot'][om] for om in oms])
    return oms,meas/meas[0],mod/mod[0]
print(f"\nMeasured vs MODEL normalised shapes (omega=1 base):")
for name in ['twin','cousin','sexyA','sexyB']:
    oms,m,p=shape(name)
    print(f"  {name:>7}:")
    print(f"    omega:    {oms}")
    print(f"    measured: {[f'{x:.3f}' for x in m]}")
    print(f"    model:    {[f'{x:.3f}' for x in p]}")
    err=100*np.abs(p-m)/m
    print(f"    max err:  {err.max():.1f}%")

# ---- emit CSV (cs_mechanism_S{K}_data.csv) ----
import csv as _csv
with open(f'cs_mechanism_S{MAXK}_data.csv','w',newline='') as _f:
    _w=_csv.writer(_f); _w.writerow(['pair','omega','measured','model','err_pct'])
    for name in ['twin','cousin','sexyA','sexyB']:
        oms,m,p=shape(name)
        for o,mm,pp in zip(oms,m,p):
            e=100*(pp-mm)/mm if mm>0 else 0
            _w.writerow([name,o,f'{mm:.3f}',f'{pp:.3f}',f'{e:.1f}'])
print(f"\n[ok] wrote cs_mechanism_S{MAXK}_data.csv")
