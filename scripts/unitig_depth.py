#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Per-unitig ONT depth from a PAF: primary (all MAPQ), unique (primary MAPQ >= 20) and secondary.
Usage: unitig_depth.py <paf.gz> <unitigs.fai> <out_prefix> [--window 10000]
Writes <out>.windows.tsv (unitig start end total unique secondary) and <out>.unitigs.tsv (medians).
Each window's depth is the bp of alignment inside it over its length, which is the mean of per-base depth, accumulated per window: memory
grows with the number of windows, not of bases (per-base arrays took 12 bytes/bp, past RAM on a multi-Gb genome).'''
import sys, gzip, argparse
import numpy as np
p=argparse.ArgumentParser(); p.add_argument("paf"); p.add_argument("fai"); p.add_argument("out"); p.add_argument("--window",type=int,default=10000); a=p.parse_args()
W_=a.window
ln={l.split("\t")[0]:int(l.split("\t")[1]) for l in open(a.fai)}
# per unitig and kind (total, unique, secondary): S = bp in each window from alignment ends; F = difference array of windows fully covered
S={u:np.zeros((3,n//W_+2)) for u,n in ln.items()}; F={u:np.zeros((3,n//W_+3),dtype=np.int64) for u,n in ln.items()}
def add(u,k,x,y):
    if y<=x: return
    k0,k1=x//W_,(y-1)//W_; s=S[u][k]
    if k0==k1: s[k0]+=y-x; return
    s[k0]+=(k0+1)*W_-x; s[k1]+=y-k1*W_
    if k1>k0+1: f=F[u][k]; f[k0+1]+=1; f[k1]-=1
with gzip.open(a.paf,"rt") as fh:
    for l in fh:
        s=l.split("\t",12); t=s[5]; x,y=int(s[7]),int(s[8]); q=int(s[11])
        if "tp:A:P" in s[12]:
            add(t,0,x,y)
            if q>=20: add(t,1,x,y)
        else: add(t,2,x,y)
w=open(a.out+".windows.tsv","w"); w.write("#unitig\tstart\tend\ttotal\tunique\tsecondary\n")
u_=open(a.out+".unitigs.tsv","w"); u_.write("#unitig\tlength\tdepth\tunique\tsecondary\tdepth_p10\tdepth_p90\n")
for u,n in ln.items():
    bp=S[u]+np.cumsum(F[u],axis=1)[:,:S[u].shape[1]]*W_
    W=[]
    for st in range(0,n,W_):
        en=min(n,st+W_)
        if en-st<W_/2 and st>0: continue
        v=[float(bp[k][st//W_])/(en-st) for k in range(3)]; W.append(v)
        w.write(f"{u}\t{st}\t{en}\t{v[0]:.1f}\t{v[1]:.1f}\t{v[2]:.1f}\n")
    W=np.array(W); t=W[:,0]
    u_.write(f"{u}\t{n}\t{np.median(t):.1f}\t{np.median(W[:,1]):.1f}\t{np.median(W[:,2]):.1f}\t{np.percentile(t,10):.1f}\t{np.percentile(t,90):.1f}\n")
