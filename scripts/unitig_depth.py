#!/usr/bin/env python3
'''Per-unitig ONT depth from a PAF: primary (all MAPQ), unique (primary MAPQ >= 20) and secondary.
Usage: unitig_depth.py <paf.gz> <unitigs.fai> <out_prefix> [--window 10000]
Writes <out>.windows.tsv (unitig start end total unique secondary) and <out>.unitigs.tsv (medians).'''
import sys, gzip, argparse
import numpy as np
p=argparse.ArgumentParser(); p.add_argument("paf"); p.add_argument("fai"); p.add_argument("out"); p.add_argument("--window",type=int,default=10000); a=p.parse_args()
ln={l.split("\t")[0]:int(l.split("\t")[1]) for l in open(a.fai)}
D={u:[np.zeros(n+1,dtype=np.int32) for _ in range(3)] for u,n in ln.items()}
with gzip.open(a.paf,"rt") as fh:
    for l in fh:
        s=l.split("\t",12); t=s[5]; x,y=int(s[7]),int(s[8]); prim="tp:A:P" in s[12]; q=int(s[11])
        d=D[t]
        if prim:
            d[0][x]+=1; d[0][y]-=1
            if q>=20: d[1][x]+=1; d[1][y]-=1
        else: d[2][x]+=1; d[2][y]-=1
w=open(a.out+".windows.tsv","w"); w.write("#unitig\tstart\tend\ttotal\tunique\tsecondary\n")
u_=open(a.out+".unitigs.tsv","w"); u_.write("#unitig\tlength\tdepth\tunique\tsecondary\tdepth_p10\tdepth_p90\n")
for u,n in ln.items():
    c=[np.cumsum(x[:-1]) for x in D[u]]
    W=[]
    for st in range(0,n,a.window):
        en=min(n,st+a.window)
        if en-st<a.window/2 and st>0: continue
        v=[float(x[st:en].mean()) for x in c]; W.append(v)
        w.write(f"{u}\t{st}\t{en}\t{v[0]:.1f}\t{v[1]:.1f}\t{v[2]:.1f}\n")
    W=np.array(W); t=W[:,0]
    u_.write(f"{u}\t{n}\t{np.median(t):.1f}\t{np.median(W[:,1]):.1f}\t{np.median(W[:,2]):.1f}\t{np.percentile(t,10):.1f}\t{np.percentile(t,90):.1f}\n")
