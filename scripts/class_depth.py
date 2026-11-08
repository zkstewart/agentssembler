#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Sequence-class depth: for each 10 kb window of a unitig, its primary ONT depth plus that of every
near-identical copy (self-alignment >= minId over >= minLen), read at the matching position.
Reads from identical copies split between them, so class depth / haplotype depth = haplotypes
carrying the sequence, and class depth vs. own depth shows how many copies share it.
Usage: class_depth.py depth.windows.tsv unitigs_self.paf.gz out.tsv [--minId 0.998 --minLen 5000]'''
import sys, gzip, argparse, collections
p=argparse.ArgumentParser(); p.add_argument("win"); p.add_argument("paf"); p.add_argument("out")
p.add_argument("--minId",type=float,default=0.998); p.add_argument("--minLen",type=int,default=5000); p.add_argument("--w",type=int,default=10000); a=p.parse_args()
D=collections.defaultdict(list)
for l in open(a.win):
    if l[0]=='#': continue
    s=l.split('\t'); D[s[0]].append((int(s[1]),float(s[3]),float(s[4]),float(s[5])))
A=collections.defaultdict(list)
with gzip.open(a.paf,'rt') as fh:
    for l in fh:
        s=l.split('\t',12)
        q,qs,qe,st,t,ts,te,m,b=s[0],int(s[2]),int(s[3]),s[4],s[5],int(s[7]),int(s[8]),int(s[9]),int(s[10])
        if q==t or b<a.minLen or m/b<a.minId: continue
        A[q].append((qs,qe,st,t,ts,te)); A[t].append((ts,te,st,q,qs,qe))
def dat(u,x):
    v=D.get(u)
    if not v: return 0.0
    i=min(len(v)-1,max(0,x//a.w)); return v[i][1]
with open(a.out,'w') as o:
    o.write("#unitig\tstart\tprimary\tunique\tsecondary\tn_copies\tclass_depth\tcopies\n")
    for u,v in D.items():
        for st,d,q,sec in v:
            mid=st+a.w//2; best={}
            for s0,e0,sd,t,ts,te in A.get(u,()):
                if s0<=mid<e0:
                    f=(mid-s0)/max(1,e0-s0); x=ts+f*(te-ts) if sd=='+' else te-f*(te-ts)
                    dv=dat(t,int(x)); best[t]=max(best.get(t,0),dv)
            o.write(f"{u}\t{st}\t{d:.1f}\t{q:.1f}\t{sec:.1f}\t{len(best)+1}\t{d+sum(best.values()):.1f}\t{','.join(best) or '.'}\n")
