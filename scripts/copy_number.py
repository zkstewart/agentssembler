#!/usr/bin/env python3
'''Non-redundant unitig set from read depth. Unitigs are taken longest first. For each 10 kb window
the sequence class depth / haplotype depth gives k, the haplotypes carrying that sequence (>= 1).
A unitig is redundant when >= --minCov of its windows already have >= k copies among the kept
unitigs (near-identical, self-alignment >= --minId over >= --minLen). Otherwise it is kept: identical
sequence is kept as often as the depth says it occurs, and unique sequence is never dropped.
Usage: copy_number.py class_depth.tsv unitigs_self.paf.gz out_prefix --fai unitigs.fasta.fai --hap H --repeat R'''
import gzip, argparse, collections
import numpy as np
p=argparse.ArgumentParser(); p.add_argument("cls"); p.add_argument("paf"); p.add_argument("out")
p.add_argument("--fai",required=True); p.add_argument("--hap",type=float,required=True); p.add_argument("--minId",type=float,default=0.998); p.add_argument("--minLen",type=int,default=5000)
p.add_argument("--minCov",type=float,default=0.9); p.add_argument("--w",type=int,default=10000); p.add_argument("--repeat",type=float,required=True); a=p.parse_args()
W=collections.defaultdict(list)
for l in open(a.cls):
    if l[0]=='#': continue
    s=l.split('\t'); W[s[0]].append((int(s[1]),float(s[2]),float(s[3]),float(s[4]),float(s[6])))
Ln={l.split('\t')[0]:int(l.split('\t')[1]) for l in open(a.fai)}
A=collections.defaultdict(list)
with gzip.open(a.paf,'rt') as fh:
    for l in fh:
        s=l.split('\t',12)
        q,qs,qe,t,ts,te,m,b=s[0],int(s[2]),int(s[3]),s[5],int(s[7]),int(s[8]),int(s[9]),int(s[10])
        if q==t or b<a.minLen or m/b<a.minId: continue
        A[q].append((qs,qe,t)); A[t].append((ts,te,q))
kept=set(); rows=[]
for u in sorted(W,key=lambda x:-Ln[x]):
    ws=W[u]; n=red=0; ks=[]
    for st,pd,qd,sd,cd in ws:
        k=max(1,int(round(cd/a.hap))); ks.append(k); mid=st+a.w//2
        have={ t for s0,e0,t in A.get(u,()) if s0<=mid<e0 and t in kept }
        n+=1; red+= len(have)>=k
    frac=red/max(1,n)
    prim=np.median([x[1] for x in ws]); uq=np.median([x[2] for x in ws]); sec=np.median([x[3] for x in ws]); cd=np.median([x[4] for x in ws])
    high=cd>=a.repeat
    status='redundant' if frac>=a.minCov else 'keep'
    if status=='keep': kept.add(u)
    rows.append((u,Ln[u],prim,uq,sec,cd,int(np.median(ks)),frac,'repeat' if high else 'normal',status))
with open(a.out+'.tsv','w') as o:
    o.write('#unitig\tlength\tprimary\tunique\tsecondary\tclass_depth\tk\tfrac_covered_by_kept\ttype\tstatus\n')
    for r in rows: o.write(f"{r[0]}\t{r[1]}\t{r[2]:.1f}\t{r[3]:.1f}\t{r[4]:.1f}\t{r[5]:.1f}\t{r[6]}\t{r[7]:.3f}\t{r[8]}\t{r[9]}\n")
tot=lambda f:sum(r[1] for r in rows if f(r))/1e6
print(f"keep {sum(r[9]=='keep' for r in rows)} unitigs {tot(lambda r:r[9]=='keep'):.1f} Mb; redundant {sum(r[9]!='keep' for r in rows)} unitigs {tot(lambda r:r[9]!='keep'):.1f} Mb")
for k in (1,2,3,4): print(f" kept k={k}: {tot(lambda r:r[9]=='keep' and r[6]==k):.1f} Mb")
print(f" kept k>=5: {tot(lambda r:r[9]=='keep' and r[6]>=5):.1f} Mb")
kp=[r for r in rows if r[9]=='keep']
for a_,b_ in ((0,10),(10,18),(18,24),(24,32),(32,60),(60,1e9)):
    print(f" kept primary {a_}-{b_}: {sum(r[1] for r in kp if a_<=r[2]<b_)/1e6:.1f} Mb, class {sum(r[1] for r in kp if a_<=r[5]<b_)/1e6:.1f} Mb")
