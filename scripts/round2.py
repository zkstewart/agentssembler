#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Round 2 inputs. From a chain self-alignment and chain read depth:
  - redundant chains: >= 90% of windows covered by a longer chain (>= 99.8% identity) with summed depth
    below 1.5 x haplotype depth; high-depth (class >= --repeat x) chains are listed separately;
  - extra links: exact dovetails between chain ends (>= 99.8%, >= 5 kb, within --endSlack bp of both ends), as
    GFA L lines between the terminal unitigs, so build_chains treats them like graph links.
Usage: round2.py chains.agp chains.fasta.fai depth.windows.tsv chains_self.paf.gz kept.txt out_dir --hap H --repeat R'''
import sys, gzip, collections, argparse
import numpy as np
p=argparse.ArgumentParser(); [p.add_argument(x) for x in 'agp fai win paf kept out'.split()]
p.add_argument('--hap',type=float,required=True); p.add_argument('--endSlack',type=int,default=200); p.add_argument('--repeat',type=float,required=True); a=p.parse_args()
L={l.split('\t')[0]:int(l.split('\t')[1]) for l in open(a.fai)}
comp=collections.defaultdict(list)
for l in open(a.agp):
    if l[0]=='#': continue
    s=l.rstrip().split('\t'); comp[s[0]].append((s[5],s[8],int(s[6]),int(s[7])))
W=collections.defaultdict(list)
for l in open(a.win):
    if l[0]=='#': continue
    s=l.split('\t'); W[s[0]].append((int(s[1]),float(s[3])))
A=collections.defaultdict(list); dt=collections.defaultdict(list)
with gzip.open(a.paf,'rt') as fh:
    for l in fh:
        s=l.split('\t',12)
        q,ql,qs,qe,sd,t,tl,ts,te,m,b=s[0],int(s[1]),int(s[2]),int(s[3]),s[4],s[5],int(s[6]),int(s[7]),int(s[8]),int(s[9]),int(s[10])
        if q==t or b<5000 or m/b<0.998: continue
        A[q].append((qs,qe,t,ts,te,sd)); A[t].append((ts,te,q,qs,qe,sd))
        S=a.endSlack
        if sd=='+':
            if ql-qe<=S and ts<=S: e1,e2=(q,'+'),(t,'-')
            elif qs<=S and tl-te<=S: e1,e2=(q,'-'),(t,'+')
            else: continue
        else:
            if ql-qe<=S and tl-te<=S: e1,e2=(q,'+'),(t,'+')
            elif qs<=S and ts<=S: e1,e2=(q,'-'),(t,'-')
            else: continue
        ov=max(qe-qs,te-ts)
        if ov>=min(ql,tl)*0.9: continue # containment, not a dovetail
        dt[frozenset((e1,e2))].append(ov)
def wdepth(c,x):
    v=W.get(c)
    if not v: return 0
    return v[min(len(v)-1,max(0,x//10000))][1]
red=[];rep=[];kept2=[]
for c,w in W.items():
    n=len(w); cov=0; cls=[]
    for st,d in w:
        m=st+5000; best={}
        for a0,e0,t,ts,te,sd in A.get(c,()):
            if a0<=m<e0 and L[t]>L[c]:
                f=(m-a0)/max(1,e0-a0); x=ts+f*(te-ts) if sd=='+' else te-f*(te-ts)
                best[t]=max(best.get(t,0),wdepth(t,int(x)))
        if best: cov+=1
        cls.append(d+sum(best.values()))
    med=float(np.median(cls))
    if med>=a.repeat: rep.append(c)
    elif cov/n>=0.9 and med<1.5*a.hap: red.append(c)
drop=set(c for c in red)|set(rep)
# terminal unitig ends; a chain's + end is its last component's leaving end
def term(c,e):
    x=comp[c][-1] if e=='+' else comp[c][0]
    if e=='+': return x[0],x[1]
    return x[0],{'+':'-','-':'+'}[x[1]]
inv={'+':'-','-':'+'}
with open(a.out+'/dovetails.gfa','w') as o:
    n=0
    for k,v in dt.items():
        e1,e2=tuple(k)
        if e1[0] in drop or e2[0] in drop: continue
        (u1,s1),(u2,s2)=term(*e1),term(*e2)
        # L: leave u1 through s1, enter u2 oriented so the end met is s2
        o.write(f"L\t{u1}\t{s1}\t{u2}\t{inv[s2]}\t{max(v)}M\n"); n+=1
units=lambda c:[x[0] for x in comp[c]]
keep=[u for c in W if c not in drop for u in units(c)]
open(a.out+'/kept2.txt','w').write('\n'.join(keep)+'\n')
open(a.out+'/redundant_chains.txt','w').write(''.join(f"{c}\t{L[c]}\t{','.join(units(c))}\n" for c in red))
open(a.out+'/repeat_chains.txt','w').write(''.join(f"{c}\t{L[c]}\t{','.join(units(c))}\n" for c in rep))
print(f'redundant chains {len(red)} {sum(L[c] for c in red)/1e6:.1f} Mb; repeat chains {len(rep)} {sum(L[c] for c in rep)/1e6:.1f} Mb; kept unitigs {len(keep)}; dovetail links {n}')
