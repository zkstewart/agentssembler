#!/usr/bin/env python3
'''How much of each chain lies in a longer chain, summarised by depth. Used to look before setting the round-2 rules.
Usage: chain_dups.py chains.fasta.fai depth.windows.tsv chains_self.paf.gz'''
import sys, gzip, collections, numpy as np
fai, wins, selfpaf = sys.argv[1:4]
L={l.split('\t')[0]:int(l.split('\t')[1]) for l in open(fai)}
win=collections.defaultdict(list)
for l in open(wins):
    if l[0]=='#': continue
    s=l.split('\t'); win[s[0]].append((int(s[1]),float(s[3]),float(s[4])))
A=collections.defaultdict(list); cnt=collections.Counter(); ex=[]
with gzip.open(selfpaf,'rt') as fh:
    for l in fh:
        s=l.split('\t',12)
        q,ql,qs,qe,t,tl,ts,te,m,b=s[0],int(s[1]),int(s[2]),int(s[3]),s[5],int(s[6]),int(s[7]),int(s[8]),int(s[9]),int(s[10])
        if q==t or b<5000 or m/b<0.998: continue
        A[q].append((qs,qe,t,ts,te)); A[t].append((ts,te,q,qs,qe))
        qend=min(qs,ql-qe)<=2000; tend=min(ts,tl-te)<=2000
        c='dovetail' if qend and tend else 'one_end' if qend or tend else 'internal'; cnt[c]+=b
print({k:round(v/1e6,1) for k,v in cnt.items()})
rows=[]
for c,w in win.items():
    n=dup=0; d=[];dn=[]
    for st,p,q in w:
        m=st+5000; n+=1; d.append(p)
        if any(a<=m<b and L[t]>L[c] for a,b,t,_,_ in A.get(c,())): dup+=1; dn.append(p)
    rows.append((c,L[c],dup/max(n,1),float(np.median(d))))
rows.sort(key=lambda r:-r[1])
print('chains >=90% in a longer chain:',sum(r[2]>=.9 for r in rows),round(sum(r[1] for r in rows if r[2]>=.9)/1e6,1),'Mb')
print('chains 50-90%:',sum(.5<=r[2]<.9 for r in rows),round(sum(r[1] for r in rows if .5<=r[2]<.9)/1e6,1),'Mb')
print('chains <50%:',sum(r[2]<.5 for r in rows),round(sum(r[1] for r in rows if r[2]<.5)/1e6,1),'Mb')
for a,b in [(0,12),(12,18),(18,24),(24,36),(36,60),(60,1e9)]:
    x=[r for r in rows if r[2]>=.9 and a<=r[3]<b]; y=[r for r in rows if r[2]<.5 and a<=r[3]<b]
    print(f'  depth {a}-{b}: contained {len(x)} chains {sum(r[1] for r in x)/1e6:.1f} Mb | not {len(y)} chains {sum(r[1] for r in y)/1e6:.1f} Mb')
print('size dist chains <50% dup:',[ (int(t),sum(1 for r in rows if r[2]<.5 and r[1]>=t)) for t in (1e4,1e5,5e5,1e6,5e6,1e7)])
