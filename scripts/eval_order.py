#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Compare an order with the reference chromosome (a check, not an input): each chain gets a median reference
position on the group's main chromosome from unitig alignments; reports chains in order with position, and
the fraction of adjacent ordered pairs moving consistently along the reference.
Usage: eval_order.py order.tsv chains.agp unitigs_vs_ref.paf.gz'''
import sys, gzip, collections
ordf,agp,paf=sys.argv[1:4]
comp=collections.defaultdict(list)
for l in open(agp):
    if l[0]=='#': continue
    s=l.rstrip().split('\t'); comp[s[5]].append(s[0])
hits=collections.defaultdict(list)
with gzip.open(paf,'rt') as fh:
    for l in fh:
        s=l.split('\t',12)
        if 'tp:A:P' in s[12] and int(s[11])>=20 and s[0] in comp:
            for c in comp[s[0]]: hits[c].append((s[5],(int(s[7])+int(s[8]))//2,int(s[3])-int(s[2])))
G=collections.defaultdict(list)
for l in open(ordf):
    if l[0]=='#': continue
    s=l.rstrip().split('\t')
    if s[1]!='.': G[s[0]].append(s)
tot=good=0
for g,v in sorted(G.items()):
    cnt=collections.Counter()
    for s in v:
        for r,_,n in hits[s[2]]: cnt[r]+=n
    main=cnt.most_common(1)[0][0] if cnt else '.'
    seq=[]
    for s in v:
        h=[(p,n) for r,p,n in hits[s[2]] if r==main]
        if h:
            h.sort(); t=sum(n for _,n in h); c=0
            for p,n in h:
                c+=n
                if c>=t/2: seq.append((s[2],p/1e6,int(s[4])/1e6,s[6],s[7],s[8])); break
        else: seq.append((s[2],None,int(s[4])/1e6,s[6],s[7],s[8]))
    pos=[x[1] for x in seq if x[1] is not None]
    up=sum(1 for a,b in zip(pos,pos[1:]) if b>a); dn=len(pos)-1-up
    print(g,main,f'{sum(int(s[4]) for s in v)/1e6:.1f}Mb',len(v),'chains; adjacent pairs up/down',up,dn,'|',' '.join(f"{x[0][5:]}:{'NA' if x[1] is None else round(x[1],1)}" for x in seq))
    tot+=len(pos)-1 if len(pos)>1 else 0; good+=max(up,dn)
print(f'collinear adjacent pairs: {good}/{tot} = {good/max(tot,1):.3f}')
