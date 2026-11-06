#!/usr/bin/env python3
'''Judge a chain grouping against independent evidence (never used to build it): reference chromosome
of each group from unitig alignments, and baseline haplotype composition.
Usage: eval_groups.py groups.tsv chains.agp unitigs_vs_ref.paf.gz baseline.agp'''
import sys, gzip, collections
g,agp,paf,base=sys.argv[1:5]
gr={};Ln={}
for l in open(g):
    if l[0]=='#': continue
    c,n,x=l.rstrip().split('\t'); gr[c]=x; Ln[c]=int(n)
u2c=collections.defaultdict(list)
for l in open(agp):
    if l[0]=='#': continue
    s=l.rstrip().split('\t'); u2c[s[5]].append(s[0])
ref=collections.defaultdict(collections.Counter)
with gzip.open(paf,'rt') as fh:
    for l in fh:
        s=l.split('\t',12)
        if 'tp:A:P' in s[12] and int(s[11])>=20 and s[0] in u2c:
            for c in u2c[s[0]]: ref[gr[c]][s[5]]+=int(s[3])-int(s[2])
bs=collections.defaultdict(collections.Counter)
for l in open(base):
    if l[0]=='#': continue
    s=l.rstrip().split('\t')
    if s[4]=='W' and s[5] in u2c:
        for c in u2c[s[5]]: bs[gr[c]][s[0]]+=int(s[7])-int(s[6])+1
sz=collections.Counter()
for c,x in gr.items(): sz[x]+=Ln[c]
print('group\tMb\tchains\tref_top\tshare\tbaseline_top\tshare')
agree=tot=0
for x in sorted(sz):
    n=sum(1 for c in gr if gr[c]==x); r=ref[x].most_common(2); b=bs[x].most_common(2)
    rt=sum(ref[x].values()) or 1; bt=sum(bs[x].values()) or 1
    print(x,f'{sz[x]/1e6:.1f}',n,r[0][0] if r else '.',f'{r[0][1]/rt:.2f}' if r else 0,b[0][0] if b else '.',f'{b[0][1]/bt:.2f}' if b else 0,sep='\t')
    if x!='unassigned' and b: agree+=b[0][1]; tot+=bt
print(f'bp in the dominant baseline haplotype of each group: {agree/tot:.3f}')
