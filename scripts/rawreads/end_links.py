#!/usr/bin/env python3
"""Reads joining a sequence end (unplaced) or gap edge/scaffold end to another. Each read: alignments q>=Q, >=5 kb, adjacent on the read, each reaching within SL of a target end or gap edge, facing each other.
Usage: end_links.py reads.paf.gz out.tsv target.fasta.fai scaffolds.agp [minMapq]"""
import gzip,collections,sys
paf,out,fai,agp=sys.argv[1:5]; Q=int(sys.argv[5]) if len(sys.argv)>5 else 20
SL=2000
L={l.split('\t')[0]:int(l.split('\t')[1]) for l in open(fai)}
G=collections.defaultdict(list)
for l in open(agp):
    if l[0]=='#':continue
    s=l.split('\t')
    if s[4] in 'NU':G[s[0]].append((int(s[1])-1,int(s[2])))
def edge(t,ts,te,st,side):
    """side: 'out' end of this aln in read direction. returns label of the edge the read leaves through"""
    # read leaves target at te if strand + going forward, ts if -
    x=te if (st=='+')==(side=='next') else ts
    right=(st=='+')==(side=='next')
    if right:
        if L[t]-x<=SL: return f'{t}:END'
        for a,b in G.get(t,[]):
            if abs(x-a)<=SL: return f'{t}:{a}L'
    else:
        if x<=SL: return f'{t}:START'
        for a,b in G.get(t,[]):
            if abs(x-b)<=SL: return f'{t}:{a}R'
    return None
R=collections.defaultdict(list)
for l in gzip.open(paf,'rt'):
    s=l.split('\t',12)
    if int(s[11])<Q or int(s[8])-int(s[7])<5000: continue
    R[s[0]].append((int(s[2]),int(s[3]),s[5],int(s[7]),int(s[8]),s[4]))
C=collections.Counter()
for q,al in R.items():
    al.sort()
    for a,b in zip(al,al[1:]):
        if b[0]-a[1]>3000: continue
        ea=edge(a[2],a[3],a[4],a[5],'next'); eb=edge(b[2],b[3],b[4],b[5],'prev')
        if ea and eb and ea!=eb and not (a[2]==b[2] and ea.rsplit(':',1)[1][:-1]==eb.rsplit(':',1)[1][:-1]):
            C[tuple(sorted((ea,eb)))]+=1
with open(out,'w') as o:
    for k,v in C.most_common(): o.write(f'{k[0]}\t{k[1]}\t{v}\n')
