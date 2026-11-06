#!/usr/bin/env python3
'''Allelic relationships between chains, without the reference: unitig self-alignments lifted
through the chains. Each aligned block >= --minLen is lifted by its ends; blocks on one chain pair are
merged on each chain. Output: chain_a chain_b bp_a bp_b median_de strand_bp+ strand_bp- blocks.tsv for synteny.
Usage: allele_table.py chains.agp unitigs_self.paf.gz out_prefix'''
import gzip, argparse, statistics, collections, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chainlib import read_agp, lift
p=argparse.ArgumentParser(); p.add_argument("agp"); p.add_argument("paf"); p.add_argument("out"); p.add_argument("--minLen",type=int,default=20000); a=p.parse_args()
P,L=read_agp(a.agp)
blocks=[]
with gzip.open(a.paf,"rt") as fh:
    for l in fh:
        s=l.rstrip("\n").split("\t")
        q,qs,qe,sd,t,ts,te,b=s[0],int(s[2]),int(s[3]),s[4],s[5],int(s[7]),int(s[8]),int(s[10])
        if b<a.minLen or q not in P or t not in P: continue
        de=next(float(x[5:]) for x in s[12:] if x.startswith("de:f:"))
        c1,x1=lift(P,q,qs); c2,x2=lift(P,q,qe-1); d1,y1=lift(P,t,ts); d2,y2=lift(P,t,te-1)
        # a block may run into a trimmed overlap: keep it only if both ends lift to one chain each
        if not (c1 and c1==c2 and d1 and d1==d2) or c1==d1: continue
        blocks.append((c1,min(x1,x2),max(x1,x2)+1,d1,min(y1,y2),max(y1,y2)+1,de,sd))
def merged(iv):
    n=e=0
    for x,y in sorted(iv):
        if y>e: n+=y-max(x,e); e=y
    return n
pairs=collections.defaultdict(list)
for b in blocks:
    k=tuple(sorted((b[0],b[3]))); pairs[k].append(b)
with open(a.out+".pairs.tsv","w") as o:
    o.write("#chain_a\tchain_b\tbp_a\tbp_b\tmedian_de\tblocks\n")
    for (x,y),v in pairs.items():
        ia=[(b[1],b[2]) if b[0]==x else (b[4],b[5]) for b in v]; ib=[(b[4],b[5]) if b[0]==x else (b[1],b[2]) for b in v]
        o.write(f"{x}\t{y}\t{merged(ia)}\t{merged(ib)}\t{statistics.median(b[6] for b in v):.4f}\t{len(v)}\n")
with open(a.out+".blocks.tsv","w") as o:
    o.write("#chain_a\tstart_a\tend_a\tchain_b\tstart_b\tend_b\tde\tstrand\n")
    for b in blocks: o.write("\t".join(map(str,b))+"\n")
print(len(blocks),"blocks",len(pairs),"chain pairs")
