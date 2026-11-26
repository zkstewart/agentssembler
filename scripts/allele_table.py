#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Allelic relationships between chains, without the reference: unitig self-alignments lifted
through the chains. Each aligned block >= --minLen is lifted by its ends; blocks on one chain pair are
merged on each chain. Output: chain_a chain_b bp_a bp_b median_de blocks in pairs.tsv, and each block in blocks.tsv with its strand
between the chains: the alignment's unitig strand, reversed once for each unitig that lies reversed in its chain (project_rescued.py
carries positions through it).
A dovetail is not an allele: a near-identical block (de <= --overlapDe) at an end of each chain, placed so that one chain runs on into
the other (end to start on the same strand, or start to start or end to end on opposite strands), and under --maxOverlapFrac of either
chain. It is the graph overlap between neighbours on a path, so it says nothing about different haplotypes: written to dovetails.tsv
(same columns as blocks.tsv) and left out of pairs.tsv and blocks.tsv. Read as alleles, such overlaps made the Hi-C rescue put the
continuation of a chain in the other haplotype, against unique Hi-C (RECORD.md). Two alleles that share
a homozygous start run parallel (start to start, same strand), and a near-full copy covers most of a chain: both stay alleles.
Usage: allele_table.py chains.agp unitigs_self.paf.gz out_prefix [--overlapDe 0.002] [--endTol 5000] [--maxOverlapFrac 0.8]'''
import math, gzip, argparse, statistics, collections, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chainlib import read_agp, lift_strand
p=argparse.ArgumentParser(); p.add_argument("agp"); p.add_argument("paf"); p.add_argument("out"); p.add_argument("--minLen",type=int,default=20000)
p.add_argument("--overlapDe",type=float,default=0.002); p.add_argument("--endTol",type=int,default=5000); p.add_argument("--maxOverlapFrac",type=float,default=0.8); a=p.parse_args()
P,L=read_agp(a.agp)
blocks=[]
with gzip.open(a.paf,"rt") as fh:
    for l in fh:
        s=l.rstrip("\n").split("\t")
        q,qs,qe,sd,t,ts,te,b=s[0],int(s[2]),int(s[3]),s[4],s[5],int(s[7]),int(s[8]),int(s[10])
        if b<a.minLen or q not in P or t not in P: continue
        de=next((float(x[5:]) for x in s[12:] if x.startswith("de:f:")),float("nan"))
        c1,x1,o1=lift_strand(P,q,qs); c2,x2,_=lift_strand(P,q,qe-1); d1,y1,o2=lift_strand(P,t,ts); d2,y2,_=lift_strand(P,t,te-1)
        # a block may run into a trimmed overlap: keep it only if both ends lift to one chain each
        if not (c1 and c1==c2 and d1 and d1==d2) or c1==d1: continue
        cs=sd if o1==o2 else ("-" if sd=="+" else "+")        # chain-level strand
        blocks.append((c1,min(x1,x2),max(x1,x2)+1,d1,min(y1,y2),max(y1,y2)+1,de,cs))
def dovetail(b):
    c,x,y,d,u,v,de,st=b
    if not de<=a.overlapDe or y-x>=a.maxOverlapFrac*L[c] or v-u>=a.maxOverlapFrac*L[d]: return False
    cs,ce,ds,dE=x<=a.endTol,y>=L[c]-a.endTol,u<=a.endTol,v>=L[d]-a.endTol
    return (ce and ds) or (cs and dE) if st=="+" else (cs and ds) or (ce and dE)
dove=[b for b in blocks if dovetail(b)]; blocks=[b for b in blocks if not dovetail(b)]
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
        o.write(f"{x}\t{y}\t{merged(ia)}\t{merged(ib)}\t{(statistics.median(de) if (de:=[b[6] for b in v if not math.isnan(b[6])]) else float('nan')):.4f}\t{len(v)}\n")
with open(a.out+".blocks.tsv","w") as o:
    o.write("#chain_a\tstart_a\tend_a\tchain_b\tstart_b\tend_b\tde\tstrand\n")
    for b in blocks: o.write("\t".join(map(str,b))+"\n")
with open(a.out+".dovetails.tsv","w") as o:
    o.write("#chain_a\tstart_a\tend_a\tchain_b\tstart_b\tend_b\tde\tstrand\n")
    for b in dove: o.write("\t".join(map(str,b))+"\n")
print(len(blocks),"blocks",len(pairs),"chain pairs;",len(dove),"dovetails (graph overlaps, not alleles) set aside")
