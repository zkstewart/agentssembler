#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Does grouping need the Hi-C allele rescue (settings.hic_allele_rescue)? Read at D04, before any grouping is built on this Hi-C.
Where two haplotypes are near-identical, most Hi-C reads on them have MAPQ 0 and stage 05 keeps none of their pairs, so Hi-C
cannot group those chains and they end unassigned. Per chain of >= --minLen with an allele partner: unique pairs per Mb (both
ends MAPQ >= 1, stage 05's pairs.tsv.gz) and allele-ambiguous ends per Mb (hic_rescue_pairs.py). A chain is Hi-C-blind when its
unique density is under --lowFrac x the length-weighted median of all such chains and its ambiguous ends outnumber its unique
ones >= --ratio times. On 292 (near-identical chr2 haplotypes) 111 chains, 25.9 Mb were rescued: unique 32 pairs/Mb (median)
against 1,577 for the chains Hi-C grouped; ambiguous / unique about 2,000 against 3.
With --genes (compleasm on the unitigs) and --agp, counts gene copies on Hi-C-blind chains.
Usage: rescue_signal.py pairs.tsv.gz pairs.rescued.tsv.gz chains.fasta.fai alleles.pairs.tsv out.tsv [--genes T --agp chains.agp]'''
import argparse, collections, gzip, os, statistics, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from chainlib import read_agp, gene_hits
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in ("unique", "rescued", "fai", "alleles", "out"): p.add_argument(x)
p.add_argument("--genes"); p.add_argument("--agp")
p.add_argument("--minLen", type=int, default=50000); p.add_argument("--lowFrac", type=float, default=0.1); p.add_argument("--ratio", type=float, default=10.0)
p.add_argument("--alarmFrac", type=float, default=0.01, help="alarm when Hi-C-blind chains hold this fraction of the chain sequence")
a = p.parse_args()
L = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.fai)}
al = set()
for l in open(a.alleles):
    if l[0] != "#": s = l.split("\t"); al |= {s[0], s[1]}
U = collections.Counter(); A = collections.Counter()
for l in gzip.open(a.unique, "rt"):
    s = l.split("\t"); U[s[0]] += 1; U[s[2]] += 1
kinds = collections.Counter()
for l in gzip.open(a.rescued, "rt"):
    s = l.rstrip("\n").split("\t"); k = s[4]; kinds[k] += 1
    if k[0] == "a": A[s[0]] += 1
    if k[1] == "a": A[s[2]] += 1
cs = [c for c in L if L[c] >= a.minLen and c in al]
ud = {c: U[c] / L[c] * 1e6 for c in cs}
w = sorted(cs, key=lambda c: ud[c]); half = sum(L[c] for c in cs) / 2; acc = 0; med = 0
for c in w:
    acc += L[c]
    if acc >= half: med = ud[c]; break
blind = [c for c in cs if ud[c] < a.lowFrac * med and A[c] >= a.ratio * max(U[c], 1)]
G = collections.Counter()
if a.genes and a.agp and os.path.exists(a.genes):
    P, _ = read_agp(a.agp)
    for g, st, c, x in gene_hits(a.genes, P):
        if c: G[c] += 1
with open(a.out, "w") as o:
    o.write("#chain\tlength\tunique_per_Mb\tambiguous_ends_per_Mb\tgene_copies\tblind\n")
    for c in sorted(cs, key=lambda c: (c not in blind, -L[c])):
        o.write(f"{c}\t{L[c]}\t{ud[c]:.0f}\t{A[c] / L[c] * 1e6:.0f}\t{G[c] if a.genes else '.'}\t{'yes' if c in blind else 'no'}\n")
tot = sum(L.values()); bb = sum(L[c] for c in blind); uu = U.total() // 2 if hasattr(U, "total") else sum(U.values()) // 2
amb = sum(v for k, v in kinds.items() if "a" in k)
print(f"Hi-C pairs kept by stage 05 (both ends unique): {uu:,}; with an allele-ambiguous end (not used unless hic_allele_rescue): {amb:,} ({amb / max(uu, 1):.2f} per unique pair).")
print(f"Chains of {a.minLen // 1000} kb or more with an allele partner: {len(cs)}; length-weighted median unique density {med:.0f} pairs/Mb.")
gtxt = f", carrying {sum(G[c] for c in blind)} gene copies" if a.genes else ""
print(f"Hi-C-blind (unique density under {a.lowFrac:g} x median, ambiguous ends >= {a.ratio:g} x unique): {len(blind)} chains, {bb / 1e6:.1f} Mb "
      f"({100 * bb / max(tot, 1):.1f}% of chain sequence){gtxt}. Table: {a.out}.")
if bb >= a.alarmFrac * tot:
    print(f"ALARM: {bb / 1e6:.1f} Mb of allele-paired chains have almost no unique Hi-C: grouping will leave them unassigned. "
          "Decide on settings.hic_allele_rescue now (D04): it changes stages 06-08, and turning it on later costs a reset of 07_group.")
