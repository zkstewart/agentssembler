#!/usr/bin/env python3
'''
Diploid fork. copy_number.py and round2.py drop a unitig whose windows are covered by near-identical kept copies at the
depth-implied k. Next to a collapsed node that is wrong: arms overlap the node by much of their length, so reads there
land on the node and the arm's depth is depressed. If a flank of a bubble carries both haplotypes, both must pass
through it, so the bubble needs two paths whatever the arm depths say.

Restores a removed unitig when it has a kept sibling sharing all its graph neighbours at both ends, and at least one
of those flanking neighbours is collapsed (primary depth from --lo to --hi). Unitigs listed as repeat chains stay out.
Usage: restore_arms.py copy_number.tsv kept2.txt repeat_chains.txt unitigs.gfa out.tsv --lo L --hi H (about 1.54 and 2.68 x haplotype depth)
'''
import argparse, collections
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in "cn kept rep gfa out".split(): p.add_argument(x)
p.add_argument("--lo", type=float, required=True); p.add_argument("--hi", type=float, required=True)
a = p.parse_args()
inv = {"+": "-", "-": "+"}
R = {l.split("\t")[0]: l.rstrip().split("\t") for l in open(a.cn) if l[0] != "#"}
L = {u: int(r[1]) for u, r in R.items()}; D = {u: float(r[2]) for u, r in R.items()}
kept = {l.strip() for l in open(a.kept) if l.strip()}
rep = set()
for l in open(a.rep):
    rep |= set(l.rstrip().split("\t")[2].split(","))
nb = collections.defaultdict(dict)
for l in open(a.gfa):
    if l[0] != "L": continue
    s = l.split("\t"); ov = int(s[5].strip()[:-1])
    nb[(s[1], s[2])][(s[3], s[4])] = ov
    nb[(s[3], inv[s[4]])][(s[1], inv[s[2]])] = ov
# a sibling enters and leaves through the same neighbours, in the same orientation
def sig(u, o):
    right = frozenset(nb[(u, o)]); left = frozenset(nb[(u, inv[o])])
    return left, right
index = collections.defaultdict(list)
for u in kept:
    for o in "+-": index[sig(u, o)].append((u, o))
flank_is_shared = lambda s: any(a.lo <= D[x] < a.hi for side in s for x, _ in side)
rows = []
for u in R:
    if u in kept or u in rep: continue
    for o in "+":
        s = sig(u, o)
        if not s[0] or not s[1]: continue               # a tip cannot be a bubble arm
        sib = [v for v, _ in index.get(s, []) + index.get((frozenset((x, inv[y]) for x, y in s[1]), frozenset((x, inv[y]) for x, y in s[0])), []) if v != u]
        if sib and flank_is_shared(s):
            ovs = sorted(nb[(u, o)].values())[-1] + sorted(nb[(u, inv[o])].values())[-1]
            rows.append((u, L[u], D[u], sib[0], L[u] - ovs))
with open(a.out, "w") as o:
    o.write("#unitig\tlength\tprimary\tsibling\tlength_minus_overlaps\n")
    for r in sorted(rows, key=lambda r: -r[1]): o.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]}\n")
print(f"restore {len(rows)} arms, {sum(r[1] for r in rows)/1e6:.2f} Mb, {sum(max(r[4],0) for r in rows)/1e6:.2f} Mb outside overlaps")
