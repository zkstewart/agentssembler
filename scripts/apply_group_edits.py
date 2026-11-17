#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Apply the hand edits of decisions D05, D06, D07 and D14 to a grouping. Every row needs a reason.
edits.tsv columns: chain, action, value, reason. Actions:
  group <gNN>       the chain joins that group
  unassigned .      the chain leaves any group, or stops being shared
  shared <cN>       the chain is collapsed sequence of chromosome cN, copied into all its haplotypes
  shared_in <gNN,gNN>  the chain is collapsed sequence carried by these haplotype groups only, all of one chromosome. For a
                    polyploid, where a chain at 2x of 4 belongs to two haplotypes. Each entry is gNN or with:<chain>. Written to
                    out_prefix.dosage.tsv, which also carries over the dosage calls of hic_group.py for chains with no edit
  pair <cN>         the group holding this chain, after all group edits, becomes a haplotype of chromosome cN. Use it when
                    edits move a haplotype between groups, since group-to-chromosome pairing is otherwise copied unchanged.
                    Here with:<chain> means that chain's chromosome after the edits. Applied after every other row.
In place of gNN or cN, with:<chain> means that chain's group or chromosome in the input grouping (except for pair, above). Prefer it:
labels are numbered by size, so any change upstream can renumber them, and an old label would then silently point at another chromosome.
Usage: apply_group_edits.py in_prefix edits.tsv out_prefix   (reads and writes prefix.tsv, prefix.shared.tsv, prefix.chroms.tsv)
'''
import sys, os, shutil
pre, ed, out = sys.argv[1:4]
G = [l.rstrip("\n").split("\t") for l in open(pre + ".tsv")]
hdrG = [r for r in G if r[0].startswith("#")]; G = [r for r in G if not r[0].startswith("#")]
S = []; hdrS = ["#chain\tlength\tchrom\tlinks\tgroups\tscore_ratio"]
if os.path.exists(pre + ".shared.tsv"):
    S = [l.rstrip("\n").split("\t") for l in open(pre + ".shared.tsv")]
    hdrS = ["\t".join(r) for r in S if r[0].startswith("#")] or hdrS; S = [r for r in S if not r[0].startswith("#")]
L = {r[0]: r[1] for r in G + S}
Grow = {r[0]: r for r in G}; Srow = {r[0]: r for r in S}
base = {r[0]: r[2] for r in G}; chrom = {l.split("\t")[0]: l.rstrip("\n").split("\t")[1] for l in open(pre + ".chroms.tsv") if l[0] != "#"} if os.path.exists(pre + ".chroms.tsv") else {}
D = {}; hdrD = "#chain\tlength\tchrom\tmedian_depth\tdosage\tgroups\tleft_out\thic_ratio\tcall"
if os.path.exists(pre + ".dosage.tsv"):
    D = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(pre + ".dosage.tsv") if l[0] != "#" and l.strip()}
def resolve(val, act, c):
    if not val.startswith("with:"): return val
    m = val[5:]; g = base.get(m)
    if g is None or not g.startswith("g"): sys.exit(f"edit for {c}: {m} is not in a group")
    return g if act == "group" else chrom[g]
pairs = []
for l in open(ed):
    if l[0] == "#" or not l.strip(): continue
    c, act, val, *why = l.rstrip("\n").split("\t")
    if not why or not why[0].strip(): sys.exit(f"edit for {c} has no reason")
    if c not in L: sys.exit(f"unknown chain {c}")
    if act == "pair": pairs.append((c, val)); continue
    D.pop(c, None)
    if act == "shared_in":
        gs = [resolve(v, "group", c) for v in val.split(",") if v]
        if not gs or any(g not in chrom for g in gs) or len({chrom[g] for g in gs}) != 1 or len(set(gs)) != len(gs): sys.exit(f"shared_in for {c}: {val} must name distinct groups of one chromosome")
        act, val = "shared", chrom[gs[0]]
        D[c] = [c, L[c], val, ".", str(len(gs)), ",".join(gs), ",".join(g for g in chrom if chrom[g] == val and g not in gs) or ".", ".", "edit"]
    else: val = resolve(val, act, c)
    if act == "group":
        Srow.pop(c, None); Grow[c] = [c, L[c], val]
    elif act == "unassigned":
        Srow.pop(c, None); Grow[c] = [c, L[c], "unassigned"]
    elif act == "shared":
        Grow.pop(c, None); Srow[c] = Srow.get(c, [c, L[c], val, "0", ".", "0"]); Srow[c][2] = val
        if c not in {r[0] for r in G if r[2] == "shared"}: Grow[c] = [c, L[c], "shared"]
    else: sys.exit(f"unknown action {act}")
    if act == "shared": Grow[c] = [c, L[c], "shared"]
with open(out + ".tsv", "w") as o:
    for h in hdrG: o.write("\t".join(h) + "\n")
    seen = set()
    for r in G:
        o.write("\t".join(Grow.get(r[0], r)) + "\n"); seen.add(r[0])
    for c, r in Grow.items():
        if c not in seen: o.write("\t".join(r) + "\n")
with open(out + ".shared.tsv", "w") as o:
    o.write("\n".join(hdrS) + "\n")
    for c, r in Srow.items(): o.write("\t".join(r) + "\n")
if D or os.path.exists(out + ".dosage.tsv"):
    with open(out + ".dosage.tsv", "w") as o:
        o.write(hdrD + "\n")
        for r in D.values(): o.write("\t".join(r) + "\n")
if not pairs:
    if os.path.abspath(pre) != os.path.abspath(out): shutil.copyfile(pre + ".chroms.tsv", out + ".chroms.tsv")
else:
    fin = {c: r[2] for c, r in Grow.items()}
    order = [l.split("\t")[0] for l in open(pre + ".chroms.tsv") if l[0] != "#" and l.strip()]
    for c, val in pairs:
        g = fin.get(c)
        if g is None or not g.startswith("g"): sys.exit(f"pair for {c}: the chain is not in a group after the edits")
        if val.startswith("with:"):
            m = val[5:]; gm = fin.get(m)
            if gm is None or gm not in chrom: sys.exit(f"pair for {c}: {m} is not in a paired group after the edits")
            val = chrom[gm]
        print(f"pair: {g} (holding {c}) {chrom.get(g, '-')} -> {val}", file=sys.stderr)
        chrom[g] = val
        if g not in order: order.append(g)
    n = {}
    for g in order: n[chrom[g]] = n.get(chrom[g], 0) + (1 if any(x == g for x in fin.values()) else 0)
    for k, v in sorted(n.items()): print(f"chromosome {k}: {v} non-empty group(s)", file=sys.stderr)
    with open(out + ".chroms.tsv", "w") as o:
        for g in order: o.write(f"{g}\t{chrom[g]}\n")
