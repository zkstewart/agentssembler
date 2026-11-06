#!/usr/bin/env python3
'''
Apply the hand edits of decisions D05, D06, D07 and D14 to a grouping. Every row needs a reason.
edits.tsv columns: chain, action, value, reason. Actions:
  group <gNN>       the chain joins that group
  unassigned .      the chain leaves any group, or stops being shared
  shared <cN>       the chain is collapsed sequence of chromosome cN, copied into all its haplotypes
In place of gNN or cN, with:<chain> means that chain's group or chromosome in the input grouping. Prefer it: labels are numbered by
size, so any change upstream can renumber them, and an old label would then silently point at another chromosome.
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
def resolve(val, act, c):
    if not val.startswith("with:"): return val
    m = val[5:]; g = base.get(m)
    if g is None or not g.startswith("g"): sys.exit(f"edit for {c}: {m} is not in a group")
    return g if act == "group" else chrom[g]
for l in open(ed):
    if l[0] == "#" or not l.strip(): continue
    c, act, val, *why = l.rstrip("\n").split("\t")
    if not why or not why[0].strip(): sys.exit(f"edit for {c} has no reason")
    if c not in L: sys.exit(f"unknown chain {c}")
    val = resolve(val, act, c)
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
if os.path.abspath(pre) != os.path.abspath(out): shutil.copyfile(pre + ".chroms.tsv", out + ".chroms.tsv")
