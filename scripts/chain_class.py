#!/usr/bin/env python3
'''
Depth class per chain, from a unitig_depth.py windows file. shared: median depth in [--sharedLo, --sharedHi). high: at or above --sharedHi. specific: below.
Also writes the shared chains to --sharedOut. A candidate list only: see docs/decisions/D14_shared_chains.md.
Usage: chain_class.py depth.windows.tsv chains.fai chain_repeat.tsv out.tsv --oneLo L --sharedLo A --sharedHi B [--sharedOut shared.txt]
'''
import argparse, collections
import numpy as np
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in ("win", "fai", "rep", "out"): p.add_argument(x)
p.add_argument("--oneLo", type=float, required=True); p.add_argument("--sharedLo", type=float, required=True); p.add_argument("--sharedHi", type=float, required=True); p.add_argument("--sharedOut")
a = p.parse_args()
W = collections.defaultdict(list)
for l in open(a.win):
    if l[0] != "#":
        s = l.split("\t"); W[s[0]].append(float(s[3]))
L = {l.split()[0]: int(l.split()[1]) for l in open(a.fai)}
R = {l.split("\t")[0]: float(l.split("\t")[2]) for l in open(a.rep) if l[0] != "#"}
sh = []
with open(a.out, "w") as o:
    o.write("#chain\tlength\tmedian_depth\tfrac_2x\tfrac_1x\trepeat\tclass\n")
    for c, v in W.items():
        v = np.array(v); m = float(np.median(v)); f2 = float(((v >= a.sharedLo) & (v < a.sharedHi)).mean()); f1 = float(((v >= a.oneLo) & (v < a.sharedLo)).mean())
        k = "shared" if a.sharedLo <= m < a.sharedHi else ("high" if m >= a.sharedHi else "specific")
        o.write(f"{c}\t{L[c]}\t{m:.1f}\t{f2:.2f}\t{f1:.2f}\t{R.get(c, 0):.2f}\t{k}\n")
        if k == "shared": sh.append(c)
if a.sharedOut: open(a.sharedOut, "w").write("\n".join(sh) + "\n")
print(f"shared {len(sh)} chains, {sum(L[c] for c in sh) / 1e6:.1f} Mb")
