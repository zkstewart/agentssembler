#!/usr/bin/env python3
'''Runs of about 2x depth inside a chain, smoothed over 5 windows. A chain's overall depth hides a collapsed stretch embedded in it, which is where
hub nominations for D14 come from when the run sits away from the chain ends. Candidates only: test each with hub_depth.py.
Usage: embedded_runs.py depth.windows.tsv chains.agp unitigs.gfa --hap H [--groups groups.tsv] [--minRun 100000]'''
import argparse, collections, statistics as st
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in "windows agp gfa".split(): p.add_argument(x)
p.add_argument("--hap", type=float, required=True); p.add_argument("--groups"); p.add_argument("--minRun", type=int, default=100000); a = p.parse_args()
lo, hi = 1.54 * a.hap, 2.68 * a.hap
W = collections.defaultdict(list)
for l in open(a.windows):
    if l[0] != "#":
        s = l.split("\t"); W[s[0]].append((int(s[1]), int(s[2]), float(s[3])))
AG = collections.defaultdict(list)
for l in open(a.agp):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t")
        if s[4] == "W": AG[s[0]].append((int(s[1]), int(s[2]), s[5]))
inv = {"+": "-", "-": "+"}
adj = collections.defaultdict(list)
for l in open(a.gfa):
    if l[0] == "L":
        s = l.split("\t"); adj[s[1], s[2]].append(s[3]); adj[s[3], inv[s[4]]].append(s[1])
G = {}
if a.groups:
    G = {l.split()[0]: l.split()[2] for l in open(a.groups) if l[0] != "#"}
print("#chain\tgroup\tchain_len\trun_start\trun_end\trun_bp\tchain_median\tunitig_at_start\tdegree_at_start\tunitig_at_end\tdegree_at_end")
n = tot = 0
for c, v in W.items():
    if len(v) < 15 or c not in AG: continue
    sm = [st.median(x[2] for x in v[max(0, i - 2):i + 3]) for i in range(len(v))]
    i = 0
    while i < len(v):
        if not lo <= sm[i] <= hi:
            i += 1; continue
        j = i
        while j + 1 < len(v) and lo <= sm[j + 1] <= hi: j += 1
        bp = v[j][1] - v[i][0]
        if bp >= a.minRun and bp < 0.9 * v[-1][1]:
            ends = []
            for pos in (v[i][0], v[j][1]):
                u = min(AG[c], key=lambda r: min(abs(r[0] - pos), abs(r[1] - pos)))[2]
                ends += [u, f"{len(adj[u, '+'])}/{len(adj[u, '-'])}"]
            n += 1; tot += bp
            print(c, G.get(c, "."), v[-1][1], v[i][0], v[j][1], bp, f"{st.median(sm):.0f}", *ends, sep="\t")
        i = j + 1
print(f"# embedded 2x runs: {n}, {tot / 1e6:.1f} Mb")
