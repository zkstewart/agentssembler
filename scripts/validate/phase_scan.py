#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Scan each scaffold for bins whose long-range Hi-C prefers a homolog over their own scaffold (a piece of the other haplotype, or a
phase switch). Specific pairs only, so shared sequence is not tested. Per --bin bin: own = pairs on the same scaffold at least --far
apart; other = pairs with the single most-contacted homolog (one partner, so ploidy does not pool homologs); depth = every own and
homolog pair of the bin. Thresholds calibrate themselves, as in utilities/hic_crosstalk.py's switch scan: two baselines, the median over
bins with other + own >= --minHom of other/own and of other/depth. A bin is flagged when other >= --minHom, other/own > --fold x its
baseline, and other/depth > --fold x its baseline. The second test keeps out bins whose own far contacts are merely missing (beside a
long stretch with no unique mapping, or near a scaffold end), which the old fixed rule (other >= 30 and other > 0.5 x own) flagged:
on the 292 rebuild 527 of 6262 bins. Flagged bins next to each other form runs, listed with the AGP components inside them.
Usage: phase_scan.py pairs.specific.tsv.gz final.agp [--bin 100000] [--far 2000000] [--minHom 30] [--fold 4]'''
import argparse, gzip, collections
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("pairs"); p.add_argument("agp")
p.add_argument("--bin", type=int, default=100000); p.add_argument("--far", type=int, default=2000000)
p.add_argument("--minHom", type=int, default=30); p.add_argument("--fold", type=float, default=4.0); a = p.parse_args()
B = a.bin
own = collections.Counter(); depth = collections.Counter(); hom = collections.defaultdict(collections.Counter)
hm = lambda s: s.rsplit("_", 1)[0]
for l in gzip.open(a.pairs, "rt"):
    s, x, t, y = l.rstrip("\n").split("\t")[:4]; ka = (s, int(x) // B); kb = (t, int(y) // B)
    if s == t:
        depth[ka] += 1; depth[kb] += 1
        if abs(int(x) - int(y)) >= a.far: own[ka] += 1; own[kb] += 1
    elif hm(s) == hm(t):
        depth[ka] += 1; depth[kb] += 1; hom[ka][t] += 1; hom[kb][s] += 1
comp = collections.defaultdict(list)
for l in open(a.agp):
    if l[0] != "#":
        f = l.rstrip("\n").split("\t")
        if f[4] == "W": comp[f[0]].append((int(f[1]), int(f[2]), f[5], f[8]))
vals = []
for k in set(own) | set(hom):
    other, n = max(hom[k].items(), key=lambda v: v[1]) if hom.get(k) else (".", 0)
    if n + own[k] >= a.minHom: vals.append((k, other, n, own[k], depth[k]))
med = lambda v: sorted(v)[len(v) // 2] if v else 0.0
base = med([n / max(o, 1) for _, _, n, o, _ in vals]); share = med([n / max(d, 1) for _, _, n, _, d in vals])
fl = sorted((v for v in vals if v[2] >= a.minHom and v[2] > a.fold * base * max(v[3], 1) and v[2] > a.fold * share * max(v[4], 1)), key=lambda v: v[0])
print(f"flagged bins {len(fl)} of {len(vals)} tested (bins with homolog + own far pairs >= {a.minHom}); baselines: homolog/own {base:.4f}, homolog/all {share:.4f}; fold {a.fold:g}")
run = []
for (sc, b), other, n, o, d in fl:
    if run and run[-1][0] == sc and run[-1][2] == b and run[-1][3] == other: run[-1][2] = b + 1; run[-1][4] += n; run[-1][5] += o; run[-1][6] += d
    else: run.append([sc, b, b + 1, other, n, o, d])
print("#scaffold start_Mb end_Mb homolog to_homolog to_self_far all components")
for sc, x, y, other, n, o, d in run:
    st, en = x * B, y * B
    cs = [f"{c[2]}{c[3]}:{(c[1] - c[0] + 1) // 1000}k" for c in comp[sc] if c[0] <= en and c[1] >= st]
    print(sc, f"{st / 1e6:.1f}", f"{en / 1e6:.1f}", other, n, o, d, " ".join(cs[:5]))
