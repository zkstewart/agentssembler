#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Per-pseudochromosome summary of how unitigs were placed and how well Hi-C supports
# their haplotype assignment (support = share of same-chromosome Hi-C weight to own
# haplotype from hic_phase.py; no-signal = < 1000 unique Hi-C reads per Mb).
import sys
from collections import defaultdict
from common import natural_key
agp, membership, support, links = sys.argv[1:5]
how = { l.split("\t")[0]: l.rstrip("\n").split("\t")[2] for l in open(membership) if not l.startswith("#") }
sup = {}
for l in open(support):
    if l.startswith("#"): continue
    sl = l.rstrip("\n").split("\t"); sup[sl[0]] = float(sl[6]) if sl[6] != "nan" else None
depth = {}
for l in open(links):
    a, b, n = l.rstrip("\n").split("\t")
    if b == ".": depth[a] = int(n)
rows = defaultdict(lambda: defaultdict(int))
for l in open(agp):
    sl = l.rstrip("\n").split("\t")
    if l.startswith("#") or sl[4] != "W": continue
    # bp of this component; Hi-C density uses the unitig's extent (components are whole
    # unitigs, or unitigs trimmed by a few overlap bp after gap filling)
    u, L, s = sl[5], int(sl[7]) - int(sl[6]) + 1, sl[0]
    span = int(sl[7])
    rows[s]["total"] += L
    h = how.get(u, "?")
    rows[s]["ref" if h == "reference+hic_phase" else "hic_placed"] += L
    if depth.get(u, 0) / span * 1e6 < 1000: rows[s]["nosignal"] += L
    elif h.startswith("gapfill") and h != "gapfill_hic": rows[s]["weak"] += L # graph-only fills: not independently phased
    elif h != "reference+hic_phase" or (sup.get(u) is not None and sup[u] >= 0.5): rows[s]["strong"] += L
    else: rows[s]["weak"] += L
print("#pseudochromosome\tlength_Mb\tref_ordered_Mb\thic_placed_Mb\tphase_strong_pct\tphase_weak_pct\tno_hic_signal_pct")
T = defaultdict(int)
for s in sorted(rows, key=natural_key):
    r = rows[s]
    for k, v in r.items(): T[k] += v
    print(f"{s}\t{r['total']/1e6:.1f}\t{r['ref']/1e6:.1f}\t{r['hic_placed']/1e6:.1f}\t{100*r['strong']/r['total']:.0f}\t{100*r['weak']/r['total']:.0f}\t{100*r['nosignal']/r['total']:.0f}")
print(f"TOTAL\t{T['total']/1e6:.1f}\t{T['ref']/1e6:.1f}\t{T['hic_placed']/1e6:.1f}\t{100*T['strong']/T['total']:.0f}\t{100*T['weak']/T['total']:.0f}\t{100*T['nosignal']/T['total']:.0f}")
