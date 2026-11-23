#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Collapsed hubs found by read depth alone, without groups. A chain classed specific by its median depth (chain_class.py) that holds a run
of >= --minRun bp where unique-read depth (MAPQ >= 20) is >= --sharedLo: both haplotypes' reads land on the core, and the arms around it, which
read 1x with reads that are mostly MAPQ 0, pull the chain's median down. chain_forks.py finds a hub only when its other exit reaches a homolog's
chain, which depends on the groups, so with the Hi-C allele rescue on, hubs surfaced one round of breaks at a time (292). This test does not.
break_left and break_right are D04 breaks (1-based last base of the left part, for edits/chain_breaks.txt): the unitig boundary nearest each end
of the core's unitigs, those lying at least half inside the run (292: the 2x run's edges in 10 kb windows fell 25-90 kb inside them).
break_chain.py splits chains only between unitigs, so a run with no such unitig gets "." and kind inside: hifiasm collapsed the hub within
one unitig. An edge at the chain's end needs no break (kind end). breakable is yes when every
inner edge has a break, and the piece the breaks free reads 2x unique in >= --minCore2x of its windows (core_2x_frac): after them the core
is its own chain, classed shared (D14). A freed piece with long 1x arms would put them on both haplotypes (292's chain0119 break did). A core that is not breakable stays in its chain; the choice
is then the chain on one haplotype (its core missing from the other) or the whole chain shared (its arms duplicated), at D14.
High-class chains are repeat arrays, not hubs, and are left out.
Allele bubbles do not pass: both alleles are separate chains at about 1x each. Connector chains (chain8*, chain9*) are left out.
Stage 05 runs it before D04, so that hub breaks go in with the chimera breaks: on stage 04's chains (07_chains2, with class "-": the class
is then the windows' median depth against --sharedLo), with --liftAgp lifting chains and positions onto the chains D04 breaks
(08_final_chains, which D03's breaks may have split). A break that lifts onto a chain end is no longer needed (both sides of the boundary are lifted, so a
chain stored reversed lifts exactly); a core whose breaks are all gone is already its own chain and is not breakable. Stage 06 runs it again on the
broken chains, as a check.
Usage: hub_cores.py depth.windows.tsv chains.agp chain_class.tsv|- out.tsv --sharedLo A [--minRun 50000] [--liftAgp chains.agp]'''
import argparse, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from chainlib import read_windows, core_runs, read_agp, lift
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in "win agp cls out".split(): p.add_argument(x)
p.add_argument("--sharedLo", type=float, required=True); p.add_argument("--minRun", type=int, default=50000); p.add_argument("--liftAgp"); p.add_argument("--minCore2x", type=float, default=0.75); a = p.parse_args()
W = read_windows(a.win)
if a.cls == "-":   # class from the windows alone: median depth against --sharedLo
    cls = {}
    for c, ws in W.items():
        m = sorted(w[2] for w in ws)[len(ws) // 2]; cls[c] = [c, str(ws[-1][1]), f"{m:.1f}", "", "", "", "specific" if m < a.sharedLo else "shared"]
else: cls = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(a.cls) if l[0] != "#"}
ROWS = {}                                                 # chain -> AGP components (chain start, end, unitig, unitig start, end, strand)
for l in open(a.agp):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t")
        if s[4] not in "NU": ROWS.setdefault(s[0], []).append((int(s[1]) - 1, int(s[2]), s[5], int(s[6]) - 1, int(s[7]), s[8]))
if a.liftAgp: P2, L2 = read_agp(a.liftAgp)
def to2(c, x):
    """0-based position x on chain c of the input AGP -> (chain, position) on --liftAgp's chains, or (None, None)"""
    for c0, c1, u, u0, u1, st in ROWS.get(c, ()):
        if c0 <= x < c1: return lift(P2, u, u0 + (x - c0) if st == "+" else u1 - 1 - (x - c0))
    return None, None
out = []
for c, ws in W.items():
    if c.startswith(("chain8", "chain9")) or c not in cls or cls[c][6] != "specific": continue
    L = int(cls[c][1])
    for s0, e0, u, t in core_runs(ws, a.sharedLo, a.minRun):
        # the core's unitigs: those lying at least half inside the run. The breaks go at their outer boundaries.
        core = [(x0, x1) for x0, x1, *_ in ROWS.get(c, ()) if min(x1, e0) - max(x0, s0) >= 0.5 * (x1 - x0)]
        if core: (bl, kl), (br, kr) = ((core[0][0], "unitig") if core[0][0] > 0 else (None, "end")), ((core[-1][1], "unitig") if core[-1][1] < L else (None, "end"))
        else: bl, kl, br, kr = None, "inside" if s0 > 0 else "end", None, "inside" if e0 < L else "end"
        f0, f1 = (bl or 0, br or L) if "inside" not in (kl, kr) else (s0, e0)          # the piece the breaks would free
        fw = [w for w in ws if w[0] >= f0 and w[1] <= f1]
        f2 = sum(w[3] >= a.sharedLo for w in fw) / len(fw) if fw else 0.0
        arm = [w for w in ws if w[1] <= s0 or w[0] >= e0]
        tu = sum(w[3] for w in arm) / max(sum(w[2] for w in arm), 1e-9) if arm else float("nan")
        out.append((c, L, cls[c][6], cls[c][2], s0, e0, (e0 - s0) // 1000, u, t, tu, bl if bl else ".", br if br else ".", f"{kl}/{kr}", round(f2, 2), "no" if "inside" in (kl, kr) or f2 < a.minCore2x else "yes"))
if a.liftAgp:
    lifted = []
    for r in out:
        r = list(r); c = r[0]
        cs, xs = to2(c, r[4]); ce, xe = to2(c, r[5] - 1)
        if cs is None or cs != ce: continue                  # the core was split or dropped by D03's breaks
        r[0], r[1] = cs, L2[cs]; r[4], r[5] = min(xs, xe), max(xs, xe) + 1
        bk = r[12].split("/")
        for i, k in ((10, 0), (11, 1)):
            if r[i] == ".": continue
            b = int(r[i])                                    # boundary between bases b-1 and b (0-based): lift both sides
            cl, xl = to2(c, b - 1); cr, xr = to2(c, b)
            if cl != cs or cr != cs: r[i] = "."; bk[k] = "end"; continue  # already a chain end after D03
            r[i] = min(xl, xr) + 1                           # 1-based last base of the left part, either strand
        if xs > xe:                                          # the chain lies reversed in the lifted AGP: swap sides
            r[10], r[11] = r[11], r[10]; bk.reverse()
        r[12] = "/".join(bk)
        if r[10] == "." and r[11] == ".": r[14] = "no"      # D03 already freed the core (or nothing is left to break)
        lifted.append(tuple(r))
    out = lifted
out.sort(key=lambda r: -r[6])
with open(a.out, "w") as o:
    o.write("#chain\tlength\tclass\tmedian_depth\trun_start\trun_end\trun_kb\trun_unique\trun_total\tarms_unique_frac\tbreak_left\tbreak_right\tbreak_kind\tcore_2x_frac\tbreakable\n")
    for r in out: o.write("\t".join(f"{x:.2f}" if isinstance(x, float) and i in (9, 13) else f"{x:.0f}" if isinstance(x, float) else str(x) for i, x in enumerate(r)) + "\n")
print(f"hub cores: {len(out)} runs in {len({r[0] for r in out})} chains, {sum(r[5] - r[4] for r in out) / 1e6:.2f} Mb; breakable {sum(r[-1] == 'yes' for r in out)}")
