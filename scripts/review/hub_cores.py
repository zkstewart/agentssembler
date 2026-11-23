#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Collapsed hubs found by read depth alone, without groups. A chain classed specific by its median depth (chain_class.py) that holds a run
of >= --minRun bp where unique-read depth (MAPQ >= 20) is >= --sharedLo: both haplotypes' reads land on the core, and the arms around it, which
read 1x with reads that are mostly MAPQ 0, pull the chain's median down. chain_forks.py finds a hub only when its other exit reaches a homolog's
chain, which depends on the groups, so with the Hi-C allele rescue on, hubs surfaced one round of breaks at a time (292). This test does not.
break_left and break_right are D04 breaks (1-based last base of the left part, for edits/chain_breaks.txt): the unitig boundary nearest each end
of the run if one lies within --snap bp of it, else the run's edge itself (break_kind inside: hifiasm collapsed the hub within one unitig);
"." when the run reaches the chain end. After the breaks the core is its own chain, classed shared (D14). Window resolution: check the edges with
the windows file before writing the breaks. High-class chains are repeat arrays, not hubs, and are left out.
Allele bubbles do not pass: both alleles are separate chains at about 1x each. Connector chains (chain8*, chain9*) are left out.
Usage: hub_cores.py depth.windows.tsv chains.agp chain_class.tsv out.tsv --sharedLo A [--minRun 50000]'''
import argparse, os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from chainlib import read_windows, core_runs
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in "win agp cls out".split(): p.add_argument(x)
p.add_argument("--sharedLo", type=float, required=True); p.add_argument("--minRun", type=int, default=50000); p.add_argument("--snap", type=int, default=20000); a = p.parse_args()
W = read_windows(a.win)
cls = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(a.cls) if l[0] != "#"}
B = {}                                                   # chain -> unitig boundaries (0-based ends of each component)
for l in open(a.agp):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if s[4] not in "NU": B.setdefault(s[0], []).append(int(s[2]))
out = []
for c, ws in W.items():
    if c.startswith(("chain8", "chain9")) or c not in cls or cls[c][6] != "specific": continue
    L = int(cls[c][1]); bd = sorted(x for x in B.get(c, []) if 0 < x < L)
    for s0, e0, u, t in core_runs(ws, a.sharedLo, a.minRun):
        def near(x):
            b = min(bd, key=lambda b: (abs(b - x), b)) if bd else None
            return (b, "unitig") if b is not None and abs(b - x) <= a.snap else (x, "inside")
        bl, kl = near(s0) if s0 > 0 else (None, "."); br, kr = near(e0) if e0 < L else (None, ".")
        arm = [w for w in ws if w[1] <= s0 or w[0] >= e0]
        tu = sum(w[3] for w in arm) / max(sum(w[2] for w in arm), 1e-9) if arm else float("nan")
        out.append((c, L, cls[c][6], cls[c][2], s0, e0, (e0 - s0) // 1000, u, t, tu, bl if bl else ".", br if br else ".", "/".join(sorted({k for k in (kl, kr) if k != "."})) or "."))
out.sort(key=lambda r: -r[6])
with open(a.out, "w") as o:
    o.write("#chain\tlength\tclass\tmedian_depth\trun_start\trun_end\trun_kb\trun_unique\trun_total\tarms_unique_frac\tbreak_left\tbreak_right\tbreak_kind\n")
    for r in out: o.write("\t".join(f"{x:.2f}" if isinstance(x, float) and i == 9 else f"{x:.0f}" if isinstance(x, float) else str(x) for i, x in enumerate(r)) + "\n")
print(f"hub cores: {len(out)} runs in {len({r[0] for r in out})} chains, {sum(r[5] - r[4] for r in out) / 1e6:.2f} Mb")
