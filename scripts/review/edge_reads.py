#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Why do reads stop at gaps? Compares gap edges with control sites on the same mappings.
For each edge, and for alignments ending within --win of it while heading towards it, classifies the end:
  read_end   the last-sequenced end of the read lies here (nothing left): stall, block or molecule end
  read_start the first-sequenced end lies here (nothing left): the molecule begins here; a pore stall cannot do this
  tail       the read carries on for >= --minTail bp but stops aligning here; destination taken from the read's other alignments
Rates are per site, per kb of window, so that arrays do not dominate; medians are reported.
Usage: edge_reads.py cmp.sites.tsv reads.tsv.gz out_prefix --paf a.paf.gz [b.paf.gz] [--suffix _h1 _h2]'''
import argparse, bisect, collections, gzip, random, statistics as st
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("sites"); p.add_argument("reads"); p.add_argument("out")
p.add_argument("--paf", nargs="+", required=True); p.add_argument("--suffix", nargs="*", default=[])
p.add_argument("--minQ", type=int, default=20); p.add_argument("--slack", type=int, default=200); p.add_argument("--minTail", type=int, default=1000)
a = p.parse_args()
BINS = [(0, 1000), (1000, 2000), (2000, 5000), (5000, 10000), (10000, 20000), (20000, 50000), (50000, 100000), (100000, 200000)]
H = open(a.sites).readline().rstrip("\n").split("\t"); S = []
for l in open(a.sites):
    if l.startswith("kind"): continue
    r = dict(zip(H, l.rstrip("\n").split("\t")))
    if r["kind"] not in ("gap", "ctl"): continue
    g0, g1 = int(r["gs"]), int(r["ge"])
    S.append((r["t"], g0, "L", r["kind"])); S.append((r["t"], g1, "R", r["kind"]))
print(f"edges: gap {sum(x[3] == 'gap' for x in S)}, control {sum(x[3] == 'ctl' for x in S)}", flush=True)
stem = lambda t: t.rsplit("_h", 1)[0] if "_h" in t else t
def analyse(paf, sites):
    by = collections.defaultdict(list)
    for s in sites: by[s[0]].append(s)
    idx = {t: (sorted(v, key=lambda x: x[1])) for t, v in by.items()}
    pos = {t: [x[1] for x in v] for t, v in idx.items()}
    reads = collections.defaultdict(list)
    for l in gzip.open(paf, "rt"):
        f = l.split("\t", 12)
        reads[f[0]].append((f[5], int(f[7]), int(f[8]), f[4], int(f[2]), int(f[3]), int(f[1]), int(f[11])))
    R = {s: dict(depth=collections.Counter(), depthq=collections.Counter(), ends=collections.Counter(), dest=collections.Counter(), ids=[]) for s in sites}
    for q, al in reads.items():
        for t, ts, te, o, qs, qe, ql, mq in al:
            if t not in pos: continue
            i = bisect.bisect_left(pos[t], ts - 200000); P = pos[t]
            while i < len(P) and P[i] <= te + 200000:
                s = idx[t][i]; i += 1; _, x, side, kind = s; r = R[s]
                # distance of each bin's middle from the edge, on the flank side
                for lo, hi in BINS:
                    m = (lo + hi) // 2; c = x - m if side == "L" else x + m
                    if ts <= c < te:
                        r["depth"][lo] += 1
                        if mq >= a.minQ: r["depthq"][lo] += 1
                d = x - te if side == "L" else ts - x
                if d < 0 or d >= 200000: continue            # passes the site, or stops too far away
                if mq < a.minQ: continue
                first = (side == "L") == (o == "-")           # the end at the site is the read's first-sequenced end
                rem = qs if first else ql - qe
                b = next(lo for lo, hi in BINS if lo <= d < hi)
                if rem <= a.slack: r["ends"][(b, "read_start" if first else "read_end")] += 1
                elif rem >= a.minTail:
                    r["ends"][(b, "tail")] += 1
                    if b < 2000:
                        tq = (0, qs) if first else (qe, ql)   # the part of the read beyond the alignment
                        best, cls = 0, "none"
                        for t2, ts2, te2, o2, qs2, qe2, _, mq2 in al:
                            ov = min(qe2, tq[1]) - max(qs2, tq[0])
                            if ov < 500 or ov <= best: continue
                            if t2 == t:
                                beyond = ts2 >= x if side == "L" else te2 <= x
                                if beyond and o2 == o: k = "across_same_strand"
                                elif o2 != o and abs(((ts2 + te2) // 2) - x) < 100000: k = "foldback"
                                else: k = "same_scaffold_elsewhere"
                            elif t2.startswith("unplaced"): k = "unplaced"
                            elif stem(t2) == stem(t): k = "homolog"
                            else: k = "other_chromosome"
                            best, cls = ov, k
                        r["dest"][cls] += 1
                if b < 2000 and rem <= a.slack: r["ids"].append(q)
                elif b < 2000 and rem >= a.minTail: r["ids"].append(q)
    return R
R = {}
if a.suffix:
    for f, suf in zip(a.paf, a.suffix): R.update(analyse(f, [s for s in S if s[0].endswith(suf)]))
else: R.update(analyse(a.paf[0], S))
want = {q for r in R.values() for q in r["ids"]}
Q = {}
for l in gzip.open(a.reads, "rt"):
    i = l.split("\t", 3)
    if i[0] in want:
        try: Q[i[0]] = float(i[2])
        except ValueError: pass
med = lambda v: st.median(v) if v else float("nan")
with open(a.out + ".sites.tsv", "w") as o:
    o.write("kind\tscaffold\tpos\tside\tdepth_q20_1kb\tdepth_all_1kb\tdepth_q20_50kb\tdepth_all_50kb\tread_end_2kb\tread_start_2kb\ttail_2kb\tdest\tq_median\n")
    for s in S:
        if s not in R: continue
        r = R[s]; e = r["ends"]
        o.write(f"{s[3]}\t{s[0]}\t{s[1]}\t{s[2]}\t{r['depthq'][0]}\t{r['depth'][0]}\t{r['depthq'][20000]}\t{r['depth'][20000]}\t{e[(0, 'read_end')] + e[(1000, 'read_end')]}\t{e[(0, 'read_start')] + e[(1000, 'read_start')]}\t{e[(0, 'tail')] + e[(1000, 'tail')]}\t{dict(r['dest'])}\t{med([Q[x] for x in r['ids'] if x in Q]):.1f}\n")
with open(a.out + ".summary.txt", "w") as o:
    def P(*x):
        print(*x); print(*x, file=o)
    for kind in ("gap", "ctl"):
        K = [s for s in S if s[3] == kind and s in R]
        P(f"\n== {kind}: {len(K)} edges")
        P("distance_from_edge\tdepth_q20\tdepth_all\tread_end/kb\tread_start/kb\ttail/kb   (medians per edge; means in brackets for ends)")
        for lo, hi in BINS:
            w = (hi - lo) / 1000
            row = [f"{lo//1000}-{hi//1000}kb", med([R[s]["depthq"][lo] for s in K]), med([R[s]["depth"][lo] for s in K])]
            for c in ("read_end", "read_start", "tail"):
                v = [R[s]["ends"][(lo, c)] / w for s in K]; row.append(f"{med(v):.2f} ({sum(v) / max(len(v), 1):.2f})")
            P("\t".join(map(str, row)))
        D = collections.Counter()
        for s in K: D.update(R[s]["dest"])
        t = sum(D.values()); P("tails within 2 kb, destination:", ", ".join(f"{k} {v} ({100 * v / max(t, 1):.0f}%)" for k, v in D.most_common()))
        P("median q of reads stopping within 2 kb:", f"{med([Q[x] for s in K for x in R[s]['ids'] if x in Q]):.1f}")
