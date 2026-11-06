#!/usr/bin/env python3
'''Reads through every sequence-merging join: adjacent AGP components with no gap, as a bridge or graph link leaves them. Read fill
counts reads over gaps, so these joins are otherwise never checked at final coordinates. For each join, primary alignments covering
it by --anchors bp on both sides, at any MAPQ and at MAPQ >= 20, with the median depth --flank bp either side for comparison.
Under 3 reads at the longest anchor flags a join. Joins inside a duplicated stretch need masked PAFs: a read cannot anchor otherwise.
Usage: join_reads.py final.agp out.tsv reads.paf.gz [more.paf.gz ...] [--anchors 2000,10000 --flank 30000]'''
import argparse, bisect, collections, gzip, statistics as st
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("agp"); p.add_argument("out"); p.add_argument("paf", nargs="+"); p.add_argument("--anchors", default="2000,10000"); p.add_argument("--flank", type=int, default=30000)
a = p.parse_args(); AN = [int(x) for x in a.anchors.split(",")]
rows = collections.defaultdict(list)
for l in open(a.agp):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); rows[s[0]].append(s)
J = collections.defaultdict(list)
for sc, r in rows.items():
    for x, y in zip(r, r[1:]):
        if x[4] == "W" and y[4] == "W": J[sc].append({"pos": int(x[2]), "left": x[5] + x[8], "right": y[5] + y[8], "n": [0] * len(AN), "q": [0] * len(AN), "dl": 0, "dr": 0})
start = {sc: [j["pos"] for j in v] for sc, v in J.items()}
mx, F = max(AN), a.flank
seen = set()
for f in a.paf:
    with gzip.open(f, "rt") as fh:
        for l in fh:
            s = l.split("\t", 12); t = s[5]
            if t not in J or "tp:A:P" not in s[12]: continue
            x, y, q = int(s[7]), int(s[8]), int(s[11]); v = J[t]
            i = bisect.bisect_left(start[t], x - F - 1)
            while i < len(v) and v[i]["pos"] <= y + F:
                j = v[i]; p0 = j["pos"]; i += 1
                for k, an in enumerate(AN):
                    if x <= p0 - an and y >= p0 + an:
                        j["n"][k] += 1
                        if q >= 20: j["q"][k] += 1
                if x <= p0 - F and y >= p0 - F: j["dl"] += 1
                if x <= p0 + F and y >= p0 + F: j["dr"] += 1
H = "\t".join(f"reads_{x}\tq20_{x}" for x in AN)
n = bad = 0
with open(a.out, "w") as o:
    o.write(f"#scaffold\tpos\tleft\tright\t{H}\tdepth_left\tdepth_right\tflag\n")
    for sc, v in J.items():
        for j in v:
            n += 1; fl = "weak" if j["n"][-1] < 3 else "no_q20" if j["q"][-1] == 0 else "ok"; bad += fl != "ok"
            c = "\t".join(f"{j['n'][k]}\t{j['q'][k]}" for k in range(len(AN)))
            o.write(f"{sc}\t{j['pos']}\t{j['left']}\t{j['right']}\t{c}\t{j['dl']}\t{j['dr']}\t{fl}\n")
print(f"{n} merged joins, {bad} flagged")
