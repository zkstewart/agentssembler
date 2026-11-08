#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Are reads really absent near gaps, or placed elsewhere? Maps a sample of raw reads to short windows of flank sequence alone, so nothing competes: each window's depth counts reads carrying that sequence, from any copy.
A read counts only if it covers at least 90% of the window at 97% identity or better: with looser rules any transposon fragment counts, and the test cannot report absence. Control sites should sit near 1-2 x haplotype depth; check that first.
Satellite windows make base-level alignment crawl: run once without --loose and without -c to get a loose table, then pass it with --maxLoose (about 10 x haplotype depth) to skip them.
Writes PREFIX.strict.tsv, and PREFIX.hits.tsv (window, read) for tracing where those reads are placed.
Usage: edge_seq_depth.py cmp.sites.tsv genome.fasta PREFIX --bam a.bam [b.bam] --frac 0.2 --minimap2 M --samtools S [--controls 150] [--loose table --maxLoose D]'''
import argparse, collections, random, statistics as st, subprocess
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("sites"); p.add_argument("fasta"); p.add_argument("out"); p.add_argument("--bam", nargs="+", required=True); p.add_argument("--frac", type=float, default=0.1)
p.add_argument("--minimap2", required=True); p.add_argument("--samtools", required=True); p.add_argument("--controls", type=int, default=150); p.add_argument("--threads", default="4"); p.add_argument("--loose"); p.add_argument("--maxLoose", type=float, default=1e12); a = p.parse_args()
LOOSE = {}
if a.loose:
    for l in open(a.loose):
        f = l.rstrip("\n").split("\t")
        if f[0] != "kind": LOOSE["|".join(f[:5])] = float(f[5])
DIST = [0, 5000, 20000, 50000, 100000, 200000]; W = 5000
H = open(a.sites).readline().rstrip("\n").split("\t"); G, C = [], []
for l in open(a.sites):
    if l.startswith("kind"): continue
    r = dict(zip(H, l.rstrip("\n").split("\t")))
    if r["kind"] == "gap": G += [(r["t"], int(r["gs"]), "L"), (r["t"], int(r["ge"]), "R")]
    elif r["kind"] == "ctl": C += [(r["t"], int(r["gs"]), "L")]
random.seed(1); C = random.sample(C, min(a.controls, len(C)))
L = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.fasta + ".fai")}
def get(t, s, e):
    x = subprocess.run([a.samtools, "faidx", a.fasta, f"{t}:{s + 1}-{e}"], capture_output=True, text=True, check=True).stdout.split("\n", 1)[1].replace("\n", "")
    return x
names = {}; bad = 0
with open(a.out + ".windows.fa", "w") as o:
    for kind, sites in (("gap", G), ("ctl", C)):
        for t, x, side in sites:
            if kind == "gap":
                probe = get(t, max(x - 1, 0), x) if side == "L" else get(t, x, x + 1)
                inside = get(t, x, x + 1) if side == "L" else get(t, x - 1, x)
                if inside.upper() != "N": bad += 1
            for d in DIST:
                s, e = (x - d - W, x - d) if side == "L" else (x + d, x + d + W)
                if s < 0 or e > L[t]: continue
                q = get(t, s, e)
                if q.upper().count("N") > 100: continue
                n = f"{kind}|{t}|{x}|{side}|{d}"
                if LOOSE and LOOSE.get(n, 0) > a.maxLoose: continue
                names[n] = 1; o.write(f">{n}\n{q}\n")
print(f"windows {len(names)}; gap sites whose inside was not N (coordinates off if many): {bad} of {len(G)}", flush=True)
stream = " ; ".join(f"{a.samtools} view -@ 2 -s {a.frac} {b}" for b in a.bam)
cmd = f"( {stream} ) | awk -F'\\t' 'length($10)>=5000{{print \">\"$1\"\\n\"$10}}' | {a.minimap2} -x lr:hq -t {a.threads} -c --secondary=yes -N 20 -p 0.5 {a.out}.windows.fa - 2>/dev/null"
P = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, text=True)
best = collections.defaultdict(dict)
for l in P.stdout:
    f = l.split("\t", 12); w, r = f[5], f[0]; b = int(f[8]) - int(f[7])
    if b >= 0.9 * W and int(f[9]) / max(int(f[10]), 1) >= 0.97: best[w][r] = W
P.wait()
D = {n: sum(best[n].values()) / W / a.frac for n in names}
with open(a.out + ".hits.tsv", "w") as o:
    for n, v in best.items():
        for r in v: o.write(f"{n}\t{r}\n")
with open(a.out + ".strict.tsv", "w") as o:
    o.write("kind\tscaffold\tpos\tside\tdist\tseq_depth\n")
    for n, v in D.items(): o.write("\t".join(n.split("|")) + f"\t{v:.1f}\n")
ctl = st.median(v for n, v in D.items() if n.startswith("ctl"))
print(f"control median sequence depth {ctl:.1f}")
print("dist_kb\tn\tp10\tq1\tmedian\tq3\tp90\t(ratio to control median)\tfrac_under_0.5x\tfrac_over_2x")
for d in DIST:
    v = sorted(x for n, x in D.items() if n.startswith("gap") and n.endswith(f"|{d}"))
    if v:
        q = lambda f: v[int(f * (len(v) - 1))]
        print(f"{d // 1000}-{d // 1000 + 5}\t{len(v)}\t{q(.1):.0f}\t{q(.25):.0f}\t{q(.5):.0f}\t{q(.75):.0f}\t{q(.9):.0f}\t{q(.5) / ctl:.2f}\t{sum(x < .5 * ctl for x in v) / len(v):.2f}\t{sum(x > 2 * ctl for x in v) / len(v):.2f}")
v = sorted(x for n, x in D.items() if n.startswith("ctl"))
q = lambda f: v[int(f * (len(v) - 1))]
print(f"control\t{len(v)}\t{q(.1):.0f}\t{q(.25):.0f}\t{q(.5):.0f}\t{q(.75):.0f}\t{q(.9):.0f}\t1.00\t{sum(x < .5 * ctl for x in v) / len(v):.2f}\t{sum(x > 2 * ctl for x in v) / len(v):.2f}")
