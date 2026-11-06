#!/usr/bin/env python3
'''
Pre-placement tests for tails from make_tails.py. See docs/decisions/D16_tails.md. The output suggests; the agent decides.

For each tail: does its sequence already exist on the receiving scaffold or on other haplotypes of the same chromosome, are its chains near-identical
to other candidate arms (one shared copy budget), and what do their depths sum to?

Usage: tail_checks.py --tails tails.tsv --chains chains.fasta --target pre_tails_target.fasta --names names.tsv --class chain_class.tsv
         --hapDepth H --ploidy N --out tail_checks.tsv [--report tail_report.tsv --minimap2 minimap2 --threads 4]
--report gives interior depth, outside graph overlaps, and adds candidates left unused to the twin test.
'''
import argparse, collections, os, subprocess, tempfile
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in ("tails", "chains", "target", "names", "class", "out"): ap.add_argument("--" + x, required=True)
ap.add_argument("--hapDepth", type=float, required=True); ap.add_argument("--ploidy", type=int, required=True); ap.add_argument("--minimap2", default="minimap2")
ap.add_argument("--threads", type=int, default=4); ap.add_argument("--report", help="make_tails.py --report output")
ap.add_argument("--minCopy", type=float, default=0.3); ap.add_argument("--minId", type=float, default=0.9)
a = ap.parse_args()
names = {l.split("\t")[0]: l.split("\t")[1] for l in open(a.names) if l[0] != "#"}
depth = {l.split("\t")[0]: l.split("\t")[2] for l in open(vars(a)["class"]) if l[0] != "#"}
T = []
for l in open(a.tails):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t"); T.append(dict(end=s[0], group=s[1], path=[p[:-1] for p in s[2].split(",")], label=s[3], pairs=s[4], bp=int(s[5])))
base = lambda c: c.split(".")[0]
arms = {base(c) for t in T for c in t["path"]}
if a.report:
    for l in open(a.report):
        if l[0] == "#": continue
        s = l.rstrip("\n").split("\t")
        if s[3] != "nan": depth[s[0]] = s[3]
        if s[7] in ("tail", "candidate_unused"): arms.add(s[0])
with tempfile.TemporaryDirectory() as d:
    q = d + "/arms.fa"
    with open(q, "w") as o:
        subprocess.run(["samtools", "faidx", a.chains] + sorted(arms), stdout=o, check=True)
    tg = subprocess.run([a.minimap2, "-x", "asm20", "-c", "-t", str(a.threads), "-N", "20", "--secondary=yes", a.target, q], capture_output=True, text=True, check=True).stdout
    ar = subprocess.run([a.minimap2, "-x", "asm20", "-c", "-t", str(a.threads), "-N", "20", "--secondary=yes", "-D", q, q], capture_output=True, text=True, check=True).stdout
def cover(text):
    C = collections.defaultdict(lambda: collections.defaultdict(list)); L = {}
    for l in text.splitlines():
        s = l.split("\t"); L[s[0]] = int(s[1])
        if s[0] == s[5] or s[0].split("_", 1)[-1] == s[5].split("_", 1)[-1] or int(s[10]) < 2000 or int(s[9]) / int(s[10]) < a.minId: continue
        C[s[0]][s[5]].append((int(s[2]), int(s[3]), int(s[9]) / int(s[10])))
    out = collections.defaultdict(dict)
    for q, v in C.items():
        for t, iv in v.items():
            iv.sort(); c = e = 0
            for x, y, _ in iv:
                if y > e: c += y - max(x, e); e = y
            out[q][t] = (c / L[q], max(i for _, _, i in iv))
    return out
CT, CA = cover(tg), cover(ar)
num = lambda v: float(v) if v not in ("nan", "") else 0.0
with open(a.out, "w") as o:
    o.write("#end\tgroup\tscaffold\tpath\tlabel\tpairs\tbp\tdepth\tcopy_on_receiver\thaps_of_chromosome_with_copy\ttwins\ttwin_depth_sum\tother_copies\tsuggest\n")
    for t in T:
        sc = names.get(t["group"], t["group"]); chrom = sc.rsplit("_", 1)[0]; rec = 0.0; haps = set(); other = []; twins = {}; d = 0.0
        for c in t["path"]:
            b = base(c); d = max(d, num(depth.get(b, "0")))
            for tt, (cv, i) in CT.get(b, {}).items():
                if cv < a.minCopy: continue
                if tt == sc: rec = max(rec, cv)
                elif tt.rsplit("_", 1)[0] == chrom: haps.add(tt)
                else: other.append(f"{tt}:{cv:.2f}@{i:.3f}")
                if tt == sc: haps.add(tt)
            for tt, (cv, i) in CA.get(b, {}).items():
                if cv >= 0.5 and i >= 0.995: twins[tt] = (cv, i)
        tsum = d + sum(num(depth.get(x.split("_", 1)[-1], "0")) for x in twins); why = []
        if rec >= a.minCopy: why.append("copy_on_receiver")
        if len(haps) >= a.ploidy: why.append("no_free_slot")
        if twins and tsum < 1.5 * a.hapDepth: why.append("twins_sum_to_one_copy")
        if t["label"] == "arbitrary": why.append("assignment_arbitrary")
        o.write("\t".join(str(x) for x in (t["end"], t["group"], sc, ",".join(t["path"]), t["label"], t["pairs"], t["bp"], round(d, 1), round(rec, 2), ";".join(sorted(haps)) or ".", ";".join(f"{k}:{v[0]:.2f}@{v[1]:.3f}" for k, v in twins.items()) or ".", round(tsum, 1), ";".join(other) or ".", "pass" if not why else ",".join(why))) + "\n")
print("wrote", a.out)
