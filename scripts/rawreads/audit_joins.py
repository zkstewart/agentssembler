#!/usr/bin/env python3
'''
Test read-derived joins against reads that were not used to make them (normally raw reads). One row per join; a suggested verdict
that the agent must read, not obey. See docs/decisions/D11_gap_closures.md.

Input: readfill reports (label=report.tsv), the AGP holding the joins, a PAF from `minimap2 -c --eqx` against that AGP's sequences, and
for the overlap test the FASTA each report's gap IDs refer to.

Columns: reads through the join with at least --flank bp each side, the same at MAPQ >= 20, those with no indel above --maxIndel within 2 kb,
reads breaking at the join, identity at the join and in the flanks, read/target length ratio, alignments per kb within 20 kb (a tandem
array shows hundreds of times normal depth), and for overlap joins the share of the trimmed sequence with a copy in the kept flank.

Each --report is label=report.tsv or label=report.tsv,pre.fasta, where pre.fasta is the sequence that report's gap IDs refer to.
Usage: audit_joins.py --agp final.agp --paf audit.paf.gz --report corrected=a.report.tsv,pre1.fasta --report raw=b.report.tsv,pre2.fasta
         [--minimap2 minimap2] --hapDepth H --out audit.tsv
'''
import argparse, collections, csv, gzip, os, re, statistics as st, subprocess, tempfile
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--agp", required=True); ap.add_argument("--paf", required=True, nargs="+"); ap.add_argument("--report", action="append", required=True)
ap.add_argument("--minimap2", default="minimap2"); ap.add_argument("--hapDepth", type=float, required=True)
ap.add_argument("--flank", type=int, default=3000); ap.add_argument("--maxIndel", type=int, default=50); ap.add_argument("--out", required=True)
a = ap.parse_args()
rows = collections.defaultdict(list)
for l in open(a.agp):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t"); rows[s[0]].append(s)
J = []
for x in a.report:
    src, p = x.split("=", 1); p, _, pre = p.partition(",")
    for r in csv.DictReader(open(p), delimiter="\t"):
        if r["filled"] != "yes": continue
        t = r["#gap_id"].split(":")[0]; lc, rc = r["left"][:-1], r["right"][:-1]; v = rows[t]; hit = None
        for i in range(len(v) - 1):
            b, c = v[i], v[i + 1]; d = v[i + 2] if i + 2 < len(v) else None
            if b[4] == "W" == c[4] and {b[5], c[5]} == {lc, rc}: hit = (int(b[2]), int(c[1]) - 1); break
            if d is not None and c[5].startswith("rf_") and {b[5], d[5]} == {lc, rc}: hit = (int(b[2]), int(d[1]) - 1); break
        J.append(dict(id=r["#gap_id"], src=src, method=r["method"], t=t, pos=hit, trimL=int(r["trim_left"]), trimR=int(r["trim_right"]), fill=int(r["fill_bp"]), reads=[], aln=0, pre=pre))
A = collections.defaultdict(list)
for j in J:
    if j["pos"]: A[j["t"]].append(j)
for path in a.paf:
    for l in gzip.open(path, "rt"):
        s = l.rstrip("\n").split("\t"); t = s[5]
        if t not in A: continue
        ts, te = int(s[7]), int(s[8])
        for j in A[t]:
            x0, x1 = j["pos"]
            if ts < x1 + 20000 and te > x0 - 20000: j["aln"] += min(te, x1 + 20000) - max(ts, x0 - 20000)
            if not (ts < x1 + 5000 and te > x0 - 5000): continue
            ql, qs, qe, mq = int(s[1]), int(s[2]), int(s[3]), int(s[11])
            lo, hi = (qs, ql - qe) if s[4] == "+" else (ql - qe, qs)
            cg = next((y[5:] for y in s[12:] if y.startswith("cg:Z:")), "")
            x = ts; eqJ = nJ = eqF = nF = rb = tb = mx = 0
            for n, o in re.findall(r"(\d+)([=XID])", cg):
                n = int(n)
                if o in "=X":
                    e = x + n; oj = max(0, min(e, x1 + 500) - max(x, x0 - 500)); of = max(0, min(e, x0 - 500) - max(x, x0 - 5000)) + max(0, min(e, x1 + 5000) - max(x, x1 + 500)); ow = max(0, min(e, x1 + a.flank) - max(x, x0 - a.flank))
                    nJ += oj; nF += of; rb += ow; tb += ow
                    if o == "=": eqJ += oj; eqF += of
                    x = e
                else:
                    inw = x0 - a.flank <= x < x1 + a.flank
                    if x0 - 2000 <= x < x1 + 2000: mx = max(mx, n)
                    if x0 - 500 <= x < x1 + 500: nJ += n
                    elif x0 - 5000 <= x < x1 + 5000: nF += n
                    if o == "D":
                        if inw: tb += n
                        x += n
                    elif inw: rb += n
            j["reads"].append(dict(thru=ts <= x0 - a.flank and te >= x1 + a.flank, mq=mq, mx=mx, br=(x0 - 500 <= te <= x1 + 500 and hi >= 2000 and te - ts >= 2000) or (x0 - 500 <= ts <= x1 + 500 and lo >= 2000 and te - ts >= 2000),
                idJ=eqJ / nJ if nJ >= 200 else None, idF=eqF / nF if nF >= 200 else None, ratio=rb / tb if tb else None))
def faidx(fa, t, x, y):
    return "".join(subprocess.run(["samtools", "faidx", fa, f"{t}:{x + 1}-{y}"], capture_output=True, text=True, check=True).stdout.split("\n")[1:])
def overlapTest(j):
    if not j["pre"] or j["method"] != "reads_overlap": return None
    t, p = j["id"].split(":"); p = int(p) - 1; tr = j["trimR"] + j["trimL"]
    if tr == 0: return None
    with tempfile.TemporaryDirectory() as d:
        open(d + "/L.fa", "w").write(">L\n" + faidx(j["pre"], t, max(0, p - tr - 20000), p) + "\n")
        open(d + "/R.fa", "w").write(">R\n" + faidx(j["pre"], t, p + 100, p + 100 + j["trimR"] if j["trimR"] else p + 100 + 1) + "\n")
        o = subprocess.run([a.minimap2, "-x", "asm20", "-c", "--secondary=no", d + "/L.fa", d + "/R.fa"], capture_output=True, text=True).stdout.splitlines()
    if not o: return (0.0, 0.0)
    b = max((x.split("\t") for x in o), key=lambda s: int(s[10]))
    return (int(b[3]) - int(b[2])) / max(1, j["trimR"]), int(b[9]) / int(b[10])
med = lambda R, k: round(st.median(x[k] for x in R if x[k] is not None), 4) if any(x[k] is not None for x in R) else None
with open(a.out, "w") as o:
    o.write("#join\tsource\tmethod\tscaffold\tposition\tfill_bp\ttrim_left\ttrim_right\tthrough\tthrough_q20\tthrough_q20_clean\tbreaks\tid_join\tid_flank\tlen_ratio\tdepth_40kb\ttrim_cover\ttrim_identity\tsuggest\n")
    for j in J:
        R = j["reads"]; T = [x for x in R if x["thru"]]; Q = [x for x in T if x["mq"] >= 20]; C = [x for x in Q if x["mx"] <= a.maxIndel]; br = sum(x["br"] for x in R)
        depth = (j["aln"] / 40000) if j["pos"] else 0; ov = overlapTest(j); ratio = med(Q, "ratio"); why = []
        if not j["pos"]: why.append("join_not_found")
        if len(C) < 2: why.append("under_2_clean_q20_reads")
        if depth > 4 * a.hapDepth: why.append("array_depth")
        if ov and (ov[0] < 0.9 or ov[1] < 0.99): why.append("trim_has_no_copy_in_kept_flank")
        if ratio is not None and abs(ratio - 1) > 0.02: why.append("length_ratio_off")
        if br > max(3, len(T)): why.append("reads_break_here")
        idJ, idF = med(Q, "idJ"), med(Q, "idF")
        if idJ is not None and idJ < 0.95: why.append("low_identity_at_join")
        p = j["pos"][0] if j["pos"] else "."
        o.write("\t".join(str(x) for x in (j["id"], j["src"], j["method"], j["t"], p, j["fill"], j["trimL"], j["trimR"], len(T), len(Q), len(C), br, idJ, idF, ratio, round(depth, 1), round(ov[0], 3) if ov else ".", round(ov[1], 4) if ov else ".", "pass" if not why else ",".join(why))) + "\n")
print(f"{len(J)} joins written to {a.out}")
