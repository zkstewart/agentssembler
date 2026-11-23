#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Where do reads go when they run off a gap edge? For each gap, takes reads aligned up to a flank's end with at least --minTail unaligned bases beyond it, maps those tails
to the whole assembly including unplaced sequence, and lists destinations per edge. A destination shared by several reads is a candidate for what belongs in the gap. Tail destinations often lie
inside a chain, past an overlap with the flank, which chain_end_links.py does not count. Candidates only: test the graph path, the homolog and the Hi-C.
also_placed: reads whose tail aligns nearly as well (secondary hit, chaining score >= --altScore of the primary, over >= --altScore of its length)
to a placed scaffold as to the destination. When it is half the reads or more, an unplaced destination duplicates sequence already placed, and
its MAPQ is low for that reason (292: chain0285, a copy of placed chain0124). A low MAPQ with also_placed 0 is a near-identical unplaced pair
instead (chain0091/chain0092).
The reads may be FASTA or FASTQ, gzipped or not.
Usage: gap_tails.py final.agp reads.fasta all_sequence.fasta out.tsv --paf a.paf.gz [b.paf.gz] [--minimap2 minimap2] [--seqkit seqkit] [--preset lr:hq] [--threads 8]'''
import argparse, collections, gzip, os, subprocess, tempfile
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in "agp reads fasta out".split(): p.add_argument(x)
p.add_argument("--paf", nargs="+", required=True); p.add_argument("--minimap2", default="minimap2"); p.add_argument("--seqkit", default="seqkit"); p.add_argument("--preset", default="lr:hq", help="minimap2 preset for the reads"); p.add_argument("--threads", default="8")
p.add_argument("--minTail", type=int, default=2000); p.add_argument("--edge", type=int, default=3000); p.add_argument("--minQ", type=int, default=20)
p.add_argument("--altScore", type=float, default=0.9); a = p.parse_args()
G = collections.defaultdict(list)
for l in open(a.agp):
    if l[0] == "#": continue
    s = l.split("\t")
    if s[4] in "NU": G[s[0]].append((int(s[1]) - 1, int(s[2])))
T = {}
for f in a.paf:
    for l in gzip.open(f, "rt"):
        s = l.split("\t", 12); t = s[5]
        if t not in G: continue
        qs, qe, ql, ts, te, mq = int(s[2]), int(s[3]), int(s[1]), int(s[7]), int(s[8]), int(s[11])
        for g0, g1 in G[t]:
            if g1 <= ts < g1 + a.edge: side, beyond = "right", (qs if s[4] == "+" else ql - qe)     # aligned to the right flank, read runs on leftward, over the gap
            elif g0 - a.edge < te <= g0: side, beyond = "left", (ql - qe if s[4] == "+" else qs)
            else: continue
            if beyond >= a.minTail:
                k = (s[0], t, g0, side)
                if k not in T or mq > T[k][3]: T[k] = (qs, qe, ql, mq, s[4])
want = collections.defaultdict(list)
for k, v in T.items(): want[k[0]].append(k)
tmp = tempfile.NamedTemporaryFile("w", suffix=".fa", delete=False)
if want:   # seqkit reads FASTA or FASTQ, gzipped or not, and keeps the file's order
    ids = tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False); ids.write("".join(n + "\n" for n in want)); ids.close()
    G1 = subprocess.Popen([a.seqkit, "grep", "-j", str(a.threads), "-f", ids.name, a.reads], stdout=subprocess.PIPE)
    G2 = subprocess.Popen([a.seqkit, "fx2tab", "-i", "-Q"], stdin=G1.stdout, stdout=subprocess.PIPE, text=True); G1.stdout.close()
    for l in G2.stdout:
        name, s = l.rstrip("\n").split("\t")[:2]
        for k in want.get(name, ()):
            qs, qe, ql, mq, o = T[k]
            t = s[:qs] if (k[3] == "right") == (o == "+") else s[qe:]
            if len(t) >= a.minTail: tmp.write(f">{name}|{k[1]}|{k[2]}|{k[3]}|{mq}\n{t}\n")
    if G2.wait() or G1.wait(): raise SystemExit(f"seqkit failed reading {a.reads}")
    os.unlink(ids.name)
tmp.close()
paf = subprocess.run([a.minimap2, "-x", a.preset, "-t", a.threads, "--secondary=yes", "-N", "5", a.fasta, tmp.name], capture_output=True, text=True, check=True).stdout
os.unlink(tmp.name)
D = collections.defaultdict(lambda: collections.defaultdict(list)); N = collections.Counter()
for k, v in T.items():
    if v[3] >= a.minQ: N[k[1], k[2], k[3]] += 1
tag = lambda s, k: next((x[5:] for x in s[12:] if x.startswith(k)), None)
P, ALT = {}, collections.defaultdict(list)     # tail -> primary (target, aligned, score); tail -> secondaries on placed scaffolds
for l in paf.splitlines():
    s = l.rstrip("\n").split("\t")
    if tag(s, "tp:A:") == "P":
        if s[0] not in P or int(s[3]) - int(s[2]) > P[s[0]][1]: P[s[0]] = (s[5], int(s[3]) - int(s[2]), int(tag(s, "s1:i:") or 0))
    elif not s[5].startswith("unplaced_"): ALT[s[0]].append((int(s[3]) - int(s[2]), int(tag(s, "s1:i:") or 0)))
def alt(q):
    t, n, sc = P[q]
    return any(m >= a.altScore * n and x >= a.altScore * sc for m, x in ALT.get(q, ()))
for l in paf.splitlines():
    s = l.rstrip("\n").split("\t")
    if tag(s, "tp:A:") != "P": continue
    q = s[0]; r, t, g0, side, mq = q.split("|")
    if int(mq) >= a.minQ and int(s[3]) - int(s[2]) >= 1000:
        D[t, int(g0), side][s[5]].append((r, int(s[3]) - int(s[2]), int(s[7]), int(s[8]), s[4], int(s[11]), alt(q)))
with open(a.out, "w") as o:
    o.write("#scaffold\tgap_start\tedge\ttails_q20\tdestination\treads\taligned_bp\tdest_start\tdest_end\tstrands\tmapq_median\talso_placed\n")
    for e in sorted(N, key=lambda x: (x[0], x[1], x[2])):
        ds = sorted(D[e].items(), key=lambda kv: -len({x[0] for x in kv[1]}))
        if not ds: o.write(f"{e[0]}\t{e[1]}\t{e[2]}\t{N[e]}\t.\t0\t0\t.\t.\t.\t.\t0\n")
        for t, v in ds:
            m = sorted(x[5] for x in v)[len(v) // 2]
            o.write(f"{e[0]}\t{e[1]}\t{e[2]}\t{N[e]}\t{t}\t{len({x[0] for x in v})}\t{sum(x[1] for x in v)}\t{min(x[2] for x in v)}\t{max(x[3] for x in v)}\t{''.join(sorted({x[4] for x in v}))}\t{m}\t{len({x[0] for x in v if x[6]})}\n")
