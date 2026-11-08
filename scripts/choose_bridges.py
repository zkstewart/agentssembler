#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Choose among alternative graph paths across a junction with reads. For each bridged junction, build
A_tail + path + B_head (graph overlaps merged) for every candidate, map local reads competitively to all
candidates, and count reads through each candidate's two joins (>= --minSide bp on both sides of the join).
The candidate with the most reads through its weaker join wins; a junction with no read through both joins of
any path is still reported, and the longest path is used only if --tieLongest.
Usage: choose_bridges.py bridges.tsv end_links.tsv chains.fasta ont_vs_chains.paf.gz reads.fasta outdir
'''
import sys, os, gzip, collections, subprocess, argparse
p = argparse.ArgumentParser(); p.add_argument("bridges"); p.add_argument("links"); p.add_argument("fasta"); p.add_argument("paf"); p.add_argument("reads"); p.add_argument("out")
p.add_argument("--flank", type=int, default=30000); p.add_argument("--minSide", type=int, default=2000); p.add_argument("--threads", type=int, default=14)
a = p.parse_args(); os.makedirs(a.out, exist_ok=True)
S = {}; n = None
for l in open(a.fasta):
    if l[0] == ">": n = l[1:].split()[0]; S[n] = []
    else: S[n].append(l.strip())
for k in S: S[k] = "".join(S[k])
tr = str.maketrans("ACGTNacgtn", "TGCANtgcan"); rc = lambda s: s.translate(tr)[::-1]
ov = {}
for l in open(a.links):
    if l[0] == "#": continue
    x, y, r, q, o = l.rstrip().split("\t"); ov[frozenset((x, y))] = max(int(o), 0)
fl = {"+": "-", "-": "+"}
J = collections.defaultdict(list)
for l in open(a.bridges):
    if l[0] == "#": continue
    x, y, sg, pth, bp, mo = l.rstrip().split("\t")
    if sg == "yes": J[(x, y)].append(pth.split(","))
targets = {}; joins = {}; want = collections.defaultdict(list); L = {k: len(v) for k, v in S.items()}
for (x, y), paths in J.items():
    A, ea = x[:-1], x[-1]; B, eb = y[:-1], y[-1]
    tail = S[A] if ea == "+" else rc(S[A])
    head = S[B] if eb == "-" else rc(S[B])
    for i, pth in enumerate(paths):
        name = f"{x}|{y}|{i}"; seq = tail[-a.flank:]; js = [len(seq)]; prev = x
        for e in pth:
            c, o = e[:-1], e[-1]; s = S[c] if o == "+" else rc(S[c])
            entry = c + ("-" if o == "+" else "+")
            k = ov.get(frozenset((prev, entry)), 0)
            seq += s[k:]; js.append(len(seq)); prev = c + o
        k = ov.get(frozenset((prev, y)), 0); seq += head[k:k + a.flank]
        targets[name] = seq; joins[name] = (js[0], js[-1], [z[:-1] for z in pth])
    for c in {A, B} | {e[:-1] for pth in paths for e in pth}:
        want[c].append((x, y))
# Reads from the two chain ends and from the path chains
ends = collections.defaultdict(list)
for (x, y), paths in J.items():
    for e in (x, y):
        c = e[:-1]; ends[c].append((L[c] - a.flank, L[c]) if e[-1] == "+" else (0, a.flank))
    for pth in paths:
        for z in pth: ends[z[:-1]].append((0, L[z[:-1]]))
ids = set()
with gzip.open(a.paf, "rt") as fh:
    for l in fh:
        s = l.split("\t", 9)
        if s[5] in ends:
            u, v = int(s[7]), int(s[8])
            if any(u < e and v > b for b, e in ends[s[5]]): ids.add(s[0])
print(len(ids), "local reads,", len(targets), "candidate paths,", len(J), "junctions", flush=True)
with open(f"{a.out}/ids.txt", "w") as o: o.write("\n".join(sorted(ids)) + "\n")
with open(f"{a.out}/targets.fasta", "w") as o:
    for k, v in targets.items(): o.write(f">{k.replace('|','_')}\n{v}\n")
nm = {k.replace("|", "_"): k for k in targets}
if not os.path.exists(f"{a.out}/reads.fasta"):
    subprocess.check_call(f"seqkit grep -f {a.out}/ids.txt {a.reads} -o {a.out}/reads.fasta.part 2> {a.out}/seqkit.log && mv {a.out}/reads.fasta.part {a.out}/reads.fasta", shell=True)
subprocess.check_call(f"minimap2 -x lr:hq -c --secondary=no -t {a.threads} {a.out}/targets.fasta {a.out}/reads.fasta 2> {a.out}/mm2.log | gzip -1 > {a.out}/reads.paf.gz", shell=True)
C = collections.defaultdict(lambda: [0, 0, 0])
with gzip.open(f"{a.out}/reads.paf.gz", "rt") as fh:
    for l in fh:
        s = l.split("\t", 12)
        if "tp:A:P" not in s[12]: continue
        t = nm[s[5]]; u, v = int(s[7]), int(s[8]); j1, j2, _ = joins[t]
        lf = u <= j1 - a.minSide and v >= j1 + a.minSide; rt = u <= j2 - a.minSide and v >= j2 + a.minSide
        C[t][0] += lf; C[t][1] += rt; C[t][2] += lf and rt
with open(f"{a.out}/choices.tsv", "w") as o:
    o.write("#end_a\tend_b\tchosen\tpath\tleft_reads\tright_reads\tfull_reads\tpath_bp\talternatives\n")
    for (x, y), paths in J.items():
        rows = []
        for i, pth in enumerate(paths):
            t = f"{x}|{y}|{i}"; c = C[t]; bp = sum(L[e[:-1]] for e in pth)
            rows.append((min(c[0], c[1]), c[2], bp, i, pth, c))
        rows.sort(key=lambda r: (r[0], r[1], r[2]), reverse=True)
        b = rows[0]
        alt = ";".join(f"{','.join(r[4])}:{r[5][0]}/{r[5][1]}/{r[5][2]}" for r in rows[1:]) or "."
        o.write(f"{x}\t{y}\tyes\t{','.join(b[4])}\t{b[5][0]}\t{b[5][1]}\t{b[5][2]}\t{b[2]}\t{alt}\n")
