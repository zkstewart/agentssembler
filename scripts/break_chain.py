#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Apply chain breaks without remapping. Each line of the breaks file: chain position (1-based last base of
the left part). The left part keeps its name and coordinates; the right part is <chain>b and gets back
the graph overlap trimmed from its first unitig. A chain may have several lines: positions are all on the
unbroken chain, and the parts from left to right are <chain>, <chain>b, <chain>c, ... Writes a new AGP and
FASTA, and lifts Hi-C pairs and links.
Usage: break_chain.py breaks.txt in.agp in.pairs.tsv.gz outdir unitigs.fasta
'''
import sys, gzip, collections, os
brk, agp, pairs, outdir, unitigs = sys.argv[1:6]
B = collections.defaultdict(set)
for l in open(brk):
    if l[0] == "#" or not l.strip(): continue
    c, p = l.split()[:2]; B[c].add(int(p))
B = {c: sorted(ps) for c, ps in B.items()}
SUF = "bcdefghijklmnopqrstuvwxyz"
for c, ps in B.items():
    if len(ps) > len(SUF): sys.exit(f"{c}: more than {len(SUF)} breaks")
U = {}; n = None
for l in open(unitigs):
    if l[0] == ">": n = l[1:].split()[0]; U[n] = []
    else: U[n].append(l.strip())
for k in U: U[k] = "".join(U[k])
rc = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
A = collections.defaultdict(list); order = []
for l in open(agp):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if s[0] not in A: order.append(s[0])
    A[s[0]].append(s)
shift = {}; new = collections.OrderedDict()
for c in order:
    if c not in B: new[c] = A[c]; continue
    rows = A[c]; ends = {int(s[2]) for s in rows}; L = int(rows[-1][2])
    for p in B[c]:
        if p not in ends or p >= L: sys.exit(f"{c}: {p} is not a junction")
    cuts = [0] + B[c] + [L]; shift[c] = []
    for k in range(len(cuts) - 1):
        a, b = cuts[k], cuts[k + 1]
        part = [s for s in rows if a < int(s[2]) <= b]
        if k == 0: new[c] = part; continue
        name = c + SUF[k - 1]
        f = part[0][:]; ulen = len(U[f[5]])
        if f[8] == "+": restored = int(f[6]) - 1; f[6] = "1"
        else: restored = ulen - int(f[7]); f[7] = str(ulen)
        out = []; pos = 0
        for i, s in enumerate([f] + part[1:]):
            m = int(s[7]) - int(s[6]) + 1
            out.append([name, str(pos + 1), str(pos + m), str(i + 1), "W", s[5], s[6], s[7], s[8]]); pos += m
        new[name] = out; shift[c].append((a, name, restored - a))
        print(f"{c}: {name} from {a + 1}, {pos} bp (restored {restored} bp of {f[5]})")
    print(f"{c}: left {cuts[1]} bp")
os.makedirs(outdir, exist_ok=True)
with open(f"{outdir}/chains.agp", "w") as o:
    o.write("##agp-version 2.1\n")
    for c, rows in new.items():
        for s in rows: o.write("\t".join(s) + "\n")
with open(f"{outdir}/chains.fasta", "w") as o:
    for c, rows in new.items():
        q = []
        for s in rows:
            x = U[s[5]][int(s[6]) - 1:int(s[7])]
            q.append(x if s[8] == "+" else x.translate(rc)[::-1])
        o.write(f">{c}\n{''.join(q)}\n")
def mv(c, p):
    for a, name, d in reversed(shift.get(c, ())):
        if p > a: return name, p + d
    return c, p
links = collections.Counter(); dep = collections.Counter()
with gzip.open(pairs, "rt") as fh, gzip.open(f"{outdir}/pairs.tsv.gz", "wt", compresslevel=1) as o:
    for l in fh:
        x, u, y, v, *rest = l.rstrip("\n").split("\t"); u = int(u); v = int(v)
        x, u = mv(x, u); y, v = mv(y, v)
        o.write("\t".join([x, str(u), y, str(v)] + rest) + "\n")
        dep[x] += 1; dep[y] += 1
        if x != y: links[tuple(sorted((x, y)))] += 1
with open(f"{outdir}/links.tsv", "w") as o:
    for (x, y), m in links.items(): o.write(f"{x}\t{y}\t{m}\n")
    for c, m in dep.items(): o.write(f"{c}\t.\t{m}\n")
