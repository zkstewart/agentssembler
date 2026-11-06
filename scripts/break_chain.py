#!/usr/bin/env python3
'''
Apply chain breaks without remapping. Each line of the breaks file: chain position (1-based last base of
the left part). The left part keeps its name and coordinates; the right part is <chain>b and gets back
the graph overlap trimmed from its first unitig. Writes a new AGP and FASTA, and lifts Hi-C pairs and links.
Usage: break_chain.py breaks.txt in.agp in.pairs.tsv.gz outdir unitigs.fasta
'''
import sys, gzip, collections, os
brk, agp, pairs, outdir, unitigs = sys.argv[1:6]
B = {}
for l in open(brk):
    if l[0] == "#" or not l.strip(): continue
    c, p = l.split()[:2]; B[c] = int(p)
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
    left = [s for s in A[c] if int(s[2]) <= B[c]]; right = [s for s in A[c] if int(s[1]) > B[c]]
    if not left or not right or int(left[-1][2]) != B[c]: sys.exit(f"{c}: {B[c]} is not a junction")
    f = right[0][:]
    ulen = len(U[f[5]])
    if f[8] == "+": restored = int(f[6]) - 1; f[6] = "1"
    else: restored = ulen - int(f[7]); f[7] = str(ulen)
    r = [f] + right[1:]
    out = []; pos = 0
    for i, s in enumerate(r):
        m = int(s[7]) - int(s[6]) + 1
        out.append([c + "b", str(pos + 1), str(pos + m), str(i + 1), "W", s[5], s[6], s[7], s[8]]); pos += m
    new[c] = left; new[c + "b"] = out; shift[c] = (B[c], c + "b", restored - B[c])
    print(f"{c}: left {B[c]} bp, right {pos} bp (restored {restored} bp of {f[5]})")
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
    if c in shift and p > shift[c][0]: return shift[c][1], p + shift[c][2]
    return c, p
links = collections.Counter(); dep = collections.Counter()
with gzip.open(pairs, "rt") as fh, gzip.open(f"{outdir}/pairs.tsv.gz", "wt", compresslevel=1) as o:
    for l in fh:
        x, u, y, v = l.split("\t"); u = int(u); v = int(v)
        x, u = mv(x, u); y, v = mv(y, v)
        o.write(f"{x}\t{u}\t{y}\t{v}\n")
        dep[x] += 1; dep[y] += 1
        if x != y: links[tuple(sorted((x, y)))] += 1
with open(f"{outdir}/links.tsv", "w") as o:
    for (x, y), m in links.items(): o.write(f"{x}\t{y}\t{m}\n")
    for c, m in dep.items(): o.write(f"{c}\t.\t{m}\n")
