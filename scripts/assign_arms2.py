#!/usr/bin/env python3
'''
Diploid bridge choice for junctions adjacent in a chain order. Replaces "assign" in assign_arms.py, whose rule rejected any
path holding a chain the other haplotype used, so a path through a true 2x node always counted as reused and the 1x arm
beside that node was given to both haplotypes.

Each chain has a copy target from its depth outside the graph overlaps at both ends: 2 when depth >= --twoX, else 1. A chain
with no sequence outside its overlaps carries no weight. Junctions sharing candidate chains are solved together. A solution
scores + interior bp for a chain used as often as its target and - interior bp per use above the target. Ties go to the most
observations, counting reads through the path's joins and Hi-C pairs between its chains and the group's own chains.

A junction's label records what settled its arm:
  only        one candidate path
  coverage    every other path scores lower, alone and with a change at one other junction
  supported   some alternatives score the same, including exchanges of arms, and observations exceed each by >= --ratio and >= --minObs
  arbitrary   some alternative scores the same and observations do not separate them

Interior depth for chains in candidate paths is written to --depth; chains that were candidates but stay unused go to --unused.

Usage: assign_arms2.py --order order.tsv --bridges graph.bridges.tsv --ranked choices_base.tsv --groups groups.tsv
         --links links.tsv --agp chains.agp --gfa unitigs.gfa --paf ont_vs_chains.paf.gz --out choices.tsv
         --report report.tsv --depth depth.tsv --unused unused.tsv
'''
import argparse, collections, gzip, itertools
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in ("order", "bridges", "ranked", "groups", "links", "agp", "gfa", "paf", "out", "report", "depth", "unused"):
    ap.add_argument("--" + x, required=True)
ap.add_argument("--twoX", type=float, required=True, help="interior depth at or above this is 2x: 1.5 x haplotype depth")
ap.add_argument("--minAllele", type=float, required=True, help="alternative chain at this depth or above is a real allele")
ap.add_argument("--maxPair", type=float, required=True, help="two alternatives whose depths sum to at most this are two alleles of one locus")
ap.add_argument("--ratio", type=float, default=2.0); ap.add_argument("--minObs", type=int, default=5)
ap.add_argument("--maxCombos", type=int, default=200000)
a = ap.parse_args()
fl = {"+": "-", "-": "+"}
base = lambda c: c.split(".g")[0]
bend = lambda e: base(e[:-1]) + e[-1]
def rev(path): return [z[:-1] + fl[z[-1]] for z in reversed(path)]

gr = {}; L = {}
for l in open(a.groups):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); L[s[0]] = int(s[1]); gr[s[0]] = s[2]
members = collections.defaultdict(set)
for c, g in gr.items():
    if g != "unassigned" and ".g" not in c: members[g].add(c)
lk = collections.defaultdict(dict)
for l in open(a.links):
    x, y, n = l.rstrip().split("\t")
    if y != ".": lk[x][y] = int(n); lk[y][x] = int(n)

# candidates per copy-level junction, each written in the junction's own direction
cand = collections.defaultdict(list)
for l in open(a.bridges):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if s[2] == "yes": cand[(s[0], s[1])].append(tuple(s[3].split(",")))
# read counts per base junction and path
rank = {}
def put(x, y, pth, c):
    rank[(x, y, pth)] = c; rank[(y, x, tuple(rev(list(pth))))] = (c[1], c[0], c[2])
for l in open(a.ranked):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    put(s[0], s[1], tuple(s[3].split(",")), (int(s[4]), int(s[5]), int(s[6])))
    if s[8] != ".":
        for alt in s[8].split(";"):
            p, c = alt.rsplit(":", 1); c = c.split("/")
            put(s[0], s[1], tuple(p.split(",")), (int(c[0]), int(c[1]), int(c[2])))

# junctions adjacent in the order
J = {}
for rows in [[l.rstrip("\n").split("\t") for l in open(a.order) if l[0] != "#"]]:
    for r, d in zip(rows, rows[1:]):
        if r[0] != d[0]: continue
        x, y = r[2] + r[3], d[2] + fl[d[3]]
        if (x, y) in cand: J[(x, y)] = (r[0], cand[(x, y)])
        elif (y, x) in cand: J[(y, x)] = (r[0], cand[(y, x)])
print(len(J), "adjacent junctions with candidates,", sum(len(v[1]) for v in J.values()), "paths")

# interior depth
need = {z[:-1] for g, ps in J.values() for p in ps for z in p}
parts = collections.defaultdict(list)
for l in open(a.agp):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if s[4] == "W": parts[s[0]].append((s[5], s[8]))
ov = collections.defaultdict(int)
for l in open(a.gfa):
    if l[0] != "L": continue
    s = l.rstrip("\n").split("\t"); o = int(s[5].rstrip("M")) if s[5] != "*" else 0
    ov[(s[1], s[2])] = max(ov[(s[1], s[2])], o); k = (s[3], fl[s[4]]); ov[k] = max(ov[k], o)
iv = {}
for c in need:
    u, o = parts[c][-1]; right = ov[(u, o)]
    u, o = parts[c][0]; left = ov[(u, fl[o])]
    iv[c] = (left, L[c] - right)
cov = collections.Counter()
for l in gzip.open(a.paf, "rt"):
    s = l.split("\t", 12); c = s[5]
    if c not in iv or "tp:A:P" not in l: continue
    lo, hi = max(iv[c][0], int(s[7])), min(iv[c][1], int(s[8]))
    if hi > lo: cov[c] += hi - lo
W = {}; T = {}; D = {}
for c in need:
    w = max(iv[c][1] - iv[c][0], 0); W[c] = w; D[c] = cov[c] / w if w else float("nan")
    T[c] = 2 if w and D[c] >= a.twoX else 1
# a chain near 2x with a real alternative in the same junction is the larger of two alleles, not one 2x chain
sib = collections.defaultdict(set)
for g, ps in J.values():
    sets = [{z[:-1] for z in p} for p in ps]
    for c in {x for t in sets for x in t}:
        for t in sets:
            if c not in t:
                keep = any(c in u and not (u & t) for u in sets)
                if keep: sib[c] |= {x for x in t if not any(c in u and x in u for u in sets)}
for c in need:
    if T[c] == 2 and any(W[x] and D[x] >= a.minAllele and D[c] + D[x] <= a.maxPair for x in sib[c]): T[c] = 1
with open(a.depth, "w") as o:
    o.write("#chain\tlength\tinterior\tdepth\ttarget\n")
    for c in sorted(need): o.write(f"{c}\t{L[c]}\t{W[c]}\t{D[c]:.1f}\t{T[c] if W[c] else '.'}\n")

# junctions sharing a candidate chain are solved together
parent = {j: j for j in J}
def find(x):
    while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
    return x
owner = {}
for j, (g, ps) in J.items():
    for p in ps:
        for z in p:
            c = z[:-1]
            if c in owner: parent[find(j)] = find(owner[c])
            else: owner[c] = j
comps = collections.defaultdict(list)
for j in J: comps[find(j)].append(j)

def obs(j, p):
    x, y = j; g = J[j][0]; r = rank.get((bend(x), bend(y), p), (0, 0, 0))
    return r[0] + r[1] + 5 * r[2] + sum(lk[z[:-1]].get(m, 0) for z in p for m in members[g])
def score(sol):
    n = collections.Counter(z[:-1] for p in sol.values() for z in p); s = 0
    for c, k in n.items():
        s += W[c] if k <= T[c] and k == T[c] else (-W[c] * (k - T[c]) if k > T[c] else 0)
    return s
def key(sol): return (score(sol), sum(obs(j, p) for j, p in sol.items()))

chosen = {}; stats = collections.Counter(); searched = collections.Counter()
for root, js in comps.items():
    js = sorted(js); opts = [J[j][1] for j in js]; n = 1
    for o in opts: n *= len(o)
    if n <= a.maxCombos:
        searched["exhaustive"] += 1
        best = max((dict(zip(js, c)) for c in itertools.product(*opts)), key=key)
    else:                                        # local search from each junction's first path
        searched["local"] += 1
        best = {j: o[0] for j, o in zip(js, opts)}; moved = True
        while moved:
            moved = False
            for j, o in zip(js, opts):
                for p in o:
                    t = dict(best); t[j] = p
                    if key(t) > key(best): best = t; moved = True
    chosen.update(best)
print("components:", dict(searched))

# label each junction: try its other paths alone and together with a change at one other junction, so that an exchange of
# arms between the haplotypes, which coverage cannot see, counts as a tie
compof = {j: find(j) for j in J}
def decided(o0, o1): return o0 >= a.minObs + o1 and o0 >= a.ratio * max(o1, 1)
labels = {}
for j, p in chosen.items():
    ps = J[j][1]
    if len(ps) == 1: labels[j] = "only"; continue
    here = score(chosen); ties = []
    others = [k for k in J if k != j and compof[k] == compof[j]]
    for q in ps:
        if q == p: continue
        t = dict(chosen); t[j] = q
        if score(t) >= here: ties.append((obs(j, p), obs(j, q)))
        for k in others:
            for q2 in J[k][1]:
                if q2 == chosen[k]: continue
                t2 = dict(t); t2[k] = q2
                if score(t2) >= here: ties.append((obs(j, p) + obs(k, chosen[k]), obs(j, q) + obs(k, q2)))
    if not ties: labels[j] = "coverage"
    elif all(decided(x, y) for x, y in ties): labels[j] = "supported"
    else: labels[j] = "arbitrary"
    stats[labels[j]] += 1
stats["only"] = sum(1 for v in labels.values() if v == "only")

with open(a.out, "w") as o, open(a.report, "w") as q:
    o.write("#end_a\tend_b\tchosen\tpath\tleft_reads\tright_reads\tfull_reads\tpath_bp\talternatives\n")
    q.write("#group\tend_a\tend_b\thow\tpath\tleft_reads\tright_reads\tfull_reads\tobservations\tcandidates\n")
    for j, p in sorted(chosen.items(), key=lambda kv: (J[kv[0]][0], kv[0])):
        x, y = j; g, ps = J[j]; r = rank.get((bend(x), bend(y), p), (0, 0, 0))
        alts = ";".join(",".join(t) for t in ps if t != p) or "."
        o.write(f"{x}\t{y}\tyes\t{','.join(p)}\t{r[0]}\t{r[1]}\t{r[2]}\t{sum(L[z[:-1]] for z in p)}\t{labels[j]}:{alts}\n")
        q.write(f"{g}\t{x}\t{y}\t{labels[j]}\t{','.join(p)}\t{r[0]}\t{r[1]}\t{r[2]}\t{obs(j, p)}\t{len(ps)}\n")

n = collections.Counter(z[:-1] for p in chosen.values() for z in p)
with open(a.unused, "w") as o:
    o.write("#chain\tlength\tinterior\tdepth\ttarget\tused\tstate\n")
    tot = collections.Counter()
    for c in sorted(need, key=lambda c: -W[c]):
        k = n[c]
        st = "no_interior" if not W[c] else "ok" if k == T[c] else "unused" if k == 0 else "under" if k < T[c] else "over"
        tot[st] += W[c]; tot[st + "_n"] += 1
        o.write(f"{c}\t{L[c]}\t{W[c]}\t{D[c]:.1f}\t{T[c]}\t{k}\t{st}\n")
print("labels", dict(stats)); print("interior bp and chains by state", dict(tot))
