#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
One-sided graph paths at gaps and scaffold ends. A bridge needs a path joining both flanks; an arm that leaves one flank and stops
before the other is left unplaced by every earlier step, although the graph says where it starts. This places such arms beside
the flank they link to, with the gap kept beyond them.

Candidates are simple paths through unassigned chains from each gap end, and from each scaffold end without an outward telomere,
over graph links with overlap >= --minOverlap. A path ends before any chain that fails a test:
  copy limit   target from depth outside graph overlaps (2 at >= --twoX, else 1), less uses in the chosen bridges
  repeat       class repeat, or depth >= --maxDepth
  chromosome   Hi-C to another chromosome at >= --minPairs and >= --ratio times this chromosome's

Ends sharing candidate chains are solved together for the most interior bp, within the copy limits. Ties go to the most Hi-C pairs
between the path's chains and the group's own chains. Labels record what settled each tail:
  only        one non-empty candidate at this end, and no other end could take it
  supported   alternatives score the same, including an exchange with another end, and pairs exceed each by >= --ratio and >= --minObs
  arbitrary   alternatives score the same and pairs do not separate them; the haplotype is not known

Usage: make_tails.py --order order.tsv --junctions pre.junctions.tsv --placement pre.placement.tsv --links end_links.tsv
         --groups groups.tsv --chroms chroms.tsv --choices choices.tsv --telomeres telomeres.tsv --hic links.tsv
         --class chain_class.tsv --agp chains.agp --gfa unitigs.gfa --paf ont_vs_chains.paf.gz --out tails.tsv --report report.tsv
'''
import os, argparse, collections, gzip, itertools
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in ("order", "junctions", "placement", "links", "groups", "chroms", "choices", "telomeres", "hic", "class", "agp", "gfa", "paf", "out", "report"):
    ap.add_argument("--" + x, required=True)
ap.add_argument("--twoX", type=float, required=True, help="1.5 x haplotype depth"); ap.add_argument("--maxDepth", type=float, required=True, help="3 x haplotype depth")
ap.add_argument("--exclude", help="chains, first column, not to walk through: connectors have no sequence of their own to offer")
ap.add_argument("--minOverlap", type=int, default=5000); ap.add_argument("--maxHops", type=int, default=8)
ap.add_argument("--minPairs", type=int, default=20); ap.add_argument("--ratio", type=float, default=2.0)
ap.add_argument("--minObs", type=int, default=5); ap.add_argument("--maxCombos", type=int, default=200000)
a = ap.parse_args()
fl = {"+": "-", "-": "+"}
base = lambda c: c.split(".g")[0]

L = {}; G = {}
for l in open(a.groups):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); L[s[0]] = int(s[1]); G[s[0]] = s[2]
chrom = {}                                                   # group -> chromosome
for l in open(a.chroms):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); chrom[s[0]] = s[1]
members = collections.defaultdict(set); cmembers = collections.defaultdict(set)
for c, g in G.items():
    if g != "unassigned":
        cmembers[chrom[g]].add(base(c))
        if ".g" not in c: members[g].add(c)
hic = collections.defaultdict(dict)
for l in open(a.hic):
    x, y, n = l.rstrip().split("\t")
    if y != ".": hic[x][y] = int(n); hic[y][x] = int(n)
cls = {}
for l in open(a.__dict__["class"]):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); cls[s[0]] = s[-1]

adj = collections.defaultdict(dict)
excl = {l.split("\t")[0] for l in open(a.exclude)} if a.exclude and os.path.exists(a.exclude) else set()
for l in open(a.links):
    if l[0] == "#": continue
    x, y, r, q, o = l.rstrip("\n").split("\t"); o = int(o)
    if o >= a.minOverlap and x[:-1] not in excl and y[:-1] not in excl: adj[x][y] = o; adj[y][x] = o

# chains kept by the scaffold builder, and the ends that face a gap or a scaffold end
placed = {}
for l in open(a.placement):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); placed[s[0]] = s[4]
tel = {}
for l in open(a.telomeres):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); tel[s[0]] = s[4]
def tel_on(end):
    c = end[:-1]; v = tel.get(c, "--")
    return v[0 if end[-1] == "-" else 1] == "T"
ends = {}                                                    # end -> group
for l in open(a.junctions):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if s[4] == "hic":
        x, y = s[2], s[3][:-1] + fl[s[3][-1]]
        ends[x] = G[x[:-1]]; ends[y] = G[y[:-1]]
O = collections.defaultdict(list)
for l in open(a.order):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t")
        if placed.get(s[2]) == "placed": O[s[0]].append((s[2], s[3]))
for g, r in O.items():
    for e in (r[0][0] + fl[r[0][1]], r[-1][0] + r[-1][1]):
        if tel_on(e): print("telomere faces out, no tail:", e)
        else: ends[e] = g

# usage in the bridges that were chosen
used = collections.Counter()
for l in open(a.choices):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t")
        for z in s[3].split(","): used[z[:-1]] += 1

# interior depth of every unassigned chain
free = {c for c, g in G.items() if g == "unassigned"}
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
for c in free:
    if not parts[c]: continue
    u, o = parts[c][-1]; right = ov[(u, o)]
    u, o = parts[c][0]; left = ov[(u, fl[o])]
    iv[c] = (left, L[c] - right)
cov = collections.Counter()
for l in gzip.open(a.paf, "rt"):
    s = l.split("\t", 12); c = s[5]
    if c not in iv or "tp:A:P" not in l: continue
    lo, hi = max(iv[c][0], int(s[7])), min(iv[c][1], int(s[8]))
    if hi > lo: cov[c] += hi - lo
W = {}; D = {}; T = {}
for c in iv:
    w = max(iv[c][1] - iv[c][0], 0); W[c] = w; D[c] = cov[c] / w if w else float("nan")
    T[c] = 2 if w and D[c] >= a.twoX else 1
cap = {c: T[c] - used[c] if W[c] else 9 for c in iv}      # a chain with nothing outside its overlaps adds no sequence, so no limit

stops = collections.Counter(); stop_notes = {}
def allowed(c, g):
    if c not in iv: return "not_free"
    if cap[c] <= 0: return "copy"
    if cls.get(c) in ("repeat", "high") or (W[c] and D[c] >= a.maxDepth): return "repeat"
    here = sum(n for m, n in hic[c].items() if m in cmembers[chrom[g]])
    for ch, ms in cmembers.items():
        if ch == chrom[g]: continue
        there = sum(n for m, n in hic[c].items() if m in ms)
        if there >= a.minPairs and there >= a.ratio * max(here, 1): return "chromosome"
    return None

def paths(e, g):
    '''all simple paths from an end, as lists of oriented chains; the last chain must have sequence outside its overlaps'''
    out = []; stack = [(e, [])]
    while stack:
        cur, p = stack.pop()
        if p and W[p[-1][:-1]] > 0: out.append(tuple(p))
        if len(p) >= a.maxHops: continue
        for n in adj[cur]:
            c = n[:-1]
            if c in {z[:-1] for z in p}: continue
            why = allowed(c, g)
            if why:
                if why != "not_free": stops[why] += 1; stop_notes.setdefault((e, c), why)
                continue
            s = fl[n[-1]]                                    # entered at end n, leaves by the other end
            stack.append((c + s, p + [c + s]))
    return out

C = {e: paths(e, g) for e, g in ends.items()}
C = {e: v for e, v in C.items() if v}
print(len(ends), "ends examined;", len(C), "with candidates;", sum(len(v) for v in C.values()), "paths; stops", dict(stops))

parent = {e: e for e in C}
def find(x):
    while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
    return x
owner = {}
for e, ps in C.items():
    for p in ps:
        for z in p:
            c = z[:-1]
            if c in owner: parent[find(e)] = find(owner[c])
            else: owner[c] = e
comps = collections.defaultdict(list)
for e in C: comps[find(e)].append(e)

def feasible(sol):
    n = collections.Counter(z[:-1] for p in sol.values() for z in p)
    return all(k <= cap[c] for c, k in n.items())
def score(sol): return sum(W[z[:-1]] for p in sol.values() for z in set(p))
def obs(e, p):
    g = ends[e]
    return sum(hic[z[:-1]].get(m, 0) for z in p for m in members[g])
def key(sol): return (score(sol), sum(obs(e, p) for e, p in sol.items()))

chosen = {}; search = collections.Counter()
for root, es in comps.items():
    es = sorted(es); opts = [[()] + C[e] for e in es]; n = 1
    for o in opts: n *= len(o)
    best = None
    if n <= a.maxCombos:
        search["exhaustive"] += 1
        for combo in itertools.product(*opts):
            sol = dict(zip(es, combo))
            if feasible(sol) and (best is None or key(sol) > key(best)): best = sol
    else:
        search["local"] += 1
        best = {e: () for e in es}; moved = True
        while moved:
            moved = False
            for e, o in zip(es, opts):
                for p in o:
                    t = dict(best); t[e] = p
                    if feasible(t) and key(t) > key(best): best = t; moved = True
    chosen.update(best)
print("components", dict(search))

compof = {e: find(e) for e in C}
def decided(x, y): return x >= a.minObs + y and x >= a.ratio * max(y, 1)
labels = {}
for e, p in chosen.items():
    if not p: labels[e] = "none"; continue
    ties = []; here = score(chosen)
    others = [k for k in C if k != e and compof[k] == compof[e]]
    for q in [()] + C[e]:
        if q == p: continue
        t = dict(chosen); t[e] = q
        if feasible(t) and score(t) >= here: ties.append((obs(e, p), obs(e, q)))
        for k in others:
            for q2 in [()] + C[k]:
                if q2 == chosen[k]: continue
                t2 = dict(t); t2[k] = q2
                if feasible(t2) and score(t2) >= here: ties.append((obs(e, p) + obs(k, chosen[k]), obs(e, q) + obs(k, q2)))
    labels[e] = "only" if not ties else "supported" if all(decided(x, y) for x, y in ties) else "arbitrary"

n = collections.Counter(z[:-1] for p in chosen.values() for z in p)
with open(a.out, "w") as o:
    o.write("#end\tgroup\tpath\tlabel\tpairs\tbp\n")
    for e, p in sorted(chosen.items(), key=lambda kv: (ends[kv[0]], kv[0])):
        if p: o.write(f"{e}\t{ends[e]}\t{','.join(p)}\t{labels[e]}\t{obs(e, p)}\t{sum(L[z[:-1]] for z in p)}\n")
with open(a.report, "w") as o:
    o.write("#chain\tlength\tinterior\tdepth\ttarget\tused_bridges\tused_tails\tstate\n")
    for c in sorted(iv, key=lambda c: -W[c]):
        if used[c] or n[c] or any(c == z[:-1] for v in C.values() for p in v for z in p):
            st = "tail" if n[c] else "bridge" if used[c] else "candidate_unused"
            o.write(f"{c}\t{L[c]}\t{W[c]}\t{D[c]:.1f}\t{T[c]}\t{used[c]}\t{n[c]}\t{st}\n")
    for (e, c), why in sorted(stop_notes.items()):
        if not n[c]: o.write(f"# stop {e} {c} {why} depth {D.get(c, float('nan')):.1f}\n")
print("labels", dict(collections.Counter(v for v in labels.values())))
print("tails", sum(1 for p in chosen.values() if p), "with", sum(L[c] * k for c, k in n.items()) / 1e6, "Mb of chains;", len(n), "chains;", sum(L[c] for c in n) / 1e6, "Mb distinct")
