#!/usr/bin/env python3
'''For scaffold ends without a telomere: terminal chain, tips from graph, and evidence linking unplaced telomeric chains
(graph links, reads split between chains, Hi-C to terminal chains). Reads come from the chain-level PAF.
Usage: tel_ends.py WC.graph.agp chains.agp order.edited.tsv names.tsv telomeres.tsv end_links.tsv groups.tsv ont_vs_chains.paf.gz pairs.tsv.gz'''
import sys, gzip, collections
agp, cagp, od, nm, te, el, gr, paf, pairs = sys.argv[1:10]
names = {}; flip = {}
for l in open(nm):
    if l[0] != "#":
        g, n, f = l.rstrip().split("\t"); names[g] = n; flip[g] = f == "flip"
O = collections.defaultdict(list)
for l in open(od):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); O[s[0]].append((s[2], s[3]))
tel = {}; cnt = {}
for l in open(te):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); tel[s[0]] = s[4]; cnt[s[0]] = (s[2], s[3])
g = {}; L = {}
for l in open(gr):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t"); L[s[0]] = int(s[1]); g[s[0]] = s[2]
adj = collections.defaultdict(list)
for l in open(el):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t"); adj[s[0]].append((s[1], s[2], s[4])); adj[s[1]].append((s[0], s[2], s[4]))
fl = {"+": "-", "-": "+"}
ends = []
for gg, seq in O.items():
    first = seq[0]; last = seq[-1]
    # outward end of the first and last chain in order, before the naming flip
    ends.append((gg, "first", first[0], fl[first[1]], first))
    ends.append((gg, "last", last[0], last[1], last))
print("chain-level telomeres at order ends (a chain-end facing out carrying a telomere is 'ok')")
want = {}
for gg, side, c, e, rec in ends:
    t = tel.get(c, "--"); i = 0 if e == "-" else 1
    out_t = t[i] == "T"; in_t = t[1 - i] == "T"
    if ".g" in c and c not in tel:
        b = c.split(".g")[0]; t = tel.get(b, "--"); out_t = t[i] == "T"; in_t = t[1 - i] == "T"
    if not out_t:
        print(f"{names[gg]:8} {gg} order-{side}: {c}{e} {L[c]/1e3:.0f}kb tel={t} {'TELOMERE_FACES_IN' if in_t else ''} links: {adj[c+e][:6]}")
        want[c + e] = names[gg]
un = [c for c in tel if tel[c] != "--" and g.get(c) == "unassigned"]
print("\nunplaced telomeric chains:", [(c, L[c], tel[c], cnt[c]) for c in un])
for c in un:
    i = "-" if tel[c][0] == "T" else "+"                    # telomere end; inner end is the other
    inner = c + fl[i]
    print(f" {c}: inner end {inner} links {adj[inner][:8]}")
Hi = {c: collections.Counter() for c in un}
for l in gzip.open(pairs, "rt"):
    s = l.split("\t"); a, b = s[0], s[2]
    if a in Hi and a != b: Hi[a][b] += 1
    if b in Hi and a != b: Hi[b][a] += 1
gsz = collections.Counter()
for c, x in g.items(): gsz[x] += L[c]
for c, cn in Hi.items():
    pg = collections.Counter()
    for t, v in cn.items(): pg[g.get(t, "?")] += v
    print(f" {c} Hi-C total {sum(cn.values())} by group per Mb: " + " ".join(f"{k}:{v/(gsz[k]/1e6):.1f}" for k, v in pg.most_common(5)))
    print("   top chains: " + " ".join(f"{k}[{g.get(k)}]:{v}" for k, v in cn.most_common(6)))
# reads: alignments on an unplaced telomeric chain and another chain, near both ends
A = collections.defaultdict(list)
for l in gzip.open(paf, "rt"):
    s = l.split("\t", 12)
    if int(s[11]) < 1: continue
    A[s[0]].append((s[5], int(s[6]), int(s[7]), int(s[8]), s[4], int(s[11])))
R = collections.defaultdict(collections.Counter)
for q, xs in A.items():
    if len(xs) < 2: continue
    for x in xs:
        if x[0] in Hi:
            for y in xs:
                if y[0] != x[0]: R[x[0]][y[0]] += 1
for c in un: print(f" {c} reads shared with: {R[c].most_common(6)}")
