#!/usr/bin/env python3
'''
Links between chain ends, from the two kinds of evidence that outrank Hi-C.
  graph: a GFA link between the terminal unitigs of two chains, leaving each chain outwards. Chaining took
         only reciprocal 1:1 links, so these are the links it left unused at forks.
  reads: a read whose consecutive alignments (MAPQ >= --minMapq, >= --minAnchor bp, on different chains) leave
         one chain end and enter another, each within --slack bp of the end.
Ends are chain+ (3' end) and chain- (5' end). Reads that stop short of a chain end say nothing.

Usage: chain_end_links.py chains.agp unitigs.gfa ont_vs_chains.paf.gz chains.fai out.tsv
Output: end_a end_b reads q_gap_median graph_overlap
'''
import argparse, collections, gzip, statistics
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("agp"); p.add_argument("gfa"); p.add_argument("paf"); p.add_argument("fai"); p.add_argument("out")
p.add_argument("--minMapq", type=int, default=20)
p.add_argument("--minAnchor", type=int, default=5000)
p.add_argument("--slack", type=int, default=5000)
a = p.parse_args()
inv = {"+": "-", "-": "+"}
L = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.fai)}

# Terminal unitigs: leaving the chain through its 3' end means leaving the last unitig in its AGP strand;
# leaving through its 5' end means leaving the first unitig against its AGP strand.
rows = collections.defaultdict(list)
for l in open(a.agp):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if s[4] == "W": rows[s[0]].append((s[5], s[8]))
exit_of = {}    # (unitig, direction leaving the unitig) -> chain end
for c, r in rows.items():
    first_u, first_o = r[0]; last_u, last_o = r[-1]
    exit_of[(last_u, last_o)] = (c, "+")
    exit_of[(first_u, inv[first_o])] = (c, "-")

graph = {}
for l in open(a.gfa):
    if l[0] != "L": continue
    s = l.split("\t")
    e1 = exit_of.get((s[1], s[2])); e2 = exit_of.get((s[3], inv[s[4]]))
    if e1 and e2 and e1[0] != e2[0]:
        ov = int(s[5].strip()[:-1])
        k = tuple(sorted((e1, e2)))
        graph[k] = max(graph.get(k, 0), ov)

alns = collections.defaultdict(list)
with gzip.open(a.paf, "rt") as fh:
    for l in fh:
        s = l.split("\t", 12)
        if "tp:A:P" not in s[12] or int(s[11]) < a.minMapq or int(s[3]) - int(s[2]) < a.minAnchor: continue
        alns[s[0]].append((int(s[2]), int(s[3]), s[4], s[5], int(s[6]), int(s[7]), int(s[8])))
gaps = collections.defaultdict(list)
for q, xs in alns.items():
    if len(xs) < 2: continue
    xs.sort()
    for x, y in zip(xs, xs[1:]):
        if x[3] == y[3]: continue
        xdone = (x[4] - x[6] <= a.slack) if x[2] == "+" else (x[5] <= a.slack)
        yready = (y[5] <= a.slack) if y[2] == "+" else (y[4] - y[6] <= a.slack)
        if xdone and yready:
            gaps[tuple(sorted(((x[3], x[2]), (y[3], inv[y[2]]))))].append(y[0] - x[1])

with open(a.out, "w") as o:
    o.write("#end_a\tend_b\treads\tq_gap_median\tgraph_overlap\n")
    for k in sorted(set(gaps) | set(graph)):
        g = gaps.get(k, [])
        o.write(f"{k[0][0]}{k[0][1]}\t{k[1][0]}{k[1][1]}\t{len(g)}\t{int(statistics.median(g)) if g else 'NA'}\t{graph.get(k, -1)}\n")
print(f"{len(graph)} graph links, {len(gaps)} read-supported end pairs, {sum(1 for k in gaps if len(gaps[k]) >= 3)} with >= 3 reads")
