#!/usr/bin/env python3
'''
Reads that join unitig ends, and support for every graph link.

A read links end A to end B when two consecutive alignments along the read (primary or
supplementary, MAPQ >= --minMapq, each >= --minAnchor bp) run off the end of A and on to the
start of B, within --slack bp of the unitig ends. A read may also confirm a graph link through a
single alignment only where the read runs past A's end and B is aligned next, so nothing is
inferred from reads that stay inside one unitig.

Ends are written as unitig+ (leaving the unitig's 3' end) and unitig- (leaving its 5' end);
a link is the pair of ends that face each other, as "A+ B-" for A then B, both forward.

Usage: read_links.py <reads.paf.gz> <unitigs.gfa> <kept.fai> <out.tsv>
Output: end_a end_b reads graph_overlap median_query_gap
'''
import gzip, argparse, statistics, collections

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("paf"); p.add_argument("gfa"); p.add_argument("fai"); p.add_argument("out")
    p.add_argument("--minMapq", type=int, default=20)
    p.add_argument("--minAnchor", type=int, default=5000)
    p.add_argument("--slack", type=int, default=2000)
    a = p.parse_args()

    L = { l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.fai) }
    inv = { "+": "-", "-": "+" }
    ov = {}
    for l in open(a.gfa):
        if l[0] != "L": continue
        s = l.split("\t")
        if s[1] not in L or s[3] not in L: continue
        # GFA link: leave s[1] in orientation s[2], enter s[3] in orientation s[4]
        e1 = (s[1], s[2]); e2 = (s[3], inv[s[4]])
        ov[tuple(sorted((e1, e2)))] = int(s[5][:-1])

    alns = collections.defaultdict(list)
    with gzip.open(a.paf, "rt") as fh:
        for l in fh:
            s = l.split("\t", 12)
            if "tp:A:P" not in s[12]: continue
            if int(s[11]) < a.minMapq or int(s[3]) - int(s[2]) < a.minAnchor: continue
            alns[s[0]].append((int(s[2]), int(s[3]), s[4], s[5], int(s[6]), int(s[7]), int(s[8])))

    gaps = collections.defaultdict(list)
    for q, xs in alns.items():
        if len(xs) < 2: continue
        xs.sort()
        for x, y in zip(xs, xs[1:]):
            if x[3] == y[3]: continue
            # end of x that the read leaves through, and end of y it enters through
            xdone = (x[4] - x[6] <= a.slack) if x[2] == "+" else (x[5] <= a.slack)
            yready = (y[5] <= a.slack) if y[2] == "+" else (y[4] - y[6] <= a.slack)
            if not (xdone and yready): continue
            e1 = (x[3], x[2]); e2 = (y[3], inv[y[2]])
            gaps[tuple(sorted((e1, e2)))].append(y[0] - x[1])

    with open(a.out, "w") as o:
        o.write("#end_a\tend_b\treads\tgraph_overlap\tmedian_query_gap\n")
        for k in sorted(set(gaps) | set(ov)):
            g = gaps.get(k, [])
            o.write(f"{k[0][0]}{k[0][1]}\t{k[1][0]}{k[1][1]}\t{len(g)}\t{ov.get(k, -1)}\t{int(statistics.median(g)) if g else 'NA'}\n")

if __name__ == "__main__":
    main()
