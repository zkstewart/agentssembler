#!/usr/bin/env python3
'''
Links that hic_order.py rewards: graph links with >= 5 kb overlap or >= 3 reads between chains of the same group, then
same-group bridges through short chains. Hand-made links are added afterwards, with a reason, not here.
Usage: make_order_links.py groups.tsv end_links.tsv graph.bridges.tsv out.tsv
'''
import sys
groups, ends, bridges, out = sys.argv[1:5]
g = {l.split("\t")[0]: l.rstrip().split("\t")[2] for l in open(groups) if l[0] != "#"}
seen = set()
with open(out, "w") as o:
    for l in open(ends):
        if l[0] == "#": continue
        a, b, r, q, ov = l.rstrip().split("\t")
        if (int(ov) >= 5000 or int(r) >= 3) and g[a[:-1]] == g[b[:-1]] != "unassigned" and frozenset((a, b)) not in seen:
            seen.add(frozenset((a, b))); o.write(f"{a}\t{b}\tdirect\n")
    for l in open(bridges):
        if l[0] == "#": continue
        a, b, s, *_ = l.rstrip().split("\t")
        if s == "yes" and frozenset((a, b)) not in seen:
            seen.add(frozenset((a, b))); o.write(f"{a}\t{b}\tbridge\n")
