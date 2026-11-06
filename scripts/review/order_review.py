#!/usr/bin/env python3
'''Per adjacent pair in an order: direct graph/read link? same-group bridge candidates? plus chains' telomere ends.
Usage: order_review.py order.tsv end_links.tsv graph.bridges.tsv telomeres.tsv'''
import sys, collections
od, el, br, te = sys.argv[1:5]
fl = {"+": "-", "-": "+"}
direct = {}
hdr = None
for l in open(el):
    s = l.rstrip("\n").split("\t")
    if l[0] == "#": hdr = s; continue
    direct[(s[0], s[1])] = s[2:]; direct[(s[1], s[0])] = s[2:]
print("end_links header", hdr)
bridge = collections.defaultdict(list)
for l in open(br):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    bridge[(s[0], s[1])].append(s[3]); bridge[(s[1], s[0])].append(s[3])
tel = collections.defaultdict(list)
for l in open(te):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t"); tel[s[0]] = s[1:]
print("telomere ex", list(tel.items())[:1])
O = collections.defaultdict(list)
for l in open(od):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t"); O[s[0]].append(s)
st = collections.Counter()
for g, rows in O.items():
    print(f"\n### {g}")
    for i, r in enumerate(rows):
        c, o, n = r[2], r[3], int(r[4])
        link = ""
        if i + 1 < len(rows):
            d = rows[i + 1]
            # right end of this chain = chain+ when orientation + ; left end of next = chain- when orientation +
            a = c + o; b = d[2] + fl[d[3]]
            if (a, b) in direct: link = "direct"
            elif (a, b) in bridge: link = f"bridge({len(bridge[(a,b)])})"
            else: link = "NONE"
            st[link.split("(")[0]] += 1
        print(f"{r[1]:>3} {c:<14}{o} {n/1e6:7.3f} links={r[5]:>6} oLLR={float(r[6]):>9.1f} pLLR={float(r[7]):>9.1f} {link:<10}")
print(st)
