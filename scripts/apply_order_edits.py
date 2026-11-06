#!/usr/bin/env python3
'''Apply a hand-edit table to a chain order and groups file. Actions: unplace (out of the order), append or prepend (added to a
group's end and assigned to it). Every row must carry a reason. Ranks are renumbered; added rows have no Hi-C statistics.
Usage: apply_order_edits.py edits.tsv order.tsv groups.tsv order_out.tsv groups_out.tsv'''
import sys, collections
ed, od, gr, oo, go = sys.argv[1:6]
drop = set(); add = collections.defaultdict(list); to = {}; L = {}
hdr = None; G = []
for l in open(gr):
    if l[0] == "#": hdr = l; continue
    s = l.rstrip("\n").split("\t"); L[s[0]] = s[1]; G.append(s)
for l in open(ed):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if not s[5].strip(): raise SystemExit(f"no reason for {s[0]}")
    if s[0] not in L: raise SystemExit(f"unknown chain {s[0]}")
    if s[1] == "unplace": drop.add(s[0])
    elif s[1] in ("append", "prepend"): add[s[2]].append((s[1], s[0], s[4])); to[s[0]] = s[2]
    else: raise SystemExit(f"unknown action {s[1]}")
R = collections.defaultdict(list); head = ""
for l in open(od):
    if l[0] == "#": head = l; continue
    s = l.rstrip("\n").split("\t")
    if s[2] not in drop: R[s[0]].append(s)
for g, v in add.items():
    for act, c, o in v:
        row = [g, "0", c, o, L[c], "0", "0", "0", "edit"]
        if act == "append": R[g].append(row)
        else: R[g].insert(0, row)
with open(oo, "w") as o:
    o.write(head)
    for g in sorted(R):
        for i, s in enumerate(R[g], 1):
            s[1] = str(i); o.write("\t".join(s) + "\n")
with open(go, "w") as o:
    if hdr: o.write(hdr)
    for s in G:
        if s[0] in to: s[2] = to[s[0]]
        o.write("\t".join(s) + "\n")
print(len(drop), "unplaced,", len(to), "added")
