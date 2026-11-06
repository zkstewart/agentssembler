#!/usr/bin/env python3
'''Merged colinear blocks of a scaffold-to-homolog PAF, >= 100 kb aligned, for a query interval. Usage: blocks.py paf qstart_Mb qend_Mb'''
import sys
paf, a, b = sys.argv[1], float(sys.argv[2]) * 1e6, float(sys.argv[3]) * 1e6
R = []
for l in open(paf):
    s = l.split("\t")
    if int(s[10]) < 20000: continue
    R.append((int(s[2]), int(s[3]), s[4], int(s[7]), int(s[8]), int(s[10]), int(s[9])))
R.sort()
out = []
for r in R:
    if r[1] < a or r[0] > b: continue
    if out and out[-1][2] == r[2] and r[0] - out[-1][1] < 500000 and ((r[2] == "+" and -500000 < r[3] - out[-1][4] < 1000000) or (r[2] == "-" and -500000 < out[-1][3] - r[4] < 1000000)):
        o = out[-1]; o[1] = max(o[1], r[1]); o[3] = min(o[3], r[3]); o[4] = max(o[4], r[4]); o[5] += r[5]; o[6] += r[6]
    else: out.append(list(r))
for o in out:
    if o[5] >= 50000: print(f"  h2 {o[0]/1e6:6.2f}-{o[1]/1e6:6.2f} {o[2]} h1 {o[3]/1e6:6.2f}-{o[4]/1e6:6.2f} aln {o[5]/1e3:.0f} kb id {o[6]/o[5]:.3f}")
