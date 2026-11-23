#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Depth outside graph overlaps for named chains. Usage: hub_depth.py chains.agp gfa paf chain ...'''
import os, sys, collections
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from chainlib import chain_interior, end_overlaps, interior_cover
agp, gfa, paf = sys.argv[1:4]; Q = sys.argv[4:]
fl = {"+": "-", "-": "+"}; parts = collections.defaultdict(list); L = {}
for l in open(agp):
    if l[0] != "#":
        s = l.rstrip("\n").split("\t")
        if s[4] == "W": parts[s[0]].append((s[5], s[8])); L[s[0]] = int(s[2])
iv = chain_interior(parts, L, end_overlaps(gfa), Q)
cov, q = interior_cover(paf, iv, q20=True)
for c in Q:
    w = iv[c][1] - iv[c][0]
    print(c, L[c], "interior", w, "depth", round(cov[c] / w, 1) if w > 0 else "NA", "q20", round(q[c] / w, 1) if w > 0 else "NA")
