#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Inverted homolog blocks merged along h2: alignments >= 50 kb on the opposite strand, chain within 300 kb, total aligned >= 400 kb. Usage: merge.py chr*.h2_vs_h1.paf'''
import sys
for f in sys.argv[1:]:
    R=[]
    for l in open(f):
        s=l.split("\t")
        if s[4]=="-" and int(s[10])>=50000: R.append((int(s[2]),int(s[3]),int(s[7]),int(s[8]),int(s[10]),s[0],s[5]))
    R.sort(); C=[]
    for r in R:
        if C and r[0]-C[-1][1]<300000 and abs(C[-1][2]-r[3])<500000: c=C[-1]; c[1]=max(c[1],r[1]); c[2]=min(c[2],r[2]); c[3]=max(c[3],r[3]); c[4]+=r[4]
        else: C.append(list(r))
    for c in C:
        if c[4]>=400000: print(f"{c[5]}:{c[0]}:{c[1]}\t{c[6]}:{c[2]}:{c[3]}\taligned {c[4]/1e3:.0f} kb")
