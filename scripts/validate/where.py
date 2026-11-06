#!/usr/bin/env python3
'''Component and chain position at breakpoints, and the nearest join. Usage: where.py scaffolds.agp unitigs.fasta.fai evidence.tsv'''
import collections, sys
agp, fai, ev = sys.argv[1:4]
rows=collections.defaultdict(list)
for l in open(agp):
    if l[0]=="#":continue
    s=l.rstrip("\n").split("\t")
    if s[4]=="W": rows[s[0]].append((int(s[1]),int(s[2]),s[5],s[8]))
L={}
for l in open(fai):
    x=l.split("\t");L[x[0]]=int(x[1])
for l in open(ev):
    if l[0]=="#":continue
    s=l.split("\t");c,x=s[0],int(s[1])
    r=[z for z in rows[c] if z[0]<=x<=z[1]]
    # nearest junction between components
    near=min((abs(z[1]-x),z[1]) for z in rows[c])
    z=r[0] if r else None
    print(c,x,s[2],"utg",z[2] if z else "gap", "utg_len",L.get(z[2]) if z else ".", "offset_in_comp",x-z[0] if z else ".", "comp_len",z[1]-z[0]+1 if z else ".", "nearest_comp_end_kb",round(near[0]/1000,1))
