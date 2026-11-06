#!/usr/bin/env python3
'''Strand blocks of scaffold-to-homolog alignments around candidate inversions, and reads spanning each transition.
Each PAF is named <scaffold>_vs_<homolog>.paf. regions.tsv: scaffold, start_Mb, end_Mb.
Usage: transitions.py regions.tsv reads.paf.gz a_vs_b.paf [c_vs_d.paf ...]'''
import gzip,collections,sys,os
regions, reads, *pafs = sys.argv[1:]
Q={l.split()[0]:(float(l.split()[1]),float(l.split()[2])) for l in open(regions) if l.strip() and l[0]!='#'}
Sp=collections.defaultdict(list)
for l in gzip.open(reads,"rt"):
    s=l.split("\t",12)
    if s[5] in Q and "tp:A:P" in s[12] and int(s[11])>=20: Sp[s[5]].append((int(s[7]),int(s[8])))
for f in sorted(pafs):
    q=os.path.basename(f).split("_vs_")[0]; lo,hi=Q[q]
    X=[]
    for l in open(f):
        s=l.split("\t")
        if int(s[3])-int(s[2])>=50000 and int(s[11])>=20: X.append((int(s[2]),int(s[3]),s[4],int(s[7]),int(s[8])))
    X.sort(); print("==",os.path.basename(f))
    # main strand and merged blocks around candidate
    bl=[]
    for x in X:
        if bl and bl[-1][2]==x[2]: bl[-1][1]=max(bl[-1][1],x[1]); bl[-1][3]+=x[1]-x[0]
        else: bl.append([x[0],x[1],x[2],x[1]-x[0]])
    for i,b in enumerate(bl):
        if b[1]/1e6>lo-3 and b[0]/1e6<hi+3: print(f"  {b[0]/1e6:.2f}-{b[1]/1e6:.2f} {b[2]} aligned {b[3]/1e6:.2f}")
    for u,v in zip(bl,bl[1:]):
        if abs(u[1]/1e6-lo)<1.5 or abs(u[1]/1e6-hi)<1.5:
            a,b=min(u[1],v[0]),max(u[1],v[0])
            n=sum(1 for t0,t1 in Sp[q] if t0<=a-5000 and t1>=b+5000)
            print(f"  transition {a/1e6:.3f}-{b/1e6:.3f} ({(b-a)/1000:.0f} kb) reads spanning whole interval with 5 kb anchors: {n}")
