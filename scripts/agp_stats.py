#!/usr/bin/env python3
'''Size, gaps, contig N50 and scaffold N50 of one or more AGPs. Usage: agp_stats.py a.agp [b.agp ...]'''
import sys
for f in sys.argv[1:]:
    gaps=gb=0; ctg=[]; sc={}; cur=0; prev=None; 
    for l in open(f):
        if l[0]=="#": continue
        s=l.rstrip("\n").split("\t")
        if s[0]!=prev:
            if cur: ctg.append(cur)
            cur=0; prev=s[0]
        sc[s[0]]=int(s[2])
        if s[4] in "NU":
            gaps+=1; gb+=int(s[5])
            if cur: ctg.append(cur)
            cur=0
        else: cur+=int(s[2])-int(s[1])+1
    if cur: ctg.append(cur)
    def n50(v):
        v=sorted(v,reverse=True);t=sum(v);c=0
        for x in v:
            c+=x
            if c*2>=t: return x
    print(f, f"scaffolds={len(sc)} size={sum(sc.values())/1e6:.1f}Mb gaps={gaps} gap_bp={gb} contigs={len(ctg)} ctgN50={n50(ctg)/1e6:.2f} scafN50={n50(list(sc.values()))/1e6:.2f}")
