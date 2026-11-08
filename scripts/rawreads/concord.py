#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Does correction move reads? Best alignment of each read in the raw and corrected mapping to one target, for reads at q20 in both.
Usage: concord.py raw.paf.gz corrected.paf.gz hap_switch.tsv'''
import gzip,collections,sys
rawP,corP,outF=sys.argv[1:4]
def best(p):
    B={}
    for l in gzip.open(p,'rt'):
        s=l.split('\t',12); i=s[0].split(':')[0]; bl=int(s[10])
        if i not in B or bl>B[i][0]: B[i]=(bl,s[5],int(s[7]),int(s[8]),s[4],int(s[11]),int(s[9])/bl)
    return B
R=best(rawP); C=best(corP)
print('raw ids',len(R),'corr ids',len(C))
n=same=hap=other=far=0; hc=collections.Counter(); idr=[];idc=[]
ex=[]
for i,c in C.items():
    r=R.get(i)
    if not r: hc['corr_not_in_raw_map']+=1; continue
    if r[5]<20 or c[5]<20: hc['lowq_either']+=1; continue
    n+=1; idr.append(r[6]); idc.append(c[6])
    ov=min(r[3],c[3])-max(r[2],c[2])
    if r[1]==c[1] and ov>0: same+=1; continue
    base=lambda t:t.rsplit('_',1)[0]
    if r[1]!=c[1] and base(r[1])==base(c[1]) and not r[1].startswith('unplaced'): hap+=1; ex.append((i,r[1],r[2],c[1],c[2]))
    elif r[1]==c[1]: far+=1
    else: other+=1
print('both q20',n,'same place',same,'other hap same chrom',hap,'same scaffold elsewhere',far,'other sequence',other,hc)
import statistics as st
print('median approx identity to target raw',round(st.median(idr),4),'corrected',round(st.median(idc),4))
open(outF,'w').write('\n'.join('\t'.join(map(str,e)) for e in ex)+'\n')
# phase balance: corrected-only vs dropped, by hap
