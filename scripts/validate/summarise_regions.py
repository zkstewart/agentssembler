#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Summarise region-to-assembly alignments: top targets of each region by strand. Usage: summarise_regions.py regions/*.paf.gz'''
import gzip,collections,sys
print("#region\tlen_Mb\ttarget\tprimary_Mb_plus\tprimary_Mb_minus\tsecondary_Mb_plus\tsecondary_Mb_minus\tmean_id")
for f in sorted(sys.argv[1:]):
    C=collections.defaultdict(lambda:[0,0,0,0,0,0]); q=None;L=0
    for l in gzip.open(f,"rt"):
        s=l.split("\t");q=s[0];L=int(s[1]);al=int(s[3])-int(s[2])
        if al<10000: continue
        pr="tp:A:P" in l
        k=(0 if pr else 2)+(0 if s[4]=="+" else 1)
        C[s[5]][k]+=al; C[s[5]][4]+=int(s[9]); C[s[5]][5]+=int(s[10])
    tops=sorted(C.items(),key=lambda t:-(t[1][0]+t[1][1]+t[1][2]+t[1][3]))[:5]
    for t,v in tops: print(q,round(L/1e6,2),t,*[round(x/1e6,2) for x in v[:4]],round(v[4]/max(v[5],1),3),sep="\t")
