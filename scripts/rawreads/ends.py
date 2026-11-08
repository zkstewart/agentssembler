#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Reads that reach a sequence end, and where the reads continue past it. Each END is target:START or target:END; each SET is label=reads.paf.gz.
Usage: ends.py target.fasta.fai SET[,SET...] END [END ...]'''
import gzip,collections,sys
fai,sets,*ends=sys.argv[1:]
L={l.split('\t')[0]:int(l.split('\t')[1]) for l in open(fai)}
E=[(x.rsplit(':',1)[0],x.rsplit(':',1)[1]) for x in ends]
want={e[0] for e in E}
R=collections.defaultdict(list)
for src,p in (x.split('=',1) for x in sets.split(',')):
    A=collections.defaultdict(list)
    for l in gzip.open(p,'rt'):
        s=l.split('\t',12); A[s[0]].append((int(s[2]),int(s[3]),s[4],s[5],int(s[7]),int(s[8]),int(s[11]),int(s[1])))
    print('==',src)
    for t,e in E:
        n=0;co=collections.Counter();reach=0;ov=[]
        for q,al in A.items():
            for a in al:
                if a[3]!=t or a[5]-a[4]<3000: continue
                near=(a[5]>=L[t]-2000) if e=='END' else (a[4]<=2000)
                if not near: continue
                reach+=1
                # bases of the read beyond the end
                beyond=(a[7]-a[1] if a[2]=='+' else a[0]) if e=='END' else (a[0] if a[2]=='+' else a[7]-a[1])
                if beyond<3000: continue
                n+=1; ov.append(beyond)
                lo,hi=(a[1],a[7]) if (a[2]=='+')==(e=='END') else (0,a[0])
                o=[b for b in al if b is not a and min(hi,b[1])-max(lo,b[0])>=3000]
                if o:
                    b=max(o,key=lambda x:x[1]-x[0]); co[f'{b[3]}:{b[4]//100000}00k:q{"hi" if b[6]>=20 else "lo"}']+=1
                else: co['none']+=1
        print(t,e,'reach',reach,'beyond3kb',n,'med',sorted(ov)[len(ov)//2] if ov else 0,co.most_common(5))
    del A
