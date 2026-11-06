#!/usr/bin/env python3
'''Is a read kept by correction, overall and among reads near gaps, by quality and length? Usage: retention.py reads.tsv.gz corrected.lens.tsv.gz neargap.ids'''
import gzip,collections,sys
readsT,corrT,idsF=sys.argv[1:4]
ids=set(open(idsF).read().split())
corr={l.split('\t')[0].split(':')[0] for l in gzip.open(corrT,'rt')}
B=collections.defaultdict(lambda:[0,0,0,0])
for l in gzip.open(readsT,'rt'):
    i,n,q=l.split('\t')
    if int(n)<5000: continue
    q=min(int(float(q))//5*5,30); len_b=min(int(n)//10000*10,50)
    for key in((q,'q'),(len_b,'len')):
        b=B[key]; b[0]+=1; b[1]+=i in corr
        if i in ids: b[2]+=1; b[3]+=i in corr
print('bin kind all_kept near_gap_n near_gap_kept')
for k in sorted(B,key=lambda x:(x[1],x[0])):
    b=B[k]; print(k[0],k[1],b[0],round(b[1]/b[0],3),b[2],round(b[3]/max(b[2],1),3))
