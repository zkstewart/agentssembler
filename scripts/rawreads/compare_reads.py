#!/usr/bin/env python3
"""Raw vs corrected reads: what correction removed overall, and at each gap of the pre-fill AGPs, with control sites.
reads.tsv.gz: id, length, mean Q. corrected.lens.tsv.gz: id (optional :N for a split fragment), length.
Each SET is agp,fai,raw.paf.gz,corrected.paf.gz, in the coordinates of that AGP; give one per mapping target.
usage: compare_reads.py outprefix reads.tsv.gz corrected.lens.tsv.gz SET [SET ...]
Corrected spans are keyed by (id, read length) so the two halves of a split read never count as one read."""
import gzip, sys, bisect, collections, statistics as st
out, readsTab, corrTab, *sets = sys.argv[1:]
SLACK, ANCH, WIN = 1000, 5000, 60000
DIST = ((0,1000),(1000,2000),(2000,5000),(5000,10000),(10000,20000),(20000,50000))
corr = collections.defaultdict(list)
for l in gzip.open(corrTab, 'rt'):
    i, n = l.split('\t'); corr[i.split(':')[0]].append(int(n))
print('corrected ids', len(corr), 'records', sum(map(len, corr.values())), 'bases', sum(map(sum, corr.values())), 'split ids', sum(len(v) > 1 for v in corr.values()))
bins = [0, 1000, 5000, 10000, 20000, 50000, 100000, 10**9]
S = [[0]*7 for _ in bins]; seen = 0
Q = collections.defaultdict(lambda: [0,0,0,0]); qual = {}
with gzip.open(readsTab, 'rt') as f:
    for l in f:
        i, n, q = l.rstrip('\n').split('\t'); n = int(n)
        b = bisect.bisect_right(bins, n)-1
        r = S[b]; r[0] += 1; r[1] += n
        c = corr.get(i)
        if n >= 5000:
            qual[i] = float(q); qb = Q[min(int(float(q))//5*5, 30)]; qb[0] += 1; qb[1] += n; qb[2] += bool(c); qb[3] += sum(c) if c else 0
        if c:
            seen += 1; r[2] += 1; r[3] += sum(c); r[4] += len(c) > 1
            if len(c) == 1 and c[0] < 0.9*n: r[5] += 1
            if len(c) == 1 and c[0] < 0.5*n: r[6] += 1
print('corrected ids found in raw', seen, 'of', len(corr))
with open(out+'.by_length.tsv', 'w') as w:
    w.write('raw_min\traw_max\treads\tbases\twith_corrected\tfrac_reads\tcorrected_bases\tcorr_over_raw_bases\tsplit\tsingle_lt90pct\tsingle_lt50pct\n')
    for k in range(len(bins)-1):
        r = S[k]
        w.write(f'{bins[k]}\t{bins[k+1]}\t{r[0]}\t{r[1]}\t{r[2]}\t{r[2]/max(r[0],1):.3f}\t{r[3]}\t{r[3]/max(r[1],1):.3f}\t{r[4]}\t{r[5]}\t{r[6]}\n')

with open(out+'.by_q.tsv', 'w') as w:
    w.write('q_bin\treads\tbases\tfrac_reads_kept\tcorr_over_raw_bases\n')
    for k in sorted(Q):
        r = Q[k]; w.write(f'{k}\t{r[0]}\t{r[1]}\t{r[2]/r[0]:.3f}\t{r[3]/r[1]:.3f}\n')
SP = open(out+'.spanners.tsv', 'w'); SP.write('t\tgs\tread\tlen\tq\tmapq_min\tin_corrected\tcorrected_spans\n')

def gaps(agp):
    pos = collections.defaultdict(list)
    for l in open(agp):
        if l[0] == '#': continue
        s = l.rstrip('\n').split('\t')
        if s[4] in 'NU': pos[s[0]].append((int(s[1])-1, int(s[2])))
    return pos
def load(paf, sites):
    """for each site, the alignments overlapping it +- WIN, keyed by read"""
    idx = {t: sorted((a, b, n) for n, (a, b) in enumerate(v)) for t, v in sites.items()}
    R = {(t, a, b): collections.defaultdict(list) for t, v in sites.items() for a, b in v}
    for l in gzip.open(paf, 'rt'):
        s = l.split('\t', 12)
        t = s[5]
        if t not in idx: continue
        ts, te = int(s[7]), int(s[8])
        for a, b, _ in idx[t]:
            if ts < b+WIN and te > a-WIN:
                ql, qs, qe = int(s[1]), int(s[2]), int(s[3])
                oqs, oqe = (qs, qe) if s[4] == '+' else (ql-qe, ql-qs)
                R[(t, a, b)][(s[0], ql)].append((ts, te, s[4], int(s[11]), oqs, oqe))
    return R
rows = []; cov = collections.Counter(); 
for agp, fai, rawpaf, corpaf in (x.split(',') for x in sets):
    h = agp
    G = gaps(agp)
    lens = {l.split('\t')[0]: int(l.split('\t')[1]) for l in open(fai)}
    sites = {}
    for t in G:
        sites[t] = [(a, b, 'gap') for a, b in G[t]]
        sites[t] += [(x, x+100, 'ctl') for x in range(200000, lens[t]-200000, 500000) if all(abs(x-a) > 120000 for a, b in G[t])]
    plain = {t: [(a, b) for a, b, _ in v] for t, v in sites.items()}
    raw = load(rawpaf, plain); cor = load(corpaf, plain)
    print(h, 'loaded', flush=True)
    for t, v in sites.items():
        for gs, ge, kind in v:
            row = dict(kind=kind, h=h, t=t, gs=gs, ge=ge); edge = {}
            for src, R in (('raw', raw[(t, gs, ge)]), ('corr', cor[(t, gs, ge)])):
                L, Rr, thru = {}, {}, {}
                for key, al in R.items():
                    for a in al:
                        if a[0] <= gs-ANCH and a[1] >= ge+ANCH: thru[key] = a
                        if gs-SLACK <= a[1] <= gs+SLACK and a[1]-a[0] >= ANCH: L[key] = a
                        if ge-SLACK <= a[0] <= ge+SLACK and a[1]-a[0] >= ANCH: Rr[key] = a
                        for sd, x0, sg in (('L', gs, -1), ('R', ge, 1)):
                            for d0, d1 in DIST:
                                lo, hi = (x0-d1, x0-d0) if sg < 0 else (x0+d0, x0+d1)
                                ov = min(a[1], hi)-max(a[0], lo)
                                if ov > 0: cov[(kind, src, d0)] += ov
                span, spanq = set(thru), {k for k in thru if thru[k][3] >= 20}
                for k in set(L) & set(Rr):
                    a, b = L[k], Rr[k]
                    if a[2] == b[2]:
                        span.add(k)
                        if a[3] >= 20 and b[3] >= 20: spanq.add(k)
                row[src+'_span'], row[src+'_spanq20'], row[src+'_L'], row[src+'_R'] = len(span), len(spanq), len(L), len(Rr)
                edge[src] = (L, Rr, span)
            if kind == 'gap':
                cs = {k[0].split(':')[0] for k in edge['corr'][2]}
                for k in edge['raw'][2]:
                    r = raw[(t, gs, ge)][k]; m = min(a[3] for a in r)
                    SP.write(f"{t}\t{gs}\t{k[0]}\t{k[1]}\t{qual.get(k[0], -1):.1f}\t{m}\t{int(k[0] in corr)}\t{int(k[0] in cs)}\n")
            cL = collections.defaultdict(list)
            for sd, D in zip('LR', edge['corr'][:2]):
                for (i, ql), a in D.items(): cL[i.split(':')[0]].append((sd, ql, a))
            over = []; nabs = nkeep = nshort = 0
            for sd, D in zip('LR', edge['raw'][:2]):
                for (i, ql), a in D.items():
                    if (i, ql) in edge['raw'][2]: continue
                    ov = ql-a[5] if sd == 'L' else a[4]
                    if ov < 2000: continue
                    over.append(ov)
                    c = [x for x in cL.get(i, []) if x[0] == sd]
                    if i not in corr: nabs += 1
                    elif not c: nshort += 1
                    else:
                        ovc = c[0][1]-c[0][2][5] if sd == 'L' else c[0][2][4]
                        if ovc >= 0.7*ov: nkeep += 1
                        else: nshort += 1
            row.update(raw_overhang_reads=len(over), overhang_med=int(st.median(over)) if over else 0, oh_absent_corr=nabs, oh_kept_corr=nkeep, oh_shortened_corr=nshort)
            ids = collections.defaultdict(set)
            for sd, D in zip('LR', edge['corr'][:2]):
                for (i, ql) in D: ids[i.split(':')[0]].add((sd, ql))
            row['corr_split_across'] = sum(len({x[0] for x in w}) == 2 and len({x[1] for x in w}) == 2 for w in ids.values())
            rows.append(row)
cols = 'kind h t gs ge raw_span raw_spanq20 corr_span corr_spanq20 raw_L raw_R corr_L corr_R raw_overhang_reads overhang_med oh_absent_corr oh_kept_corr oh_shortened_corr corr_split_across'.split()
with open(out+'.sites.tsv', 'w') as w:
    w.write('\t'.join(cols)+'\n')
    for r in rows: w.write('\t'.join(str(r[c]) for c in cols)+'\n')
with open(out+'.taper.tsv', 'w') as w:
    w.write('kind\tdist\traw_bases\tcorr_bases\tcorr_over_raw\n')
    for kind in ('gap', 'ctl'):
        for d, _ in DIST:
            a, b = cov[(kind,'corr',d)], cov[(kind,'raw',d)]
            w.write(f'{kind}\t{d}\t{b}\t{a}\t{a/max(b,1):.3f}\n')
