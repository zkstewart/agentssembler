#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Hi-C log-likelihood of explicit chain orders, to compare hypotheses directly (same decay model as hic_order.py).
Each hypothesis: name, then oriented chains, e.g. "A chain1+ chain2- chain3+".
Pairs involving chains listed in the optional fourth argument (chain_stats.tsv, within < 700/Mb) are ignored, as in hic_order.py.
Usage: score_orders.py pairs.tsv.gz chains.fai hypotheses.txt [chain_stats.tsv]
'''
import sys, gzip, collections, importlib.util, os
import numpy as np
sp = importlib.util.spec_from_file_location("ho", os.path.join(os.path.dirname(os.path.abspath(__file__)), "hic_order.py"))
ho = importlib.util.module_from_spec(sp); sp.loader.exec_module(ho)
pairs, fai, hyp = sys.argv[1:4]; GAP = 20000
low = {l.split("\t")[0] for l in open(sys.argv[4]) if l[0] != "#" and float(l.split("\t")[2]) < 700} if len(sys.argv) > 4 else set()
L = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(fai)}
H = [l.split() for l in open(hyp) if l.strip() and l[0] != "#"]
want = {x[:-1] for h in H for x in h[1:]}
within = []; P = collections.defaultdict(list)
with gzip.open(pairs, "rt") as fh:
    for l in fh:
        x, u, y, v = l.split("\t"); u = int(u); v = int(v)
        if x == y: within.append(abs(u - v))
        elif x in want and y in want:
            if x < y: P[(x, y)].append((u, v))
            else: P[(y, x)].append((v, u))
logp = ho.decay_model(within)
A = {k: np.asarray(v, dtype=float) for k, v in P.items()}
def score(order):
    st = []; pos = 0
    for c, o in order: st.append(pos); pos += L[c] + GAP
    s = 0.0; per = {}
    for m in range(len(order)):
        for n in range(m + 1, len(order)):
            i, oi = order[m]; j, oj = order[n]
            k = (i, j) if i < j else (j, i)
            if k not in A or i in low or j in low: continue
            arr = A[k]; xi, xj = (arr[:, 0], arr[:, 1]) if i < j else (arr[:, 1], arr[:, 0])
            di = L[i] - xi if oi == "+" else xi
            dj = xj if oj == "+" else L[j] - xj
            v = float(logp(di + dj + st[n] - st[m] - L[i]).sum()); s += v; per[(i, j)] = v
    return s, per
PAIRS_INFO = None
res = []
for h in H:
    s, per = score([(x[:-1], x[-1]) for x in h[1:]]); res.append((h[0], s, per))
best = max(r[1] for r in res)
for n, s, per in res:
    print(f"{n}\t{s:.1f}\tdelta_from_best {s-best:.1f}\t" + " ".join(f"{i[-4:]}-{j[-4:]}:{v:.0f}" for (i, j), v in sorted(per.items())))
