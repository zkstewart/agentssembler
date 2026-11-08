#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Order and orient the chains of each haplotype group from evidence about the sample only.

Likelihood: every Hi-C pair between two chains implies a genomic distance given their order and
orientation; log P(distance) comes from the decay measured on pairs within chains, mixed with a
small flat background so that stray pairs (e.g. homolog crosstalk) do not dominate.

Prior: a link between two chain ends (--endLinks: graph overlaps, or reads joining the ends) adds
--linkBonus to the likelihood when those ends face each other in the order, because graph and read
evidence outranks Hi-C.

Telomeres (--telomeres, telomeres.py output on the chains) add --telBonus for each scaffold end where a chain's
telomeric end faces outwards: a telomere can only be at a chromosome end. Chains with a telomere or an end link are
always ordered, however few pairs they have.

Low-information chains (--stats, within-chain density below --minWithin; mostly satellite arrays) are handled in a
second stage: their few uniquely placed pairs do not follow linear distance, yet a few hundred of them can outweigh the
order of megabases of informative chains. Stage 1 orders with their pairs ignored (links and telomeres still count).
Stage 2 counts their pairs but moves only blocks made entirely of low-information chains, so the backbone is fixed.

Synteny (--synteny, synteny_links.py) adds --synBonus per homolog that supports a junction: chains next to each other on a
haplotype usually have alleles next to each other on its homologs. It is deliberately small, to settle what Hi-C leaves
near a tie without overriding strong Hi-C, which is what a real single-haplotype SV looks like.

Search: greedy insertion of chains longest first (every position and orientation), then local moves
(move any contiguous block to any position in either orientation, which includes reinserting one chain,
flipping and reversing a segment) to convergence. This is repeated from --starts random orders, because a
single block move can be blocked when two blocks must move together (e.g. both carry a telomere bonus);
the best result over all starts is kept. Chains shorter than --minLen or with
fewer than --minPairs links to the group are appended unordered as 'unplaced_in_group'.

Output columns: group, rank, chain, orientation, length, links to group, orientation_llr, position_llr,
end_link (whether the junction with the previous chain has a link). The LLRs are the best alternative's
loss in log-likelihood: small values mean Hi-C does not settle that chain's orientation or position.

Usage: hic_order.py groups.tsv pairs.tsv.gz chains.fai out.tsv [--endLinks links.tsv]
'''
import gzip, math, argparse, collections, sys
import numpy as np

def decay_model(pairs_within, bg=1e-3):
    d = np.asarray(pairs_within, dtype=float); d = d[d > 0]
    edges = np.logspace(2, 9, 141)
    h, _ = np.histogram(d, bins=edges)
    w = np.diff(edges); dens = (h + 0.5) / w / max(h.sum(), 1)
    c = np.sqrt(edges[:-1] * edges[1:]); lc, ld = np.log(c), np.log(dens)
    # beyond the observed range, extend with the slope of the last well-populated decade
    ok = np.where(h >= 50)[0]; last = ok[-1]; first = ok[max(0, len(ok) - 20)]
    slope = (ld[last] - ld[first]) / (lc[last] - lc[first] + 1e-9)
    ld2 = ld.copy(); ld2[last + 1:] = ld[last] + slope * (lc[last + 1:] - lc[last])
    floor = np.log(bg) + ld[ok].min()
    def logp(x):
        y = np.interp(np.log(np.maximum(x, 100.0)), lc, ld2)
        return np.logaddexp(y, floor)
    return logp

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("groups"); p.add_argument("pairs"); p.add_argument("fai"); p.add_argument("out")
    p.add_argument("--endLinks", help="end_a end_b weight, ends as chain+ / chain-")
    p.add_argument("--linkBonus", type=float, default=200.0)
    p.add_argument("--starts", type=int, default=40)
    p.add_argument("--synteny", help="synteny_links.py output")
    p.add_argument("--synBonus", type=float, default=1000.0)
    p.add_argument("--stats", help="chain_stats.tsv")
    p.add_argument("--minWithin", type=float, default=700)
    p.add_argument("--telomeres", help="telomeres.py output for the chains")
    p.add_argument("--telBonus", type=float, default=20000.0)
    p.add_argument("--gap", type=int, default=20000, help="assumed gap between adjacent chains, bp")
    p.add_argument("--minLen", type=int, default=20000)
    p.add_argument("--minPairs", type=int, default=10)
    p.add_argument("--maxPairs", type=int, default=3000, help="subsample per chain pair (scaled back in the score)")
    a = p.parse_args()
    L = { l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.fai) }
    gr = collections.defaultdict(list); of = {}
    for l in open(a.groups):
        if l[0] == "#": continue
        c, n, g = l.rstrip().split("\t")
        if g != "unassigned": gr[g].append(c); of[c] = g
    within = []; P = collections.defaultdict(list)
    with gzip.open(a.pairs, "rt") as fh:
        for l in fh:
            x, u, y, v = l.split("\t"); u = int(u); v = int(v)
            if x == y: within.append(abs(u - v)); continue
            if x in of and y in of and of[x] == of[y]:
                if x < y: P[(x, y)].append((u, v))
                else: P[(y, x)].append((v, u))
    logp = decay_model(within)
    rng = np.random.default_rng(1)
    A = {}; W = {}
    for k, v in P.items():
        arr = np.asarray(v, dtype=float); W[k] = 1.0
        if len(arr) > a.maxPairs:
            W[k] = len(arr) / a.maxPairs; arr = arr[rng.choice(len(arr), a.maxPairs, replace=False)]
        A[k] = arr
    bonus = collections.defaultdict(float)
    if a.endLinks:
        for l in open(a.endLinks):
            if l[0] == "#" or not l.strip(): continue
            s = l.split(); bonus[frozenset(((s[0][:-1], s[0][-1]), (s[1][:-1], s[1][-1])))] = max(bonus[frozenset(((s[0][:-1], s[0][-1]), (s[1][:-1], s[1][-1])))], a.linkBonus)
    syn = {}
    if a.synteny:
        for l in open(a.synteny):
            if l[0] == "#" or not l.strip(): continue
            x, y, n = l.split(); syn[frozenset(((x[:-1], x[-1]), (y[:-1], y[-1])))] = int(n) * a.synBonus
    tel = {}
    if a.telomeres:
        for l in open(a.telomeres):
            if l[0] == "#": continue
            r = l.rstrip().split("\t"); tel[r[0]] = r[4]
    def tel_score(order):
        if not order or not tel: return 0.0
        s = 0.0
        c, o = order[0]; call = tel.get(c, "--")          # the first chain's left end faces out
        if (call[0] if o == "+" else call[1]) == "T": s += a.telBonus
        c, o = order[-1]; call = tel.get(c, "--")         # the last chain's right end faces out
        if (call[1] if o == "+" else call[0]) == "T": s += a.telBonus
        return s
    lowinfo = set(); stage2 = [False]
    if a.stats:
        lowinfo = { l.split("\t")[0] for l in open(a.stats) if l[0] != "#" and float(l.split("\t")[2]) < a.minWithin }
    def pair_score(i, oi, j, oj, between):
        '''i before j in the scaffold, orientations oi/oj, `between` bp between them'''
        k = (i, j) if i < j else (j, i)
        if k not in A or (not stage2[0] and (i in lowinfo or j in lowinfo)): return 0.0
        arr = A[k]; xi, xj = (arr[:, 0], arr[:, 1]) if i < j else (arr[:, 1], arr[:, 0])
        di = L[i] - xi if oi == "+" else xi     # distance from the read to i's end that faces j
        dj = xj if oj == "+" else L[j] - xj     # distance from j's end that faces i to the read
        return float(W[k] * logp(di + dj + between).sum())
    cache = {}
    def ps(i, oi, j, oj, between):
        b = int(round(between / 5000.0))        # memoise on 5 kb of spacing: enough for this decay
        key = (i, oi, j, oj, b, stage2[0])
        if key not in cache: cache[key] = pair_score(i, oi, j, oj, b * 5000.0)
        return cache[key]
    def total(order):
        s = 0.0; pos = 0; st = []
        for c, o in order: st.append(pos); pos += L[c] + a.gap
        for m in range(len(order)):
            for n in range(m + 1, len(order)):
                i, oi = order[m]; j, oj = order[n]
                s += ps(i, oi, j, oj, st[n] - st[m] - L[i])
            if m + 1 < len(order):
                i, oi = order[m]; j, oj = order[m + 1]
                e1 = (i, oi); e2 = (j, "-" if oj == "+" else "+")
                s += bonus.get(frozenset((e1, e2)), 0.0)
        if syn:     # junctions between consecutive informative chains, skipping low-information spacers
            inf = [ x for x in order if x[0] not in lowinfo ]
            for (i, oi), (j, oj) in zip(inf, inf[1:]):
                s += syn.get(frozenset(((i, oi), (j, "-" if oj == "+" else "+"))), 0.0)
        return s + tel_score(order)
    def links_to(c, cs):
        return sum(len(P.get((c, x) if c < x else (x, c), ())) for x in cs if x != c)
    out = open(a.out, "w")
    out.write("#group\trank\tchain\torientation\tlength\tlinks\torientation_llr\tposition_llr\tend_link\n")
    for g in sorted(gr):
        cs = sorted(gr[g], key=lambda c: -L[c])
        # a chain that is its whole group has nothing to link to, and is already in order
        linked = { e[0] for k in bonus for e in k }
        use = cs if len(cs) == 1 else [ c for c in cs if (L[c] >= a.minLen and links_to(c, cs) >= a.minPairs) or c in linked or "T" in tel.get(c, "") ]
        rest = [ c for c in cs if c not in use ]
        def climb(order, only_low=False):
            cur = total(order); improved = True; rounds = 0
            n = len(order)
            while improved and rounds < 200:
                improved = False; rounds += 1
                for size in range(1, n + 1):
                    for x in range(n - size + 1):
                        seg = order[x:x + size]; base = order[:x] + order[x + size:]
                        if only_low and any(c not in lowinfo for c, _ in seg): continue
                        rev = [ (c, "-" if o == "+" else "+") for c, o in reversed(seg) ]
                        for pos in range(len(base) + 1):
                            for sg in (seg, rev):
                                t = base[:pos] + sg + base[pos:]
                                if t == order: continue
                                s = total(t)
                                if s > cur + 1e-6: cur, order, improved = s, t, True; break
                            if improved: break
                        if improved: break
                    if improved: break
            return cur, order, rounds
        order = []
        for c in use:
            best = None
            for pos in range(len(order) + 1):
                for o in "+-":
                    t = order[:pos] + [(c, o)] + order[pos:]; s = total(t)
                    if best is None or s > best[0]: best = (s, t)
            order = best[1]
        cur, order, it = climb(order)
        greedy = cur; wins = 0
        if len(order) > 2:
            for st in range(a.starts):
                t = [ (use[i], "+-"[int(rng.integers(2))]) for i in rng.permutation(len(use)) ]
                s, t, _ = climb(t)
                if s > cur + 1e-6: cur, order = s, t; wins += 1
        if cur > greedy + 1e-6: print(f"# {g}: restarts improved the score by {cur - greedy:.1f}", file=sys.stderr)
        stage2[0] = True
        cur, order, _ = climb(order, only_low=True)
        n = len(order)
        for r, (c, o) in enumerate(order):
            flip = order[:r] + [(c, "-" if o == "+" else "+")] + order[r + 1:]
            ollr = cur - total(flip)
            base = order[:r] + order[r + 1:]
            pl = [ total(base[:q] + [(c, o2)] + base[q:]) for q in range(n) if q != r for o2 in "+-" ]
            pllr = cur - max(pl) if pl else float("inf")
            el = "."
            if r:
                pc, po = order[r - 1]
                el = "yes" if frozenset(((pc, po), (c, "-" if o == "+" else "+"))) in bonus else "no"
            out.write(f"{g}\t{r+1}\t{c}\t{o}\t{L[c]}\t{links_to(c, cs)}\t{ollr:.1f}\t{pllr:.1f}\t{el}\n")
        for c in rest: out.write(f"{g}\t.\t{c}\t.\t{L[c]}\t{links_to(c, cs)}\t.\t.\tunplaced_in_group\n")
        stage2[0] = False
        print(f"# {g}: {n} ordered, {len(rest)} not ordered, iterations {it}", file=sys.stderr)

if __name__ == "__main__":
    main()
