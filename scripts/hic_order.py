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

Search: greedy insertion of chains longest first (every position and orientation), then local moves
(move any contiguous block to any position in either orientation, which includes reinserting one chain,
flipping and reversing a segment) to convergence. This is repeated from random orders, because a
single block move can be blocked when two blocks must move together (e.g. both carry a telomere bonus);
the best result over all starts is kept. Restarts stop after --patience in a row that do not beat the best (at most
--starts): the search spends effort where restarts still find better orders, and a group whose restarts keep failing
stops early (on 292, chr2 with 52 and 53 chains per group: 6 min instead of about 2 h, same order). Each group's restarts draw from a generator seeded by
the group's name, so a group's result does not depend on the others. The log gives, per group, the restarts run, those
that improved, and those that reached the best stage-1 order (in either direction): few reaching it means local optima, and an
order to read with its LLRs. Chains shorter than --minLen or with
fewer than --minPairs links to the group are appended unordered as 'unplaced_in_group'.

Output columns: group, rank, chain, orientation, length, links to group, orientation_llr, position_llr,
end_link (whether the junction with the previous chain has a link). The LLRs are the best alternative's
loss in log-likelihood: small values mean Hi-C does not settle that chain's orientation or position.

Usage: hic_order.py groups.tsv pairs.tsv.gz[,more.tsv.gz] chains.fai out.tsv [--endLinks links.tsv] [--decayIn decay.tsv]
       hic_order.py groups.tsv pairs.tsv.gz[,more] chains.fai - --decayOut decay.tsv
'''
import gzip, math, argparse, collections, sys, zlib
import numpy as np

EDGES = np.logspace(2, 9, 141)
def decay_hist(pairs_within):
    """Counts of within-chain pair distances in the decay model's log bins."""
    d = np.asarray(pairs_within, dtype=float); d = d[d > 0]
    return np.histogram(d, bins=EDGES)[0]
def decay_model(pairs_within, bg=1e-3, lengths=None, hist=None):
    """log P(contact | distance), up to a constant, from pairs within chains. With lengths (the chains the pairs came from), each bin is
    divided by the bp of chain at which that distance can be seen, sum(max(L - d, 0)) averaged over the bin: near the chain lengths few
    positions are that far apart, so raw counts fall away and understate long-range contact. The average is exact per chain (0 for
    L <= a, L - (a + b) / 2 for L >= b, (L - a)^2 / 2(b - a) between), so a bin the longest chain ends in is never divided by ~0."""
    edges = EDGES
    h = decay_hist(pairs_within) if hist is None else np.asarray(hist)
    w = np.diff(edges); dens = (h + 0.5) / w / max(h.sum(), 1)
    if lengths is not None:
        Ls = np.asarray(lengths, dtype=float)[:, None]; lo, hi = edges[None, :-1], edges[None, 1:]
        opp = np.where(Ls >= hi, Ls - (lo + hi) / 2, np.where(Ls > lo, (Ls - lo) ** 2 / (2 * (hi - lo)), 0.0)).sum(axis=0)
        dens = dens / np.maximum(opp, 1.0)
    c = np.sqrt(edges[:-1] * edges[1:]); lc, ld = np.log(c), np.log(dens)
    # beyond the observed range, extend with the slope of the last well-populated decade (bins with >= 50 pairs, or any pairs in a sparse
    # library); contact never rises with distance, so the slope is at most 0
    ok = np.where(h >= 50)[0]
    if len(ok) < 2: ok = np.where(h > 0)[0]
    if len(ok) == 0: return lambda x: np.zeros_like(np.asarray(x, dtype=float))
    last = ok[-1]; first = ok[max(0, len(ok) - 20)]
    slope = min((ld[last] - ld[first]) / (lc[last] - lc[first] + 1e-9), 0.0)
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
    p.add_argument("--starts", type=int, default=40, help="most random restarts per group")
    p.add_argument("--patience", type=int, default=3, help="stop after this many restarts in a row that do not beat the best")
    p.add_argument("--stats", help="chain_stats.tsv")
    p.add_argument("--minWithin", type=float, default=700)
    p.add_argument("--telomeres", help="telomeres.py output for the chains")
    p.add_argument("--telBonus", type=float, default=20000.0)
    p.add_argument("--gap", type=int, default=20000, help="assumed gap between adjacent chains, bp")
    p.add_argument("--minLen", type=int, default=20000)
    p.add_argument("--minPairs", type=int, default=10)
    p.add_argument("--maxPairs", type=int, default=3000, help="subsample per chain pair (scaled back in the score)")
    p.add_argument("--decayNorm", action="store_true", help="divide the distance histogram by the bp at which each distance can be seen")
    p.add_argument("--decayOut", help="write the within-chain distance histogram (and the chains it came from) here and stop: the decay "
                   "model is the same for every group, so stage 09 builds it once instead of once per chromosome job")
    p.add_argument("--decayIn", help="read the histogram --decayOut wrote, instead of collecting within-chain pairs")
    a = p.parse_args()
    L = { l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.fai) }
    gr = collections.defaultdict(list); of = {}
    for l in open(a.groups):
        if l[0] == "#": continue
        c, n, g = l.rstrip().split("\t")
        if g != "unassigned": gr[g].append(c); of[c] = g
    within = []; wchains = set(); P = collections.defaultdict(list); H = None
    if a.decayIn:
        for l in open(a.decayIn):
            k, *v = l.rstrip("\n").split("\t")
            if k == "hist": H = [int(x) for x in v]
            elif k == "chain": wchains.add(v[0])
    for l in (l for f in a.pairs.split(",") for l in gzip.open(f, "rt")):
        x, u, y, v = l.split("\t")
        if x == y:
            if H is None: within.append(abs(int(u) - int(v))); wchains.add(x)
            continue
        if a.decayOut: continue
        if x in of and y in of and of[x] == of[y]:
            u = int(u); v = int(v)
            if x < y: P[(x, y)].append((u, v))
            else: P[(y, x)].append((v, u))
    if a.decayOut:
        with open(a.decayOut, "w") as o:
            o.write("hist\t" + "\t".join(str(int(x)) for x in decay_hist(within)) + "\n")
            for c in sorted(wchains): o.write(f"chain\t{c}\n")
        return
    logp = decay_model(within, lengths=[L[c] for c in wchains] if a.decayNorm else None, hist=H)
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
    # partners with pairs, per stage: a chain pair without pairs scores exactly 0.0, so skipping it leaves every total
    # unchanged (same terms, same order of addition) and the search identical; with many chains most pairs are empty
    NB = [collections.defaultdict(set), collections.defaultdict(set)]
    for i, j in A:
        NB[1][i].add(j); NB[1][j].add(i)
        if i not in lowinfo and j not in lowinfo: NB[0][i].add(j); NB[0][j].add(i)
    def total(order):
        s = 0.0; pos = 0; st = []
        for c, o in order: st.append(pos); pos += L[c] + a.gap
        nb_stage = NB[1] if stage2[0] else NB[0]
        for m in range(len(order)):
            i, oi = order[m]; nb = nb_stage[i]; sm = st[m] + L[i]
            if nb:
                for n in range(m + 1, len(order)):
                    j, oj = order[n]
                    if j in nb: s += ps(i, oi, j, oj, st[n] - sm)
            if m + 1 < len(order):
                i, oi = order[m]; j, oj = order[m + 1]
                e1 = (i, oi); e2 = (j, "-" if oj == "+" else "+")
                s += bonus.get(frozenset((e1, e2)), 0.0)
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
        greedy = cur; wins = 0; runs = []; miss = 0
        canon = lambda o: min(tuple(o), tuple((c, "-" if x == "+" else "+") for c, x in reversed(o)))   # an order and its reverse are one scaffold
        if len(order) > 2:
            rg = np.random.default_rng(zlib.crc32(g.encode()))
            for st in range(a.starts):
                t = [ (use[i], "+-"[int(rg.integers(2))]) for i in rg.permutation(len(use)) ]
                s, t, _ = climb(t); runs.append(canon(t))
                if s > cur + 1e-6: cur, order = s, t; wins += 1; miss = 0
                else: miss += 1
                if miss >= a.patience: break
        if cur > greedy + 1e-6: print(f"# {g}: restarts improved the score by {cur - greedy:.1f}", file=sys.stderr)
        if runs: print(f"# {g}: restarts {len(runs)} of at most {a.starts} (patience {a.patience}), improved {wins}, reached the best order {sum(r == canon(order) for r in runs)}", file=sys.stderr)
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
