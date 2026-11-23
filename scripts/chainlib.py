# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Shared helpers: lift unitig coordinates through a chain AGP; sort keys for chains that survive renumbering.'''
import collections, gzip, re
def read_agp(path):
    '''unitig -> list of (chain, chain_start0, comp_start0, comp_end0_exclusive, strand)'''
    P = collections.defaultdict(list); L = {}
    for l in open(path):
        if l[0] == "#": continue
        s = l.rstrip("\n").split("\t")
        L[s[0]] = max(L.get(s[0], 0), int(s[2]))
        if s[4] in ("U", "N"): continue
        P[s[5]].append((s[0], int(s[1]) - 1, int(s[6]) - 1, int(s[7]), s[8]))
    return P, L
def lift(P, u, x):
    '''unitig 0-based position -> (chain, position), or (None, None) outside the kept range'''
    for c, cs, a, b, st in P.get(u, ()):
        if a <= x < b: return c, (cs + x - a if st == "+" else cs + b - 1 - x)
    return None, None
def lift_strand(P, u, x):
    '''as lift, plus the unitig's strand in that chain: (chain, position, strand), or (None, None, None)'''
    for c, cs, a, b, st in P.get(u, ()):
        if a <= x < b: return c, (cs + x - a if st == "+" else cs + b - 1 - x), st
    return None, None, None
def stable_key(parts, c):
    '''A sort key for chain c that renumbering chains does not change: its unitig names, sorted, then any copy suffix (.gNN).
    parts: base chain -> [(unitig, strand), ...] from the chain AGP. Ties broken by this key give the same choice when a re-run
    renumbers chains (connectors after a reset, say), which chain names would not.'''
    b, _, suf = c.partition(".g")
    return (tuple(sorted(u for u, _ in parts.get(b, ()))) or (b,), suf)
def stable_end(parts, e): return (stable_key(parts, e[:-1]), e[-1])
def stable_path(parts, p): return tuple(stable_end(parts, z) for z in p)
def gene_hits(table, P):
    '''compleasm full_table_busco_format.tsv on the unitigs -> [(gene, status, chain or None, chain position)], one row per copy
    (Single, Duplicated, Fragmented), each copy placed by its midpoint. P from read_agp. A copy on a unitig no chain kept
    (redundant or repeat) gets chain None.'''
    out = []
    for l in open(table):
        if l[0] == "#": continue
        s = l.rstrip("\n").split("\t")
        if len(s) < 5 or s[1] not in ("Single", "Complete", "Duplicated", "Fragmented"): continue
        c, x = lift(P, s[2], (int(s[3]) + int(s[4])) // 2)
        out.append((s[0], s[1], c, x))
    return out
def read_windows(path):
    '''unitig_depth.py windows file -> chain -> [(start, end, total, unique)], sorted'''
    W = collections.defaultdict(list)
    for l in open(path):
        if l[0] == "#": continue
        s = l.split("\t"); W[s[0]].append((int(s[1]), int(s[2]), float(s[3]), float(s[4])))
    for v in W.values(): v.sort()
    return W
def core_runs(wins, lo, minRun):
    '''Runs of adjacent windows whose unique-read depth (MAPQ >= 20 reads) is >= lo, at least minRun bp long:
    the core of a collapsed hub, where both haplotypes' reads land and all of them are unique. -> [(start, end, median unique, median total)].
    Total depth alone misses these: the hub's arms read 1x with reads that are MAPQ 0, which pulls the chain's median below the shared class.'''
    runs, cur = [], None
    for a, b, t, u in wins:
        if u >= lo and cur and a == cur[1]: cur[1] = b; cur[2].append(u); cur[3].append(t)
        elif u >= lo: cur = [a, b, [u], [t]]; runs.append(cur)
        else: cur = None
    med = lambda v: sorted(v)[len(v) // 2]
    return [(a, b, med(us), med(ts)) for a, b, us, ts in runs if b - a >= minRun]
def first_undecided(tops, decided):
    '''tops: [((reads, hic), solution), ...] in preference order. The first solution that no other decides against, where reads the same
    on both sides cancel (the comparison is then on Hi-C alone) and otherwise reads + Hi-C are compared. Against one with the same reads only
    the most Hi-C matters; against the rest only the most reads + Hi-C. The solution with the most reads + Hi-C is never decided against.'''
    mh, mt = {}, {}
    for (r, h), _ in tops: mh[r] = max(mh.get(r, h), h); mt[r] = max(mt.get(r, r + h), r + h)
    def beaten(r, h): return decided(mh[r], h) or any(decided(v, r + h) for r2, v in mt.items() if r2 != r)
    return next(t for (r, h), t in tops if not beaten(r, h))
def gfa_overlap(cig):
    '''GFA L-line overlap CIGAR -> bp of the link's overlap on the segment entered: M, I, S, = and X consume it. hifiasm writes one op
    ("12345M"); other assemblers can write several ("10M1I5M"), which int(cig[:-1]) misreads. "*" (no overlap given) -> 0.'''
    cig = cig.strip()
    if cig in ("*", ""): return 0
    return sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cig) if op in "MIS=X")
FLIP = {"+": "-", "-": "+"}
def end_overlaps(gfa):
    '''(unitig, end) -> longest overlap of any GFA link at that end; end "+" is the 3' end, "-" the 5' end'''
    ov = collections.defaultdict(int)
    for l in open(gfa):
        if l[0] != "L": continue
        s = l.rstrip("\n").split("\t"); o = gfa_overlap(s[5])
        k = (s[1], s[2]); ov[k] = max(ov[k], o); k = (s[3], FLIP[s[4]]); ov[k] = max(ov[k], o)
    return ov
def chain_interior(parts, L, ov, chains):
    '''chain -> (start, end) of its sequence outside graph overlaps at both ends. parts: chain -> [(unitig, strand), ...] in order'''
    iv = {}
    for c in chains:
        if not parts.get(c): continue
        u, o = parts[c][-1]; right = ov[(u, o)]
        u, o = parts[c][0]; left = ov[(u, FLIP[o])]
        iv[c] = (left, L[c] - right)
    return iv
def interior_cover(paf, iv, q20=False):
    '''bp of primary alignment inside each chain's interior (iv from chain_interior); with q20, also bp from MAPQ >= 20 alignments'''
    cov, q = collections.Counter(), collections.Counter()
    for l in gzip.open(paf, "rt"):
        s = l.split("\t", 12); c = s[5]
        if c not in iv or "tp:A:P" not in l: continue
        lo, hi = max(iv[c][0], int(s[7])), min(iv[c][1], int(s[8]))
        if hi > lo:
            cov[c] += hi - lo
            if int(s[11]) >= 20: q[c] += hi - lo
    return (cov, q) if q20 else cov
def utg_resolver(agp):
    '''A function that turns "utg:<unitig>[.gNN][+|-]" into the name of the chain holding that unitig in this AGP, keeping any copy
    suffix; an end sign is the end the unitig's own end faces, so it survives a chain written the other way round. Other names pass
    unchanged. Chain numbers change with any change upstream (TRAPS T16); a unitig name does not, so an edit written this way can be
    carried across a reset. Stops on a unitig in no chain or in more than one.'''
    where = collections.defaultdict(list)
    for l in open(agp):
        if l[0] == "#": continue
        s = l.rstrip("\n").split("\t")
        if s[4] not in "NU": where[s[5]].append((s[0], s[8]))
    def res(tok):
        m = re.fullmatch(r"utg:(.+?)((?:\.g\d+)?)([+-]?)", tok)
        if not m: return tok
        u, suf, end = m.groups(); w = where.get(u, [])
        if len({c for c, _ in w}) != 1:
            raise SystemExit(f"{tok}: unitig {u} is in {len({c for c, _ in w})} chains of {agp}" + (f" ({', '.join(sorted({c for c, _ in w}))})" if w else ""))
        c, st = w[0]
        return c + suf + (end if not end or st == "+" else FLIP[end])
    return res
