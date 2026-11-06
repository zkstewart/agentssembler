#!/usr/bin/env python3
'''
Hi-C check of gap fills: does each fill's allele-specific Hi-C agree with its own haplotype?

Pairs (pairs.tsv.gz, on the unitigs) are lifted through the gap-filled AGP. For every filled
gap that added sequence, reads in the fill whose mate lies on the same chromosome (this
scaffold or a homolog: same name before the last '_') are counted as own or homolog. The
same is counted for the --flank bp either side of the fill. A fill from the wrong haplotype
has a clearly lower own share than its flanks: flagged 'conflict' when it has >= --minPairs
informative pairs, an own share below 0.5, and below the flanks' share minus 0.25.
Fills with fewer informative pairs are 'no_signal' (typically homozygous sequence, where no
read maps uniquely).

Usage: gapfill_validate.py <gapfilled.agp> <gapfill.report.tsv> <pairs.tsv.gz> > validate.tsv
'''
import sys, gzip, bisect, argparse, collections

def key(s): return s.rsplit("_", 1)[0]

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("agp"); p.add_argument("report"); p.add_argument("pairs")
    p.add_argument("--flank", type=int, default=500000)
    p.add_argument("--minPairs", type=int, default=10)
    a = p.parse_args()

    place = collections.defaultdict(list) # unitig -> [(cStart, cEnd, scaf, sStart, strand)]
    fills = [] # (scaf, start, end, gap_id)
    report = { l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(a.report) if not l.startswith("#") }
    filled_units = {}
    for gid, r in report.items():
        if r[4] == "yes" and r[6] != ".":
            for x in r[6].split(";"): filled_units.setdefault(x[:-1], []).append(gid)
    pieces = collections.defaultdict(list)
    for l in open(a.agp):
        if l.startswith("#"): continue
        s = l.rstrip("\n").split("\t")
        if s[4] != "W": continue
        place[s[5]].append((int(s[6]), int(s[7]), s[0], int(s[1]), s[8]))
        pieces[s[0]].append((int(s[1]), int(s[2]), s[5]))
    # fill intervals: consecutive pieces of the fill's unitigs on the gap's scaffold
    for gid, r in report.items():
        if r[4] != "yes" or r[6] == ".": continue
        scaf = gid.split(":")[0]; us = { x[:-1] for x in r[6].split(";") }
        iv = [ (st, en) for st, en, u in pieces[scaf] if u in us ]
        # several fills may use the same unitig on one scaffold only if shared; take the
        # interval nearest the original gap position
        g = int(gid.split(":")[1])
        iv.sort(key=lambda x: abs(x[0] - g))
        if iv:
            # merge the pieces adjacent to the nearest one
            st, en = iv[0]
            for s2, e2 in sorted(iv):
                if s2 <= en + 1 and e2 >= st - 1: st, en = min(st, s2), max(en, e2)
            fills.append((scaf, st, en, gid, r[3]))
    byscaf = collections.defaultdict(list)
    for f in fills: byscaf[f[0]].append(f)
    for v in byscaf.values(): v.sort(key=lambda x: x[1])
    starts = { k: [ f[1] for f in v ] for k, v in byscaf.items() }

    def lift(u, x):
        for cs, ce, sc, ss, st in place.get(u, ()):
            if cs <= x <= ce: return sc, (ss + (x - cs) if st == "+" else ss + (ce - x))
        return None, None
    def which(sc, x, flank):
        '''index of fill containing x (flank=0), or within flank of it'''
        if sc not in byscaf: return []
        out = []
        i = bisect.bisect_right(starts[sc], x + flank) - 1
        while i >= 0 and byscaf[sc][i][2] + flank >= x - flank * 3:
            f = byscaf[sc][i]
            if flank == 0 and f[1] <= x <= f[2]: out.append(i)
            elif flank and (f[1] - flank <= x < f[1] or f[2] < x <= f[2] + flank): out.append(i)
            i -= 1
        return out
    cnt = collections.defaultdict(lambda: [0, 0, 0, 0]) # (scaf, i) -> fill own, fill homolog, flank own, flank homolog
    with gzip.open(a.pairs, "rt") as fh:
        for l in fh:
            ua, pa, ub, pb = l.split("\t")
            sa, xa = lift(ua, int(pa)); sb, xb = lift(ub, int(pb))
            if sa is None or sb is None or key(sa) != key(sb): continue
            for s1, x1, s2 in ((sa, xa, sb), (sb, xb, sa)):
                for i in which(s1, x1, 0):
                    cnt[(s1, i)][0 if s2 == s1 else 1] += 1
                for i in which(s1, x1, a.flank):
                    cnt[(s1, i)][2 if s2 == s1 else 3] += 1
    print("#gap_id\tmethod\tscaffold\tfill_start\tfill_end\tfill_bp\tfill_own\tfill_homolog\tfill_own_share\tflank_own\tflank_homolog\tflank_own_share\tverdict")
    summary = collections.Counter()
    for sc in sorted(byscaf):
        for i, f in enumerate(byscaf[sc]):
            fo, fh_, ko, kh = cnt[(sc, i)]
            fs = fo / (fo + fh_) if fo + fh_ else float("nan"); ks = ko / (ko + kh) if ko + kh else float("nan")
            if fo + fh_ < a.minPairs: v = "no_signal"
            elif fs < 0.5 and fs < ks - 0.25: v = "conflict"
            else: v = "consistent"
            summary[v] += 1
            print(f"{f[3]}\t{f[4]}\t{sc}\t{f[1]}\t{f[2]}\t{f[2]-f[1]+1}\t{fo}\t{fh_}\t{fs:.3f}\t{ko}\t{kh}\t{ks:.3f}\t{v}")
    print("# " + ", ".join(f"{k}: {n}" for k, n in sorted(summary.items())), file=sys.stderr)

if __name__ == "__main__":
    main()
