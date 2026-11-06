#!/usr/bin/env python3
'''
Haplotype-aware gap filling from long reads mapped competitively to all haplotypes.

The reads must be mapped (minimap2, PAF, no secondary alignments) to the whole assembly: every
pseudochromosome of every haplotype plus the unplaced sequences, so that each read lands on the
haplotype it matches best and its MAPQ says how uniquely. A gap is filled only from reads that

  1. anchor >= --minAnchor bp in BOTH flanks with MAPQ >= --minMapq (the read belongs to this
     haplotype at both ends, not to a homolog), reaching the gap edge (within --slack bp), in a
     consistent order and orientation; and
  2. agree with each other: the fill each read spells between fixed cut points in the flanks is
     >= --minIdentity identical to the medoid fill for >= --minAgree of the reads (and
     >= --minReads reads).

The fill is the medoid read's sequence between the cut points. The cut points sit at the flank
ends, moved inwards (trimming the flank by <= --maxTrim bp) only where the reads stop matching a
flank before its end. Where the reads show the flanks overlap, the overlap is removed instead.
Each fill is checked by realigning the spanning reads to the new local sequence: >= --minReads
must align through the fill +/- 1 kb without an indel > 50 bp.

Reads anchored uniquely at a gap edge that continue elsewhere (another scaffold, or far away on
this one) are counted; when they outnumber the spanning reads the join is reported as
'read_conflict' (a candidate misjoin for curation) and not filled.

Outputs (--out prefix):
  .agp            the input AGP with the filled gaps replaced (fill components rf_*)
  .fills.fasta    the fill sequences (components of the AGP besides the unitigs)
  .report.tsv     one row per gap of the input AGP
  .membership.tsv fill component -> scaffold, how = readfill
  .excised.tsv    components removed as false duplications (see below)
  .links.tsv      for every gap edge, where uniquely anchored reads continue instead (curation hints)

False duplications. When uniquely anchored reads leave the left edge of one gap and continue, in
the same direction, at the right edge of a later gap on the same scaffold, the segment between is
not in this haplotype. It is excised if >= --minDupCover of it also has a copy elsewhere in
--target at >= --minDupIdentity (a homozygous or other-haplotype copy placed twice); the two
gaps are then bridged from those reads like any other gap (or left as one gap). Segments without
such a copy are unique sequence in the wrong place: they are only reported.

Usage: readfill.py --agp in.agp --fasta in.fasta --paf reads.paf.gz --reads reads.fasta --out prefix
'''
import os, re, sys, gzip, bisect, argparse, tempfile, statistics, subprocess, collections, concurrent.futures
try:
    import parasail
except ImportError:
    parasail = None

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def rc(s): return s.translate(COMP)[::-1]

Aln = collections.namedtuple("Aln", "q ql qs qe strand t tl ts te mapq oqs oqe")
LAln = collections.namedtuple("LAln", "q qs qe strand t ts te cigar")

def parse_paf(line):
    s = line.split("\t", 12)
    ql, qs, qe = int(s[1]), int(s[2]), int(s[3])
    # oriented query coordinates: positions on the read as it runs along the target's + strand
    oqs, oqe = (qs, qe) if s[4] == "+" else (ql - qe, ql - qs)
    return Aln(s[0], ql, qs, qe, s[4], s[5], int(s[6]), int(s[7]), int(s[8]), int(s[11]), oqs, oqe)

def read_agp(path):
    rows = collections.defaultdict(list)
    for l in open(path):
        if l.startswith("#"): continue
        s = l.rstrip("\n").split("\t")
        rows[s[0]].append(s)
    return rows

class Fasta:
    '''random access through the .fai'''
    def __init__(self, fasta):
        self.f = open(fasta, "rb"); self.idx = {}
        for l in open(fasta + ".fai"):
            s = l.split("\t"); self.idx[s[0]] = (int(s[1]), int(s[2]), int(s[3]), int(s[4]))
    def get(self, name, start=0, end=None): # 0-based half-open
        ln, off, lb, by = self.idx[name]
        start = max(0, start); end = ln if end is None else min(ln, end)
        if end <= start: return ""
        a = off + (start // lb) * by + start % lb
        b = off + (end // lb) * by + end % lb
        return os.pread(self.f.fileno(), b - a, a).decode().replace("\n", "").replace("\r", "").upper()

CIG = re.compile(r"(\d+)([MIDNSHP=X])")
def run_minimap2(mm2, preset, ref, qry, extra=()):
    out = subprocess.run([mm2, "-x", preset, "-c", "-t", "1", "--secondary=no", *extra, ref, qry],
                         capture_output=True, text=True, check=True).stdout
    alns = []
    for l in out.splitlines():
        s = l.split("\t")
        cg = next(x[5:] for x in s[12:] if x.startswith("cg:Z:"))
        alns.append(LAln(s[0], int(s[2]), int(s[3]), s[4], s[5], int(s[7]), int(s[8]),
                         [ (int(n), op) for n, op in CIG.findall(cg) ]))
    return alns

def q_at(a, t):
    '''query position at the reference boundary t (bases [ts, t) of the reference lie before it)'''
    qp, tp = a.qs, a.ts
    if t <= tp: return qp - (tp - t)
    for n, op in a.cigar:
        if op in "M=X":
            if tp + n >= t: return qp + (t - tp)
            qp += n; tp += n
        elif op in "DN":
            if tp + n >= t: return qp
            tp += n
        elif op == "I": qp += n
    return qp + (t - tp)

def t_at(a, q):
    '''reference position at the query boundary q'''
    qp, tp = a.qs, a.ts
    if q <= qp: return tp - (qp - q)
    for n, op in a.cigar:
        if op in "M=X":
            if qp + n >= q: return tp + (q - qp)
            qp += n; tp += n
        elif op in "DN": tp += n
        elif op == "I":
            if qp + n >= q: return tp
            qp += n
    return tp + (q - qp)

def max_indel(a, t0, t1):
    '''largest indel of alignment a within reference interval [t0, t1)'''
    tp, m = a.ts, 0
    for n, op in a.cigar:
        if op in "M=X": tp += n
        elif op in "DN":
            if tp < t1 and tp + n > t0: m = max(m, n)
            tp += n
        elif op == "I" and t0 <= tp <= t1: m = max(m, n)
    return m

def identity(a, b, mm2, tmp):
    '''matching bases / longer length'''
    la, lb = len(a), len(b)
    if max(la, lb) == 0: return 1.0
    if min(la, lb) == 0: return 0.0
    if min(la, lb) / max(la, lb) < 0.9: return min(la, lb) / max(la, lb) * 0.9 # cannot be similar
    if max(la, lb) <= 3000 and parasail is not None:
        r = parasail.nw_stats_striped_16(a, b, 5, 2, parasail.dnafull)
        return r.matches / max(la, lb)
    if max(la, lb) <= 300: # no parasail: plain edit distance
        prev = list(range(lb + 1))
        for i in range(1, la + 1):
            cur = [i] + [0] * lb
            for j in range(1, lb + 1):
                cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (a[i - 1] != b[j - 1]))
            prev = cur
        return 1 - prev[lb] / max(la, lb)
    fa, fb = os.path.join(tmp, "a.fa"), os.path.join(tmp, "b.fa")
    open(fa, "w").write(f">a\n{a}\n"); open(fb, "w").write(f">b\n{b}\n")
    out = subprocess.run([mm2, "-x", "asm20", "-c", "-t", "1", fa, fb], capture_output=True, text=True, check=True).stdout
    m = sum(int(s.split("\t")[9]) for s in out.splitlines() if "tp:A:P" in s and s.split("\t")[4] == "+")
    return min(1.0, m / max(la, lb))

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--agp", required=True, help="AGP over the unitigs (e.g. the gapfill.py output)")
    p.add_argument("--fasta", required=True, help="the AGP's scaffold FASTA (with .fai)")
    p.add_argument("--paf", required=True, help="reads mapped to the whole assembly + unplaced (minimap2, --secondary=no)")
    p.add_argument("--reads", required=True, help="the reads (FASTA/FASTQ)")
    p.add_argument("--out", required=True)
    p.add_argument("--gapfillReport", help="gapfill.py report: annotate each gap with its graph outcome")
    p.add_argument("--minimap2", default="minimap2")
    p.add_argument("--seqkit", default="seqkit")
    p.add_argument("--preset", default="lr:hq", help="minimap2 preset for the reads (lr:hq corrected/Q20 ONT, map-hifi, map-ont)")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--minMapq", type=int, default=20)
    p.add_argument("--minAnchor", type=int, default=5000)
    p.add_argument("--slack", type=int, default=1000, help="anchors must reach within this of the gap edge")
    p.add_argument("--maxTrim", type=int, default=2000)
    p.add_argument("--maxOverlap", type=int, default=50000, help="largest flank overlap accepted from the reads")
    p.add_argument("--window", type=int, default=60000, help="flank length used for the local realignment")
    p.add_argument("--minReads", type=int, default=2)
    p.add_argument("--minAgree", type=float, default=0.75)
    p.add_argument("--minIdentity", type=float, default=0.98)
    p.add_argument("--target", help="the FASTA the reads were mapped to (assembly + unplaced): enables excision of "
                   "segments that reads skip and that are duplicated elsewhere")
    p.add_argument("--maxExcise", type=int, default=2000000, help="longest segment that may be excised")
    p.add_argument("--minDupCover", type=float, default=0.9, help="excise only if >= this fraction of the segment ...")
    p.add_argument("--minDupIdentity", type=float, default=0.99, help="... has a copy elsewhere at >= this identity")
    a = p.parse_args()

    rows = read_agp(a.agp)
    gaps = [] # (scaf, row index, g0, g1) 0-based half-open
    for sc, rs in rows.items():
        for i, r in enumerate(rs):
            if r[4] in ("N", "U") and 0 < i < len(rs) - 1 and rs[i-1][4] not in ("N", "U") and rs[i+1][4] not in ("N", "U"):
                gaps.append((sc, i, int(r[1]) - 1, int(r[2])))
    gid = { (g[0], g[1]): f"{g[0]}:{g[2]+1}" for g in gaps }
    byscaf = collections.defaultdict(list)
    for g in gaps: byscaf[g[0]].append(g)
    for v in byscaf.values(): v.sort(key=lambda g: g[2])
    starts = { k: [ g[2] for g in v ] for k, v in byscaf.items() }
    print(f"# {len(gaps)} gaps", file=sys.stderr)

    def near(t, ts, te):
        '''gaps of scaffold t whose edge this alignment touches'''
        if t not in byscaf: return []
        v, out = byscaf[t], []
        i = bisect.bisect_right(starts[t], te + a.slack) - 1
        while i >= 0 and v[i][3] >= ts - a.slack - 10**6:
            g = v[i]
            if (te >= g[2] - a.slack and ts <= g[2] - a.minAnchor) or (ts <= g[3] + a.slack and te >= g[3] + a.minAnchor):
                out.append(g)
            i -= 1
        return out

    # pass 1: reads touching a gap edge; pass 2: all their alignments
    op = (lambda f: gzip.open(f, "rt")) if a.paf.endswith(".gz") else open
    touch = set()
    with op(a.paf) as fh:
        for l in fh:
            x = parse_paf(l)
            if near(x.t, x.ts, x.te): touch.add(x.q)
    alns = collections.defaultdict(list)
    with op(a.paf) as fh:
        for l in fh:
            q = l.split("\t", 1)[0]
            if q in touch: alns[q].append(parse_paf(l))
    print(f"# {len(touch)} reads touch a gap edge", file=sys.stderr)

    # classify reads per gap
    span = collections.defaultdict(dict)     # gap -> read -> (strand, mapq)
    lowq = collections.defaultdict(set)
    elsewhere = collections.defaultdict(lambda: collections.Counter()) # gap -> target -> reads
    deadend = collections.Counter()
    links = []                               # gap edge -> where uniquely anchored reads continue instead
    skiplinks = []                           # ... the subset continuing along the same scaffold
    for q, xs in alns.items():
        for g in { g for x in xs for g in near(x.t, x.ts, x.te) }:
            sc, _, g0, g1 = g
            L = [ x for x in xs if x.t == sc and x.te >= g0 - a.slack and x.ts <= g0 - a.minAnchor ]
            R = [ x for x in xs if x.t == sc and x.ts <= g1 + a.slack and x.te >= g1 + a.minAnchor ]
            best = None
            for l in L:
                for r in R:
                    if l.strand != r.strand: continue
                    if l is r or -a.maxOverlap <= r.oqs - l.oqe <= 10**6:
                        mq = min(l.mapq, r.mapq)
                        if best is None or mq > best[1]: best = (l.strand, mq)
            if best:
                if best[1] >= a.minMapq: span[g][q] = best
                else: lowq[g].add(q)
                continue
            # not spanning: does a unique edge anchor continue elsewhere?
            for side, A in (("L", L), ("R", R)):
                for x in A:
                    if x.mapq < a.minMapq: continue
                    lo, hi = (x.oqe, x.ql) if side == "L" else (0, x.oqs)
                    if hi - lo < a.minAnchor: continue
                    other = [ y for y in xs if y is not x and y.mapq >= a.minMapq
                              and min(hi, y.oqe) - max(lo, y.oqs) >= a.minAnchor ]
                    if other:
                        y = max(other, key=lambda y: y.qe - y.qs)
                        elsewhere[g][y.t if y.t != sc else f"{sc}:{y.ts//50000*50}kb"] += 1
                        # where the read joins y, and whether it runs along y in the same direction
                        same = y.strand == x.strand
                        jpos = (y.ts if same else y.te) if side == "L" else (y.te if same else y.ts)
                        links.append((gid[(sc, g[1])], side, q, y.t, jpos, "same" if same else "opposite", min(x.mapq, y.mapq)))
                        if y.t == sc and same: skiplinks.append((g, side, q, jpos, x.strand, min(x.mapq, y.mapq)))
                    else: deadend[g] += 1

    # Joins that skip a segment: reads leave the left edge of gap g1 and continue, in the same
    # direction, at the right edge of a later gap g2 (and/or the reverse). The segment between is
    # excised when it also has a near-identical copy elsewhere (a false duplication).
    skip = collections.defaultdict(dict)
    for g, side, q, j, strand, mq in skiplinks:
        sc = g[0]
        for g2 in byscaf[sc]:
            if side == "L" and g2[2] > g[2] and abs(j - g2[3]) <= a.slack: skip[(g, g2)][q] = (strand, mq)
            if side == "R" and g2[2] < g[2] and abs(j - g2[2]) <= a.slack: skip[(g2, g)][q] = (strand, mq)
    skipc = sorted(((k, v) for k, v in skip.items() if len(v) >= a.minReads and k[1][2] - k[0][3] <= a.maxExcise),
                   key=lambda kv: -len(kv[1]))
    dup = {}
    if skipc and a.target:
        d = tempfile.mkdtemp(prefix="readfill.excise.", dir=os.path.dirname(os.path.abspath(a.out)))
        asm0 = Fasta(a.fasta)
        # 10 kb chunks: a copy split across places, or one a whole-segment alignment ranks as a
        # weak secondary, still counts
        own, chunks = {}, collections.defaultdict(list)
        with open(os.path.join(d, "seg.fa"), "w") as fh:
            for (g1, g2), _ in skipc:
                k = gid[(g1[0], g1[1])]; own[k] = (g1[0], g1[3], g2[2])
                seq = asm0.get(g1[0], g1[3], g2[2])
                for o in range(0, len(seq), 10000):
                    c = seq[o:o + 10000]
                    if len(c) - c.count("N") >= 1000:
                        fh.write(f">{k}|{o}\n{c}\n"); chunks[k].append((o, len(c) - c.count("N")))
        out = subprocess.run([a.minimap2, "-x", "asm5", "-c", "-t", str(a.threads), "--secondary=yes", "-N", "50", "-p", "0.5",
                              a.target, os.path.join(d, "seg.fa")], capture_output=True, text=True, check=True).stdout
        cov = collections.defaultdict(list)
        for l in out.splitlines():
            x = l.split("\t"); k, o = x[0].split("|")
            sc, s0, s1 = own[k]
            if x[5] == sc and int(x[7]) < s1 and int(x[8]) > s0: continue # itself
            if int(x[9]) / int(x[10]) >= a.minDupIdentity: cov[x[0]].append((int(x[2]), int(x[3])))
        for k, cs in chunks.items():
            c = 0
            for o, n in cs:
                end = 0
                for x0, x1 in sorted(cov[f"{k}|{o}"]):
                    if x1 > end: c += x1 - max(x0, end); end = x1
            dup[k] = min(1.0, c / max(1, sum(n for _, n in cs)))
        with open(a.out + ".skipped.tsv", "w") as fh:
            fh.write("#join_id\tscaffold\tstart\tend\tskip_reads\tdup_cover\n")
            for (g1, g2), v in skipc:
                k = gid[(g1[0], g1[1])]
                fh.write(f"{k}\t{g1[0]}\t{g1[3]+1}\t{g2[2]}\t{len(v)}\t{dup[k]:.3f}\n")
        subprocess.run(["rm", "-rf", d])
    joins, inside, excise = [], set(), {}
    for (g1, g2), v in skipc:
        if dup.get(gid[(g1[0], g1[1])], 0) < a.minDupCover: continue
        rr = set(range(g1[1], g2[1] + 1))
        if any((g1[0], i) in inside for i in rr): continue
        inside |= { (g1[0], i) for i in rr }
        J = (g1[0], g1[1], g2[1], g1[2], g2[3]); joins.append(J); span[J] = v; excise[J] = (g1, g2)
    for g in gaps:
        if (g[0], g[1]) not in inside: joins.append((g[0], g[1], g[1], g[2], g[3])); span[joins[-1]] = span[g]
    for J in joins:
        g = next(g for g in gaps if (g[0], g[1]) == (J[0], J[1]))
        if J not in excise: lowq[J] = lowq[g]; elsewhere[J] = elsewhere[g]; deadend[J] = deadend[g]
    print(f"# {len(skipc)} segments skipped by >= {a.minReads} reads; {len(excise)} duplicated elsewhere -> excised", file=sys.stderr)

    # read sequences needed for the fills
    cand = [ J for J in joins if len(span[J]) >= a.minReads ]
    need = { q for g in cand for q in span[g] }
    tmp = tempfile.mkdtemp(prefix="readfill.", dir=os.path.dirname(os.path.abspath(a.out)))
    names = os.path.join(tmp, "names.txt"); rf = os.path.join(tmp, "reads.fasta")
    open(names, "w").write("".join(q + "\n" for q in sorted(need)))
    if need:
        subprocess.run(f"{a.seqkit} grep -j {a.threads} -f {names} {a.reads} | {a.seqkit} seq -w 0 > {rf} && {a.seqkit} faidx {rf}",
                       shell=True, check=True)
    else:
        open(rf, "w").close(); open(rf + ".fai", "w").close()
    reads = Fasta(rf)
    asm = Fasta(a.fasta)
    print(f"# {len(cand)} joins with >= {a.minReads} uniquely anchored spanning reads ({len(need)} reads)", file=sys.stderr)

    def build(g):
        '''-> (status, note, cutL, cutR, fill, n_agree, ident_median, n_support) in window coordinates'''
        sc, _, _, g0, g1 = g
        d = tempfile.mkdtemp(dir=tmp)
        Lw = asm.get(sc, max(0, g0 - a.window), g0); Rw = asm.get(sc, g1, g1 + a.window)
        # each flank on its own: where the flanks overlap, a joint reference would give the
        # overlapping read stretch to only one of them
        refL = os.path.join(d, "L.fa"); open(refL, "w").write(f">L\n{Lw}\n")
        refR = os.path.join(d, "R.fa"); open(refR, "w").write(f">R\n{Rw}\n")
        qf = os.path.join(d, "reads.fa")
        with open(qf, "w") as fh:
            for q, (strand, _) in span[g].items():
                s = reads.get(q); fh.write(f">{q}\n{s if strand == '+' else rc(s)}\n")
        rseq = { q: (reads.get(q) if st == "+" else rc(reads.get(q))) for q, (st, _) in span[g].items() }
        la = collections.defaultdict(dict)
        for x in run_minimap2(a.minimap2, a.preset, refL, qf) + run_minimap2(a.minimap2, a.preset, refR, qf):
            if x.strand != "+": continue
            cur = la[x.q].get(x.t)
            if x.t == "L" and (cur is None or x.te > cur.te): la[x.q]["L"] = x
            if x.t == "R" and (cur is None or x.ts < cur.ts): la[x.q]["R"] = x
        ok = { q: v for q, v in la.items() if "L" in v and "R" in v and v["L"].te >= len(Lw) - a.maxTrim
               and v["R"].ts <= a.maxTrim and v["L"].qs < v["R"].qe }
        if len(ok) < a.minReads:
            both = [ v for v in la.values() if "L" in v and "R" in v ]
            obs = (f"; the reads stop matching the left flank a median {statistics.median(len(Lw) - v['L'].te for v in both):.0f} bp "
                   f"before the gap and start on the right flank {statistics.median(v['R'].ts for v in both):.0f} bp after it") if both else ""
            return ("local_realign_failed", f"{len(ok)} reads realign to both flank ends within --maxTrim{obs}", 0, 0, "", 0, 0, 0)
        cutL = min(v["L"].te for v in ok.values()); cutR = max(v["R"].ts for v in ok.values())
        qa = { q: q_at(v["L"], cutL) for q, v in ok.items() }
        qb = { q: q_at(v["R"], cutR) for q, v in ok.items() }
        over = [ q for q in ok if qb[q] < qa[q] ]
        note = ""
        if over:
            if len(over) < len(ok) * a.minAgree:
                return ("inconsistent", f"{len(over)} of {len(ok)} reads overlap the flanks", cutL, cutR, "", 0, 0, 0)
            ts_ = sorted(t_at(ok[q]["R"], qa[q]) for q in over)
            cutR = ts_[len(ts_) // 2] # median position in R that follows the left cut
            agree = [ q for q in over if abs(t_at(ok[q]["R"], qa[q]) - cutR) <= 50 ]
            fills = { q: "" for q in agree }; med = 1.0
            note = f"flanks overlap: right flank trimmed by {cutR} bp"
            best = None
        else:
            fills = { q: rseq[q][qa[q]:qb[q]] for q in ok }
            qs_ = list(fills)
            M = { (i, j): 1.0 for i in qs_ for j in qs_ if i == j }
            for i in range(len(qs_)):
                for j in range(i + 1, len(qs_)):
                    M[(qs_[i], qs_[j])] = M[(qs_[j], qs_[i])] = identity(fills[qs_[i]], fills[qs_[j]], a.minimap2, d)
            best = max(qs_, key=lambda i: (sum(M[(i, j)] for j in qs_), -abs(len(fills[i]) - statistics.median(len(f) for f in fills.values()))))
            agree = [ j for j in qs_ if M[(best, j)] >= a.minIdentity ]
            others = [ M[(best, j)] for j in qs_ if j != best ]
            med = statistics.median(others) if others else 1.0
        if len(agree) < a.minReads or len(agree) < len(ok) * a.minAgree:
            return ("inconsistent", f"{len(agree)} of {len(ok)} reads agree with the medoid fill", cutL, cutR, "", len(agree), med, 0)
        fill = fills[best] if best else ""
        # check: the spanning reads realign through the new local sequence
        new = Lw[:cutL] + fill + Rw[cutR:]
        nref = os.path.join(d, "new.fa"); open(nref, "w").write(f">N\n{new}\n")
        t0, t1 = cutL - 1000, cutL + len(fill) + 1000
        sup = set()
        for x in run_minimap2(a.minimap2, a.preset, nref, qf):
            if x.strand == "+" and x.ts <= t0 and x.te >= t1 and max_indel(x, cutL - 500, cutL + len(fill) + 500) <= 50:
                sup.add(x.q)
        if len(sup) < a.minReads:
            return ("check_failed", f"only {len(sup)} reads align through the filled sequence", cutL, cutR, fill, len(agree), med, len(sup))
        return ("filled", note, cutL, cutR, fill, len(agree), med, len(sup))

    results = {}
    with concurrent.futures.ThreadPoolExecutor(a.threads) as ex:
        futs = { ex.submit(build, g): g for g in cand }
        for f in concurrent.futures.as_completed(futs):
            g = futs[f]
            try: results[g] = f.result()
            except subprocess.CalledProcessError as e: results[g] = ("error", str(e)[:200], 0, 0, "", 0, 0, 0)

    # graph outcome of each gap (from gapfill.py), keyed by its flanking unitigs
    graph = {}
    if a.gapfillReport:
        for l in open(a.gapfillReport):
            if l.startswith("#"): continue
            s = l.rstrip("\n").split("\t")
            graph[(s[1][:-1], s[2][:-1])] = s[3]

    # new AGP
    out_agp = open(a.out + ".agp", "w"); out_agp.write("##agp-version 2.1\n")
    out_fa = open(a.out + ".fills.fasta", "w")
    exc = open(a.out + ".excised.tsv", "w")
    exc.write("#join_id\tscaffold\tstart\tend\tcomponent\tcomp_start\tcomp_end\tstrand\tskip_reads\tdup_cover\n")
    memb = open(a.out + ".membership.tsv", "w"); memb.write("#contig\tscaffold\thow\n")
    rep = open(a.out + ".report.tsv", "w")
    rep.write("#gap_id\tleft\tright\tgraph_outcome\tmethod\tfilled\tspan_unique\tspan_lowmapq\telsewhere\tdeadend\t"
              "agree\tidentity_median\tsupport_after\ttrim_left\ttrim_right\tfill_bp\tnote\n")
    summary = collections.Counter(); added = trimmed = 0; excised_bp = [0]
    for sc in sorted(rows, key=lambda s: [ int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s) ]):
        rs = [ list(r) for r in rows[sc] ]
        edits = {} # gap row -> fill component row (or None)
        nfill = 0
        for g in sorted((J for J in joins if J[0] == sc), key=lambda J: J[1]):
            _, i, i2, g0, g1 = g
            left, right = rs[i-1], rs[i2+1]
            if g in excise:
                for k in range(i, i2 + 1):
                    edits[k] = None
                    r = rs[k]
                    if r[4] not in ("N", "U"):
                        exc.write(f"{gid[(sc, i)]}\t{sc}\t{r[1]}\t{r[2]}\t{r[5]}\t{r[6]}\t{r[7]}\t{r[8]}\t{len(span[g])}\t{dup[gid[(sc, i)]]:.3f}\n")
                    elif k != i:
                        rep.write(f"{gid[(sc, k)]}\t.\t.\t.\tmerged\t.\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\tmerged into {gid[(sc, i)]} by excision\n")
                edits[i] = rs[i] # a gap unless filled below
                xs = f"excised {sc}:{g[3]+1}-{g[4]} ({g[4]-g[3]} bp incl. gaps, {dup[gid[(sc, i)]]:.0%} duplicated elsewhere)"
                summary["excised"] += 1; excised_bp[0] += g[4] - g[3]
            lu, ru = left[5] + left[8], right[5] + right[8]
            go = graph.get((left[5], right[5]), ".")
            ns, nl = len(span[g]), len(lowq[g])
            ew = elsewhere[g]; ewn = sum(ew.values())
            res = results.get(g)
            if res and res[0] == "filled":
                st, note, cutL, cutR, fill, agree, med, sup = res
                Lw = min(a.window, g0)
                tl, tr = Lw - cutL, cutR
                lL, lR = int(left[7]) - int(left[6]) + 1, int(right[7]) - int(right[6]) + 1
                if tl >= lL or tr >= lR:
                    res = ("trim_too_long", "trim would remove a whole component", *res[2:])
            if res and res[0] == "filled":
                if tl:
                    if left[8] == "+": left[7] = str(int(left[7]) - tl)
                    else: left[6] = str(int(left[6]) + tl)
                if tr:
                    if right[8] == "+": right[6] = str(int(right[6]) + tr)
                    else: right[7] = str(int(right[7]) - tr)
                if int(left[7]) < int(left[6]) or int(right[7]) < int(right[6]):
                    sys.exit(f"ERROR: trimming emptied a component at {gid[(sc, i)]}")
                if fill:
                    nfill += 1; name = f"rf_{sc}_{nfill}"
                    out_fa.write(f">{name}\n{fill}\n")
                    edits[i] = [sc, "0", "0", "0", "W", name, "1", str(len(fill)), "+"]
                    memb.write(f"{name}\t{sc}\treadfill\n")
                else:
                    edits[i] = None
                method = ("reads_overlap" if not fill else "reads") + ("_excised" if g in excise else "")
                if g in excise: note = xs + ("; " + note if note else "")
                summary[method] += 1; added += len(fill); trimmed += tl + tr
                rep.write(f"{gid[(sc, i)]}\t{lu}\t{ru}\t{go}\t{method}\tyes\t{ns}\t{nl}\t{ewn}\t{deadend[g]}\t{agree}\t{med:.4f}\t{sup}\t{tl}\t{tr}\t{len(fill)}\t{note or '.'}\n")
                continue
            if res: method, note = res[0], res[1]
            if g in excise: method, note = "excised_gap", xs + "; fill " + (f"{method}: {note}" if res else "not attempted")
            elif res: pass
            elif ew and ew.most_common(1)[0][1] >= max(a.minReads, ns): method, note = "read_conflict", ",".join(f"{k}:{n}" for k, n in ew.most_common(3))
            elif ns: method, note = "too_few_reads", "."
            else: method, note = "no_span", "."
            if method != "read_conflict" and ewn: note = (note if note != "." else "") + (";" if note != "." else "") + "elsewhere=" + ",".join(f"{k}:{n}" for k, n in ew.most_common(3))
            summary[method] += 1
            agree = res[5] if res else 0; med = res[6] if res else 0; sup = res[7] if res else 0
            rep.write(f"{gid[(sc, i)]}\t{lu}\t{ru}\t{go}\t{method}\tno\t{ns}\t{nl}\t{ewn}\t{deadend[g]}\t{agree}\t{med:.4f}\t{sup}\t0\t0\t0\t{note}\n")
        # write the scaffold with new coordinates
        pos, part = 0, 0
        for i, r in enumerate(rs):
            if i in edits:
                if edits[i] is None: continue
                r = edits[i]
            n = int(r[7]) - int(r[6]) + 1 if r[4] not in ("N", "U") else int(r[5])
            part += 1
            r = [sc, str(pos + 1), str(pos + n), str(part)] + r[4:]
            pos += n
            out_agp.write("\t".join(r) + "\n")
    print("# " + ", ".join(f"{k}: {n}" for k, n in sorted(summary.items())), file=sys.stderr)
    nf = sum(n for k, n in summary.items() if k.startswith("reads"))
    print(f"# filled {nf} joins ({summary['excised']} after excising a duplicated segment, {excised_bp[0]} bp); "
          f"{added} bp of read sequence added, {trimmed} bp of flank trimmed", file=sys.stderr)
    with open(a.out + ".links.tsv", "w") as fh:
        fh.write("#gap_id\tedge\tread\ttarget\tjunction\torientation\tmapq\n")
        for r in sorted(links): fh.write("\t".join(map(str, r)) + "\n")
    subprocess.run(["rm", "-rf", tmp])

if __name__ == "__main__":
    main()
