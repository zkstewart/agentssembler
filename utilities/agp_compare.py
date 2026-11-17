#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
agp_compare.py -- compare two AGP scaffoldings built from the same set of
components (contigs / unitigs), e.g. agentssembler against HapHiC or YaHS.

The two AGPs are referred to as A (the "reference" naming, e.g. chr1_h1) and
B (the scaffolding to be matched onto A, e.g. HapHiC group1).  Neither is
treated as the truth; A only supplies the names and the x axis.  No reads,
FASTA or reference genome are needed.

Method
------
1. Markers.  Every component is cut at every component boundary used by either
   AGP, giving "atomic" segments (component, start, end).  A segment used by
   both AGPs is a shared marker, weighted by its length.  This copes with one
   tool using a whole unitig while the other uses a trimmed range or splits it.
2. Scaffold matching.  A shared-bp contingency table between A scaffolds and
   B scaffolds gives each B scaffold its best A scaffold (and vice versa), a
   one-to-one optimal assignment, and purity / completeness figures.
3. Adjacency.  Each AGP's scaffolds are projected onto the markers that are
   scaffolded in BOTH AGPs.  Consecutive markers form junctions; a junction is
   conserved when the same two component ends are adjacent in the other AGP.
4. Blocks.  Runs of conserved junctions form collinear blocks.  Per matched
   scaffold pair, a weighted longest-increasing-subsequence over the blocks
   gives the collinear backbone, and each block is classed as collinear,
   inverted (in place, opposite orientation), relocated (same scaffold, out of
   order) or other_scaffold (placed in a B scaffold that belongs to a
   different A scaffold).
5. Gaps.  Each gap is located in the other AGP via the junction it sits in:
   is the join also made there, and with or without a gap?

Outputs (in --outdir, prefixed with --prefix)
---------------------------------------------
  summary.txt          human-readable report (also printed)
  matches_B_to_A.tsv   one row per B scaffold: best A scaffold and pair metrics
  matches_A_to_B.tsv   one row per A scaffold: how B covers it
  contingency.tsv      bp shared between every A scaffold and B scaffold
  markers.tsv          every shared marker with its placement in both AGPs
  blocks.tsv           collinear blocks with their class
  junctions.tsv        every projected junction of A and of B, and its fate
  gaps.tsv             every gap of A and of B, located in the other AGP
  dotplot.png          genome-wide block dot plot (needs matplotlib)

Example
-------
  agp_compare.py agentssembler.agp haphic.agp \\
      --labels agentssembler haphic -o out --prefix 13Q010
'''

import argparse
import bisect
import os
import re
import sys
from collections import Counter, defaultdict

COMPONENT_TYPES = set("ADFGOPW")
GAP_TYPES = set("NU")
UNSCAFFOLDED = "(unscaffolded)"


# --------------------------------------------------------------------------
# AGP parsing
# --------------------------------------------------------------------------
class Comp:
    __slots__ = ("obj", "obj_beg", "obj_end", "cid", "beg", "end", "strand")

    def __init__(self, obj, obj_beg, obj_end, cid, beg, end, strand):
        self.obj, self.obj_beg, self.obj_end = obj, obj_beg, obj_end
        self.cid, self.beg, self.end, self.strand = cid, beg, end, strand

    @property
    def length(self):
        return self.end - self.beg + 1

    def describe(self):
        return f"{self.cid}:{self.beg}-{self.end}({'+' if self.strand > 0 else '-'})"


class Gap:
    __slots__ = ("obj", "obj_beg", "obj_end", "length", "kind", "gap_type",
                 "linkage", "evidence")

    def __init__(self, obj, obj_beg, obj_end, length, kind, gap_type, linkage,
                 evidence):
        self.obj, self.obj_beg, self.obj_end = obj, obj_beg, obj_end
        self.length, self.kind, self.gap_type = length, kind, gap_type
        self.linkage, self.evidence = linkage, evidence


class Agp:
    def __init__(self, path, label, min_components=2, scaffold_regex=None):
        self.path, self.label = path, label
        self.entries = defaultdict(list)     # object -> [Comp | Gap] in order
        self.warnings = []
        self.unknown_orient = 0
        self._parse()
        self.obj_len = {o: max(e.obj_end for e in es)
                        for o, es in self.entries.items()}
        self.comps = {o: [e for e in es if isinstance(e, Comp)]
                      for o, es in self.entries.items()}
        self.gaps = {o: [e for e in es if isinstance(e, Gap)]
                     for o, es in self.entries.items()}
        if scaffold_regex:
            rx = re.compile(scaffold_regex)
            self.scaffolds = {o for o in self.entries if rx.search(o)}
        else:
            self.scaffolds = {o for o, cs in self.comps.items()
                              if len(cs) >= min_components}

    def _parse(self):
        with open(self.path) as fh:
            for n, line in enumerate(fh, 1):
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                if len(f) < 8:
                    sys.exit(f"{self.path}:{n}: expected >=8 tab-separated "
                             f"columns, found {len(f)}")
                obj, obeg, oend, ctype = f[0], int(f[1]), int(f[2]), f[4]
                if ctype in GAP_TYPES:
                    e = Gap(obj, obeg, oend, int(f[5]), ctype, f[6], f[7],
                            f[8] if len(f) > 8 else "na")
                elif ctype in COMPONENT_TYPES:
                    orient = f[8] if len(f) > 8 else "+"
                    if orient not in "+-" or orient == "":
                        self.unknown_orient += 1
                    e = Comp(obj, obeg, oend, f[5], int(f[6]), int(f[7]),
                             -1 if orient == "-" else 1)
                    if e.length != oend - obeg + 1:
                        self.warnings.append(
                            f"line {n}: object span != component span")
                else:
                    sys.exit(f"{self.path}:{n}: unknown component type {ctype!r}")
                self.entries[obj].append(e)
        for obj, es in self.entries.items():
            es.sort(key=lambda e: e.obj_beg)
            pos = 0
            for e in es:
                if e.obj_beg != pos + 1:
                    self.warnings.append(
                        f"{obj}: coordinates not contiguous at {e.obj_beg}")
                pos = e.obj_end

    def is_scaffold(self, obj):
        return obj in self.scaffolds


class Place:
    """Placement of one marker in one AGP."""
    __slots__ = ("obj", "beg", "end", "strand", "comp")

    def __init__(self, obj, beg, end, strand, comp):
        self.obj, self.beg, self.end = obj, beg, end
        self.strand, self.comp = strand, comp


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def natural_key(s):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


def fmt_bp(n):
    return f"{n:,}"


def pct(a, b):
    return f"{100.0 * a / b:.2f}%" if b else "n/a"


def frac(a, b):
    return f"{a / b:.4f}" if b else "NA"


def key_str(key):
    return f"{key[0]}:{key[1]}-{key[2]}"


def key_len(key):
    return key[2] - key[1] + 1


def n50(lengths):
    total, run = sum(lengths), 0
    for x in sorted(lengths, reverse=True):
        run += x
        if run * 2 >= total:
            return x
    return 0


def weighted_lis(values, weights):
    """Max-weight strictly increasing subsequence. Returns (weight, indices)."""
    n = len(values)
    if n == 0:
        return 0, []
    best, prev = list(weights), [-1] * n
    for i in range(n):
        for j in range(i):
            if values[j] < values[i] and best[j] + weights[i] > best[i]:
                best[i], prev[i] = best[j] + weights[i], j
    i = max(range(n), key=best.__getitem__)
    total, idx = best[i], []
    while i >= 0:
        idx.append(i)
        i = prev[i]
    return total, idx[::-1]


def optimal_assignment(rows, cols, W):
    """One-to-one assignment maximising shared bp. Returns ({row: col}, method)."""
    pairs = {}
    try:
        import numpy as np
        from scipy.optimize import linear_sum_assignment
        M = np.array([[W[r].get(c, 0) for c in cols] for r in rows], dtype=float)
        if M.size:
            ri, ci = linear_sum_assignment(M, maximize=True)
            pairs = {rows[i]: cols[j] for i, j in zip(ri, ci) if M[i, j] > 0}
        return pairs, "Hungarian (scipy)"
    except ImportError:
        cand = sorted(((W[r].get(c, 0), r, c) for r in rows for c in cols),
                      reverse=True)
        used = set()
        for w, r, c in cand:
            if w > 0 and r not in pairs and c not in used:
                pairs[r] = c
                used.add(c)
        return pairs, "greedy (scipy not available)"


# --------------------------------------------------------------------------
# core steps
# --------------------------------------------------------------------------
def atomise(agps):
    """Cut components at all boundaries used by any AGP; place each piece."""
    cuts = defaultdict(set)
    for agp in agps:
        for cs in agp.comps.values():
            for c in cs:
                cuts[c.cid].update((c.beg, c.end + 1))
    cuts = {cid: sorted(s) for cid, s in cuts.items()}
    placements = []
    for agp in agps:
        pl, dup = {}, 0
        for cs in agp.comps.values():
            for c in cs:
                pts = cuts[c.cid]
                i = bisect.bisect_left(pts, c.beg)
                while pts[i] <= c.end:
                    s, e = pts[i], pts[i + 1] - 1
                    key = (c.cid, s, e)
                    if key in pl:
                        dup += e - s + 1
                    else:
                        ob = c.obj_beg + (s - c.beg if c.strand > 0 else c.end - e)
                        pl[key] = Place(c.obj, ob, ob + e - s, c.strand, c)
                    i += 1
        if dup:
            agp.warnings.append(
                f"{fmt_bp(dup)} bp of component sequence is used more than "
                f"once; only the first use is compared")
        placements.append(pl)
    return placements


def project(agp, placement, markers):
    """Per scaffold, the shared-scaffolded markers in object order."""
    proj = defaultdict(list)
    for key in markers:
        proj[placement[key].obj].append(key)
    for obj in proj:
        proj[obj].sort(key=lambda k: placement[k].beg)
    return proj


def end_right(key, place):
    return (key, "h" if place.strand > 0 else "t")


def end_left(key, place):
    return (key, "t" if place.strand > 0 else "h")


def build_junctions(src, agp, placement, proj):
    """Junctions between consecutive projected markers of one AGP."""
    juncs = []
    for obj, keys in proj.items():
        gaps, comps = agp.gaps[obj], agp.comps[obj]
        for i in range(len(keys) - 1):
            kl, kr = keys[i], keys[i + 1]
            pl, pr = placement[kl], placement[kr]
            between_gaps = [g for g in gaps
                            if g.obj_beg > pl.end and g.obj_end < pr.beg]
            between_comps = [c for c in comps
                             if c.obj_beg > pl.end and c.obj_end < pr.beg]
            gap_bp = sum(g.length for g in between_gaps)
            juncs.append({
                "src": src, "obj": obj, "left": kl, "right": kr,
                "oriented": frozenset((end_right(kl, pl), end_left(kr, pr))),
                "unoriented": frozenset((kl, kr)),
                "gaps": between_gaps, "n_comps": len(between_comps),
                "seq_bp": pr.beg - pl.end - 1 - gap_bp,
                "status": None, "other": None,
            })
    return juncs


def classify_junctions(juncs, other_juncs, other_placement):
    """Label each junction by what the other AGP does with the same two ends."""
    by_oriented = {j["oriented"]: j for j in other_juncs}
    by_unoriented = {j["unoriented"]: j for j in other_juncs}
    for j in juncs:
        if j["oriented"] in by_oriented:
            j["status"], j["other"] = "adjacent", by_oriented[j["oriented"]]
        elif j["unoriented"] in by_unoriented:
            j["status"] = "adjacent_orientation_differs"
            j["other"] = by_unoriented[j["unoriented"]]
        elif other_placement[j["left"]].obj == other_placement[j["right"]].obj:
            j["status"] = "same_scaffold_not_adjacent"
        else:
            j["status"] = "different_scaffolds"


def build_blocks(projA, juncsA, plA, plB):
    """Maximal runs of markers whose junctions are conserved in B."""
    jmap = {(j["left"], j["right"]): j for j in juncsA}
    blocks = []
    for obj in sorted(projA, key=natural_key):
        keys = projA[obj]
        run = [keys[0]]
        for prev, cur in zip(keys, keys[1:]):
            if jmap[(prev, cur)]["status"] == "adjacent":
                run.append(cur)
            else:
                blocks.append(run)
                run = [cur]
        blocks.append(run)
    out = []
    for run in blocks:
        a0, b0 = plA[run[0]], plB[run[0]]
        out.append({
            "a_obj": a0.obj,
            "a_beg": min(plA[k].beg for k in run),
            "a_end": max(plA[k].end for k in run),
            "b_obj": b0.obj,
            "b_beg": min(plB[k].beg for k in run),
            "b_end": max(plB[k].end for k in run),
            "rel": a0.strand * b0.strand,       # +1 same strand, -1 opposite
            "n": len(run), "bp": sum(key_len(k) for k in run),
            "keys": run, "cls": None, "detail": "",
        })
    return out


# --------------------------------------------------------------------------
# main analysis
# --------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("agp_a", help="AGP A (supplies scaffold names, e.g. chr1_h1)")
    ap.add_argument("agp_b", help="AGP B (scaffolds are matched onto A)")
    ap.add_argument("-o", "--outdir", default=".", help="output directory [.]")
    ap.add_argument("-p", "--prefix", default="agp_compare",
                    help="output file prefix [agp_compare]")
    ap.add_argument("--labels", nargs=2, default=["A", "B"], metavar=("A", "B"),
                    help="short names for the two AGPs [A B]")
    ap.add_argument("--min_components", type=int, default=2,
                    help="an object with at least this many components is a "
                         "scaffold; smaller objects count as unscaffolded [2]")
    ap.add_argument("--scaffold_regex_a", help="instead of --min_components, "
                    "objects of A matching this regex are scaffolds")
    ap.add_argument("--scaffold_regex_b", help="as above, for B")
    ap.add_argument("--chrom_regex", default=r"^(.+)_h(?:ap)?\d+[A-Za-z]?$",
                    help="regex whose group 1 is the chromosome of an A "
                         "scaffold name, used to tell haplotype swaps from "
                         "inter-chromosome joins [%(default)s]")
    ap.add_argument("--no_plot", action="store_true", help="skip the dot plot")
    args = ap.parse_args()

    la, lb = args.labels
    A = Agp(args.agp_a, la, args.min_components, args.scaffold_regex_a)
    B = Agp(args.agp_b, lb, args.min_components, args.scaffold_regex_b)
    plA, plB = atomise([A, B])
    os.makedirs(args.outdir, exist_ok=True)
    out = lambda name: os.path.join(args.outdir, f"{args.prefix}.{name}")

    lines = []

    def say(text=""):
        lines.append(text)

    # ---- input overview --------------------------------------------------
    say("AGP comparison")
    say("==============")
    for tag, X in (("A", A), ("B", B)):
        ncomp = sum(len(X.comps[o]) for o in X.scaffolds)
        cbp = sum(c.length for o in X.scaffolds for c in X.comps[o])
        ngap = sum(len(X.gaps[o]) for o in X.scaffolds)
        uns = [o for o in X.entries if o not in X.scaffolds]
        ubp = sum(c.length for o in uns for c in X.comps[o])
        say(f"{tag} = {X.label}: {X.path}")
        say(f"    scaffolds {len(X.scaffolds)}, components in scaffolds {ncomp} "
            f"({fmt_bp(cbp)} bp), gaps {ngap}")
        say(f"    unscaffolded objects {len(uns)} ({fmt_bp(ubp)} bp)")
        if X.unknown_orient:
            say(f"    WARNING: {X.unknown_orient} components have unknown "
                f"orientation; treated as '+'")
        for w in X.warnings[:10]:
            say(f"    WARNING: {w}")
        if len(X.warnings) > 10:
            say(f"    ... and {len(X.warnings) - 10} more warnings")
    say()

    # ---- component membership -------------------------------------------
    idsA = {c.cid for cs in A.comps.values() for c in cs}
    idsB = {c.cid for cs in B.comps.values() for c in cs}

    def where(key, X, pl):
        if key not in pl:
            ids = idsA if X is A else idsB
            return "part not used" if key[0] in ids else "absent"
        return "scaffolded" if X.is_scaffold(pl[key].obj) else "unscaffolded"

    say("1. Component content")
    say("--------------------")
    say(f"Component IDs: {len(idsA & idsB)} in both AGPs, {len(idsA - idsB)} "
        f"only in A, {len(idsB - idsA)} only in B")
    for only, name in ((idsA - idsB, "A"), (idsB - idsA, "B")):
        if 0 < len(only) <= 10:
            say(f"    only in {name}: {', '.join(sorted(only, key=natural_key))}")
    for tag, X, pl, Y, plo in (("A", A, plA, B, plB), ("B", B, plB, A, plA)):
        tally = Counter()
        for key, p in pl.items():
            if X.is_scaffold(p.obj):
                tally[where(key, Y, plo)] += key_len(key)
        tot = sum(tally.values())
        other = "B" if tag == "A" else "A"
        say(f"Sequence scaffolded in {tag} ({fmt_bp(tot)} bp) is, in {other}:")
        for k in ("scaffolded", "unscaffolded", "part not used", "absent"):
            say(f"    {k:<14}{fmt_bp(tally[k]):>16} bp  {pct(tally[k], tot):>8}")
    say("  part not used = the component is in the other AGP, but not this part of it")
    say("  (e.g. a trimmed overlap); absent = the whole component is missing there.")
    markers = [k for k in plA if k in plB
               and A.is_scaffold(plA[k].obj) and B.is_scaffold(plB[k].obj)]
    total_bp = sum(key_len(k) for k in markers)
    say(f"Shared markers scaffolded in both: {len(markers)} "
        f"({fmt_bp(total_bp)} bp). All metrics below use these.")
    say()
    if not markers:
        say("No sequence is scaffolded in both AGPs; nothing to compare.")
        print("\n".join(lines))
        return

    # ---- contingency & matching -----------------------------------------
    W = defaultdict(lambda: defaultdict(int))    # W[a][b] shared bp
    Wt = defaultdict(lambda: defaultdict(int))   # transposed
    full = defaultdict(lambda: defaultdict(int))  # incl. unscaffolded/absent
    for key in set(plA) | set(plB):
        a = plA[key].obj if key in plA and A.is_scaffold(plA[key].obj) else UNSCAFFOLDED
        b = plB[key].obj if key in plB and B.is_scaffold(plB[key].obj) else UNSCAFFOLDED
        full[a][b] += key_len(key)
        if a != UNSCAFFOLDED and b != UNSCAFFOLDED:
            W[a][b] += key_len(key)
            Wt[b][a] += key_len(key)
    a_names = sorted(A.scaffolds, key=natural_key)
    b_names = sorted(B.scaffolds, key=natural_key)
    best_a = {b: max(Wt[b], key=lambda a: (Wt[b][a], a)) for b in b_names if Wt[b]}
    best_b = {a: max(W[a], key=lambda b: (W[a][b], b)) for a in a_names if W[a]}
    assign, assign_method = optimal_assignment(a_names, b_names, W)
    assign_rev = {b: a for a, b in assign.items()}

    rx = re.compile(args.chrom_regex)
    chrom = {}
    for a in a_names:
        m = rx.match(a)
        if m:
            chrom[a] = m.group(1)
    use_chrom = (len(chrom) == len(a_names)
                 and 1 < len(set(chrom.values())) < len(a_names))

    with open(out("contingency.tsv"), "w") as fh:
        cols = b_names + [UNSCAFFOLDED]
        fh.write(f"{la}\\{lb}\t" + "\t".join(cols) + "\n")
        for a in a_names + [UNSCAFFOLDED]:
            fh.write(a + "\t" + "\t".join(str(full[a].get(b, 0)) for b in cols) + "\n")

    say("2. Scaffold matching (shared bp)")
    say("--------------------------------")
    recip = [b for b in b_names if b in best_a and best_b.get(best_a[b]) == b]
    pur_b = sum(max(Wt[b].values()) for b in b_names if Wt[b])
    pur_a = sum(max(W[a].values()) for a in a_names if W[a])
    asg_bp = sum(W[a][b] for a, b in assign.items())
    say(f"B scaffolds with a best A scaffold: {len(best_a)} of {len(b_names)}; "
        f"reciprocal best pairs: {len(recip)}")
    say(f"A scaffolds that are no B scaffold's best match: "
        f"{len(set(a_names) - set(best_a.values()))}")
    say(f"Purity of B scaffolds (bp in their best A scaffold):   {pct(pur_b, total_bp)}")
    say(f"Purity of A scaffolds (bp in their best B scaffold):   {pct(pur_a, total_bp)}")
    say(f"bp on the one-to-one assignment [{assign_method}]: {pct(asg_bp, total_bp)}")
    if use_chrom:
        cw = defaultdict(lambda: defaultdict(int))
        for b in b_names:
            for a, w in Wt[b].items():
                cw[b][chrom[a]] += w
        cpur = sum(max(cw[b].values()) for b in cw)
        mixed_chr = [b for b in cw if len([c for c, w in cw[b].items() if w]) > 1]
        mixed_hap = [b for b in b_names if len(Wt[b]) > 1]
        say(f"Chromosome-level purity of B scaffolds:                {pct(cpur, total_bp)}")
        say(f"B scaffolds mixing >1 A scaffold: {len(mixed_hap)}; "
            f"mixing >1 chromosome: {len(mixed_chr)}")
    say()

    # ---- projection, junctions, blocks ----------------------------------
    projA, projB = project(A, plA, markers), project(B, plB, markers)
    jA = build_junctions("A", A, plA, projA)
    jB = build_junctions("B", B, plB, projB)
    classify_junctions(jA, jB, plB)
    classify_junctions(jB, jA, plA)
    blocks = build_blocks(projA, jA, plA, plB)

    # per-pair orientation and collinear backbone (pairs = B scaffold -> best A)
    flip, pair_stats = {}, {}
    for b in b_names:
        if b not in best_a:
            continue
        a = best_a[b]
        pb = [blk for blk in blocks if blk["a_obj"] == a and blk["b_obj"] == b]
        wts = [blk["bp"] for blk in pb]
        fw, fidx = weighted_lis([blk["b_beg"] for blk in pb], wts)
        rw, ridx = weighted_lis([-blk["b_beg"] for blk in pb], wts)
        vote = sum(blk["rel"] * blk["bp"] for blk in pb)
        f = 1 if fw > rw or (fw == rw and vote >= 0) else -1
        flip[b] = f
        backbone = set(fidx if f > 0 else ridx)
        for i, blk in enumerate(pb):
            ok = blk["rel"] == f
            if i in backbone:
                blk["cls"] = "collinear" if ok else "inverted"
            else:
                blk["cls"] = "relocated"
                blk["detail"] = "orientation_kept" if ok else "orientation_flipped"
        pair_stats[b] = {
            "n_markers": sum(blk["n"] for blk in pb), "bp": sum(wts),
            "n_blocks": len(pb), "collinear_bp": max(fw, rw),
            "orient_bp": sum(blk["bp"] for blk in pb if blk["rel"] == f),
            "inv": [blk for blk in pb if blk["cls"] == "inverted"],
            "rel": [blk for blk in pb if blk["cls"] == "relocated"],
        }
    for blk in blocks:
        if blk["cls"] is None:
            blk["cls"] = "other_scaffold"
            owner = best_a.get(blk["b_obj"])
            if use_chrom and owner:
                blk["detail"] = ("other_haplotype"
                                 if chrom[owner] == chrom[blk["a_obj"]]
                                 else "other_chromosome")

    # ---- adjacency summary ----------------------------------------------
    say("3. Order and orientation (junctions between shared markers)")
    say("----------------------------------------------------------")
    statuses = ("adjacent", "adjacent_orientation_differs",
                "same_scaffold_not_adjacent", "different_scaffolds")
    for tag, js, other in (("A", jA, "B"), ("B", jB, "A")):
        c = Counter(j["status"] for j in js)
        say(f"{len(js)} junctions in {tag}; in {other} the same two markers are:")
        for s in statuses:
            say(f"    {s:<30}{c[s]:>6}  {pct(c[s], len(js)):>8}")
    nA, nB = len(jA), len(jB)
    shared_or = sum(1 for j in jA if j["status"] == "adjacent")
    shared_un = shared_or + sum(1 for j in jA
                                if j["status"] == "adjacent_orientation_differs")
    say(f"Ordering agreement    (adjacent, any orientation): "
        f"{pct(shared_un, nA)} of A junctions, {pct(shared_un, nB)} of B junctions")
    say(f"Orientation agreement (of those, same ends joined): {pct(shared_or, shared_un)}")
    say(f"Breakpoints: {nA - shared_or} A junctions are not in B, "
        f"{nB - shared_or} B junctions are not in A")
    say()
    say(f"Collinear blocks: {len(blocks)} (N50 {fmt_bp(n50([b['bp'] for b in blocks]))} "
        f"bp of shared sequence)")
    order = ("collinear", "inverted", "relocated", "other_scaffold")
    for cls in order:
        sel = [b for b in blocks if b["cls"] == cls]
        bp = sum(b["bp"] for b in sel)
        extra = ""
        if cls == "other_scaffold" and use_chrom and sel:
            d = Counter()
            for b in sel:
                d[b["detail"]] += b["bp"]
            extra = "  (" + ", ".join(f"{k} {fmt_bp(v)} bp" for k, v in sorted(d.items())) + ")"
        say(f"    {cls:<15}{len(sel):>5} blocks {fmt_bp(bp):>15} bp  {pct(bp, total_bp):>8}{extra}")
    say("  collinear = on the pair's backbone, same orientation; inverted = in place")
    say("  but flipped; relocated = same scaffold pair, out of order; other_scaffold =")
    say("  in a B scaffold whose best match is a different A scaffold.")
    say()

    # ---- gaps -----------------------------------------------------------
    say("4. Gaps")
    say("-------")
    gap_rows = []
    for tag, X, plX, projX, js, Y, plY in (("A", A, plA, projA, jA, B, plB),
                                           ("B", B, plB, projB, jB, A, plA)):
        gaps = [g for o in X.scaffolds for g in X.gaps[o]]
        types = Counter((g.kind, g.gap_type, g.evidence) for g in gaps)
        say(f"{tag}: {len(gaps)} gaps, {fmt_bp(sum(g.length for g in gaps))} bp; "
            + "; ".join(f"{n} x {k}/{t}/{e}" for (k, t, e), n in types.most_common()))
        jmap = defaultdict(list)
        for j in js:
            jmap[j["obj"]].append(j)
        tally = Counter()
        for obj in sorted(X.scaffolds, key=natural_key):
            entries = X.entries[obj]
            for i, g in enumerate(entries):
                if not isinstance(g, Gap):
                    continue
                lc = next((e for e in reversed(entries[:i]) if isinstance(e, Comp)), None)
                rc = next((e for e in entries[i + 1:] if isinstance(e, Comp)), None)
                j = next((j for j in jmap[obj] if g in j["gaps"]), None)
                row = {"src": X.label, "obj": obj, "beg": g.obj_beg, "end": g.obj_end,
                       "len": g.length, "type": f"{g.kind}/{g.gap_type}/{g.evidence}",
                       "left": lc.describe() if lc else "NA",
                       "right": rc.describe() if rc else "NA"}
                if j is None:
                    row.update(status="outside_shared_markers", lm="NA", rm="NA",
                               oobj="NA", ogaps="NA", ocomps="NA", obp="NA",
                               verdict="not_assessable")
                else:
                    o = j["other"]
                    oobjs = sorted({plY[j["left"]].obj, plY[j["right"]].obj},
                                   key=natural_key)
                    row.update(status=j["status"], lm=key_str(j["left"]),
                               rm=key_str(j["right"]), oobj=",".join(oobjs),
                               ogaps=len(o["gaps"]) if o else "NA",
                               ocomps=o["n_comps"] if o else "NA",
                               obp=o["seq_bp"] if o else "NA")
                    if j["status"] == "adjacent":
                        row["verdict"] = ("gap_in_both" if o["gaps"]
                                          else "join_gapless_in_other")
                    elif j["status"] == "adjacent_orientation_differs":
                        row["verdict"] = "join_made_orientation_differs"
                    else:
                        row["verdict"] = "join_not_made_in_other"
                tally[row["verdict"]] += 1
                gap_rows.append(row)
        other = "B" if tag == "A" else "A"
        for v in ("gap_in_both", "join_gapless_in_other",
                  "join_made_orientation_differs", "join_not_made_in_other",
                  "not_assessable"):
            if tally[v]:
                say(f"    {v:<30}{tally[v]:>6}  {pct(tally[v], len(gaps)):>8}")
        filled = [r for r in gap_rows if r["src"] == X.label
                  and r["verdict"] in ("gap_in_both", "join_gapless_in_other")]
        if filled:
            ncomp = sum(r["ocomps"] for r in filled)
            nbp = sum(r["obp"] for r in filled)
            withc = sum(1 for r in filled if r["ocomps"] > 0)
            say(f"    where {other} makes the same join, it places {ncomp} extra "
                f"components at {withc} of those {len(filled)} gaps")
            say(f"    ({fmt_bp(nbp)} bp of sequence between the flanking markers in {other})")
    cons = [j for j in jA if j["status"] == "adjacent"]
    tab = Counter((bool(j["gaps"]), bool(j["other"]["gaps"])) for j in cons)
    say(f"Of the {len(cons)} junctions conserved in both (same ends joined):")
    say(f"    gap in both {tab[(True, True)]}, gap only in A {tab[(True, False)]}, "
        f"gap only in B {tab[(False, True)]}, gapless in both {tab[(False, False)]}")
    say("  join_gapless_in_other = the other AGP joins the same markers with no gap;")
    say("  join_not_made_in_other = the flanking markers are not adjacent there;")
    say("  not_assessable = the gap lies outside the outermost shared markers.")
    say()

    # ---- per-scaffold tables --------------------------------------------
    a_cbp = {a: sum(c.length for c in A.comps[a]) for a in a_names}
    b_cbp = {b: sum(c.length for c in B.comps[b]) for b in b_names}
    hdr = [f"{lb}_scaffold", f"{lb}_length", f"{lb}_components", f"best_{la}",
           "shared_bp", "purity", f"{la}_completeness", f"{la}_recall_of_scaffolded",
           "reciprocal_best", f"one_to_one_{la}", f"other_{la}_hits",
           "chromosome_purity", "relative_orientation", "markers", "blocks",
           "collinear_bp_frac", "orientation_agree_bp_frac", "inverted_blocks",
           "inverted_bp", "relocated_blocks", "relocated_bp",
           f"bp_not_scaffolded_in_{la}"]
    say(f"5. Per-scaffold matches ({lb} -> {la})")
    say("--------------------------------------")
    wb = max([len(lb)] + [len(b) for b in b_names]) + 2
    wa = max([len(la) + 5] + [len(a) for a in a_names]) + 2
    say(f"{lb:<{wb}}{'best ' + la:<{wa}}{'shared Mb':>10}{'purity':>8}{'compl.':>8}"
        f"{'ori':>4}{'blocks':>7}{'collin.':>8}{'orient.':>8}{'inv':>4}{'reloc':>6}  other hits (Mb)")
    with open(out("matches_B_to_A.tsv"), "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for b in b_names:
            if b not in best_a:
                fh.write("\t".join([b, str(B.obj_len[b]), str(len(B.comps[b]))]
                                   + ["NA"] * (len(hdr) - 4) + [str(b_cbp[b])]) + "\n")
                say(f"{b:<{wb}}{'-':<{wa}} no sequence scaffolded in {la}")
                continue
            a, tot = best_a[b], sum(Wt[b].values())
            w, ps = Wt[b][a], pair_stats[b]
            others = sorted(((v, k) for k, v in Wt[b].items() if k != a), reverse=True)
            ostr = ";".join(f"{k}:{v}" for v, k in others) or "."
            cp = frac(max(cw[b].values()), tot) if use_chrom else "NA"
            fh.write("\t".join(map(str, [
                b, B.obj_len[b], len(B.comps[b]), a, w, frac(w, tot),
                frac(w, a_cbp[a]), frac(w, sum(W[a].values())),
                "yes" if best_b.get(a) == b else "no", assign_rev.get(b, "NA"),
                ostr, cp, "+" if flip[b] > 0 else "-", ps["n_markers"],
                ps["n_blocks"], frac(ps["collinear_bp"], ps["bp"]),
                frac(ps["orient_bp"], ps["bp"]), len(ps["inv"]),
                sum(x["bp"] for x in ps["inv"]), len(ps["rel"]),
                sum(x["bp"] for x in ps["rel"]), b_cbp[b] - tot])) + "\n")
            short = ", ".join(f"{k} {v / 1e6:.1f}" for v, k in others[:3]) or "."
            say(f"{b:<{wb}}{a:<{wa}}{w / 1e6:>10.2f}{w / tot:>8.3f}{w / a_cbp[a]:>8.3f}"
                f"{'+' if flip[b] > 0 else '-':>4}{ps['n_blocks']:>7}"
                f"{ps['collinear_bp'] / ps['bp']:>8.3f}{ps['orient_bp'] / ps['bp']:>8.3f}"
                f"{len(ps['inv']):>4}{len(ps['rel']):>6}  {short}")
    say("  purity = share of the B scaffold's shared bp in its best A scaffold;")
    say("  compl. = share of the A scaffold's sequence held by this B scaffold;")
    say("  collin./orient. = bp fraction on the collinear backbone / agreeing in strand.")
    say()

    with open(out("matches_A_to_B.tsv"), "w") as fh:
        fh.write("\t".join([
            f"{la}_scaffold", f"{la}_length", f"{la}_components", f"{la}_gaps",
            f"bp_scaffolded_in_{lb}", f"bp_unscaffolded_or_absent_in_{lb}",
            f"best_{lb}", "shared_bp", "frac_of_scaffold", "reciprocal_best",
            f"one_to_one_{lb}", f"{lb}_scaffolds_hit", f"other_{lb}_hits",
            "blocks", "junctions", "junctions_conserved"]) + "\n")
        jc = Counter(j["obj"] for j in jA)
        jcons = Counter(j["obj"] for j in jA if j["status"] == "adjacent")
        nblk = Counter(b["a_obj"] for b in blocks)
        for a in a_names:
            tot = sum(W[a].values())
            b = best_b.get(a)
            others = sorted(((v, k) for k, v in W[a].items() if k != b), reverse=True)
            fh.write("\t".join(map(str, [
                a, A.obj_len[a], len(A.comps[a]), len(A.gaps[a]), tot,
                a_cbp[a] - tot, b or "NA", W[a][b] if b else 0,
                frac(W[a][b], a_cbp[a]) if b else "NA",
                "yes" if b and best_a.get(b) == a else "no",
                assign.get(a, "NA"), len(W[a]),
                ";".join(f"{k}:{v}" for v, k in others) or ".",
                nblk[a], jc[a], jcons[a]])) + "\n")

    # ---- detail tables --------------------------------------------------
    block_of = {}
    for i, blk in enumerate(blocks, 1):
        blk["id"] = i
        for k in blk["keys"]:
            block_of[k] = blk
    with open(out("markers.tsv"), "w") as fh:
        fh.write("\t".join([
            "component", "comp_beg", "comp_end", "length", f"{la}_scaffold",
            f"{la}_beg", f"{la}_end", f"{la}_strand", f"{lb}_scaffold",
            f"{lb}_beg", f"{lb}_end", f"{lb}_strand", "block", "class"]) + "\n")
        for a in a_names:
            for k in projA.get(a, []):
                pa, pb, blk = plA[k], plB[k], block_of[k]
                fh.write("\t".join(map(str, [
                    k[0], k[1], k[2], key_len(k), pa.obj, pa.beg, pa.end,
                    "+" if pa.strand > 0 else "-", pb.obj, pb.beg, pb.end,
                    "+" if pb.strand > 0 else "-", blk["id"], blk["cls"]])) + "\n")
    with open(out("blocks.tsv"), "w") as fh:
        fh.write("\t".join([
            "block", f"{la}_scaffold", f"{la}_beg", f"{la}_end", f"{lb}_scaffold",
            f"{lb}_beg", f"{lb}_end", "relative_strand", "markers", "shared_bp",
            "class", "detail", f"{lb}_scaffold_best_{la}", "first_marker",
            "last_marker"]) + "\n")
        for blk in blocks:
            fh.write("\t".join(map(str, [
                blk["id"], blk["a_obj"], blk["a_beg"], blk["a_end"], blk["b_obj"],
                blk["b_beg"], blk["b_end"], "+" if blk["rel"] > 0 else "-",
                blk["n"], blk["bp"], blk["cls"], blk["detail"] or ".",
                best_a.get(blk["b_obj"], "NA"), key_str(blk["keys"][0]),
                key_str(blk["keys"][-1])])) + "\n")
    with open(out("junctions.tsv"), "w") as fh:
        fh.write("\t".join([
            "source", "scaffold", "left_marker", "right_marker", "gaps",
            "components_between", "sequence_bp_between", "status_in_other",
            "other_scaffold_left", "other_scaffold_right", "other_gaps",
            "other_components_between", "other_sequence_bp_between"]) + "\n")
        for js, X, plY in ((jA, A, plB), (jB, B, plA)):
            for j in js:
                o = j["other"]
                fh.write("\t".join(map(str, [
                    X.label, j["obj"], key_str(j["left"]), key_str(j["right"]),
                    len(j["gaps"]), j["n_comps"], j["seq_bp"], j["status"],
                    plY[j["left"]].obj, plY[j["right"]].obj,
                    len(o["gaps"]) if o else "NA", o["n_comps"] if o else "NA",
                    o["seq_bp"] if o else "NA"])) + "\n")
    with open(out("gaps.tsv"), "w") as fh:
        fh.write("\t".join([
            "source", "scaffold", "gap_beg", "gap_end", "gap_length", "gap_type",
            "left_component", "right_component", "left_shared_marker",
            "right_shared_marker", "junction_status_in_other", "verdict",
            "other_scaffold", "other_gaps", "other_components_between",
            "other_sequence_bp_between"]) + "\n")
        for r in gap_rows:
            fh.write("\t".join(map(str, [
                r["src"], r["obj"], r["beg"], r["end"], r["len"], r["type"],
                r["left"], r["right"], r["lm"], r["rm"], r["status"],
                r["verdict"], r["oobj"], r["ogaps"], r["ocomps"], r["obp"]])) + "\n")

    # ---- plot -----------------------------------------------------------
    if not args.no_plot:
        try:
            dotplot(out("dotplot.png"), A, B, a_names, b_names, blocks, best_a,
                    flip, la, lb)
            say(f"Dot plot: {out('dotplot.png')}")
        except ImportError:
            say("Dot plot skipped (matplotlib not available).")
    say(f"Tables written to {args.outdir}/ with prefix '{args.prefix}.'")

    report = "\n".join(lines) + "\n"
    with open(out("summary.txt"), "w") as fh:
        fh.write(report)
    sys.stdout.write(report)


# --------------------------------------------------------------------------
# plotting
# --------------------------------------------------------------------------
CLASS_STYLE = [  # class, colour, legend label
    ("collinear", "#2a78d6", "collinear"),
    ("inverted", "#eb6834", "inverted in place"),
    ("relocated", "#1baf7a", "relocated within scaffold"),
    ("other_scaffold", "#eda100", "in another scaffold's group"),
]


def dotplot(path, A, B, a_names, b_names, blocks, best_a, flip, la, lb):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D

    ink, muted, grid, surface = "#0b0b0b", "#898781", "#e1e0d9", "#fcfcfb"
    a_index = {a: i for i, a in enumerate(a_names)}
    a_off, pos = {}, 0
    for a in a_names:
        a_off[a] = pos
        pos += A.obj_len[a]
    a_total = pos

    # B scaffolds follow their best A scaffold, then position along it
    mid = defaultdict(list)
    for blk in blocks:
        if best_a.get(blk["b_obj"]) == blk["a_obj"]:
            mid[blk["b_obj"]].append((blk["bp"], blk["a_beg"]))

    def b_key(b):
        if b not in best_a:
            return (len(a_names), 0, natural_key(b))
        heavy = max(mid[b])[1] if mid[b] else 0
        return (a_index[best_a[b]], heavy, natural_key(b))

    b_order = sorted(b_names, key=b_key)
    b_off, pos = {}, 0
    for b in b_order:
        b_off[b] = pos
        pos += B.obj_len[b]
    b_total = pos

    segs = defaultdict(list)
    for blk in blocks:
        b, f = blk["b_obj"], flip.get(blk["b_obj"], 1)
        y0, y1 = blk["b_beg"], blk["b_end"]
        if f < 0:
            y0, y1 = B.obj_len[b] - y1 + 1, B.obj_len[b] - y0 + 1
        y0, y1 = y0 + b_off[b], y1 + b_off[b]
        x0, x1 = blk["a_beg"] + a_off[blk["a_obj"]], blk["a_end"] + a_off[blk["a_obj"]]
        if blk["rel"] * f < 0:
            y0, y1 = y1, y0
        segs[blk["cls"]].append([(x0 / 1e6, y0 / 1e6), (x1 / 1e6, y1 / 1e6)])

    fig, ax = plt.subplots(figsize=(12, 12), dpi=150)
    fig.patch.set_facecolor(surface)
    ax.set_facecolor(surface)
    for a in a_names[1:]:
        ax.axvline(a_off[a] / 1e6, color=grid, lw=0.6, zorder=1)
    for b in b_order[1:]:
        ax.axhline(b_off[b] / 1e6, color=grid, lw=0.6, zorder=1)
    handles = []
    for z, (cls, colour, label) in enumerate(CLASS_STYLE):
        if not segs[cls]:
            continue
        ax.add_collection(LineCollection(segs[cls], colors=colour, linewidths=2,
                                         capstyle="round", zorder=5 - z))
        handles.append(Line2D([0], [0], color=colour, lw=2,
                              label=f"{label} ({len(segs[cls])} blocks)"))
    ax.set_xlim(0, a_total / 1e6)
    ax.set_ylim(0, b_total / 1e6)
    ax.set_xticks([(a_off[a] + A.obj_len[a] / 2) / 1e6 for a in a_names])
    ax.set_xticklabels(a_names, rotation=90, fontsize=7, color=ink)
    ax.set_yticks([(b_off[b] + B.obj_len[b] / 2) / 1e6 for b in b_order])
    ax.set_yticklabels([f"{b} ({'+' if flip.get(b, 1) > 0 else '-'})"
                        for b in b_order], fontsize=7, color=ink)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_color("#c3c2b7")
    ax.set_xlabel(f"{la} scaffolds", color=ink, fontsize=10)
    ax.set_ylabel(f"{lb} scaffolds, ordered and oriented to match "
                  f"(- = reversed)", color=ink, fontsize=10)
    ax.set_title(f"Shared-marker blocks: {la} vs {lb}", color=ink, fontsize=12,
                 loc="left")
    leg = ax.legend(handles=handles, loc="upper left", fontsize=8, frameon=True,
                    facecolor=surface, edgecolor=grid)
    for t in leg.get_texts():
        t.set_color(ink)
    ax.text(1.0, 1.005, "each line = one block of markers with identical "
            "order and orientation in both AGPs", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=7, color=muted)
    fig.tight_layout()
    fig.savefig(path, facecolor=surface)
    plt.close(fig)


if __name__ == "__main__":
    main()
