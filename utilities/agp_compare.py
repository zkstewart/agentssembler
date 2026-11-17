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

Reference mode (--ref_paf).  One AGP is compared with a reference genome
instead of a second AGP.  The AGP's components are aligned to the reference
beforehand, and the alignments (PAF) say where the reference has each piece
of each component.  That layout takes the role of A, and the AGP is B:
    minimap2 -x asm5 -t 16 reference.fa components.fa > components_vs_ref.paf
    agp_compare.py my.agp --ref_paf components_vs_ref.paf -o out -p name
The PAF queries must be the components named in the AGP (e.g. the unitigs),
not the scaffolds.  Secondary alignments are dropped, alignments shorter than
--min_aln are ignored, where two alignments of a component overlap the better
one keeps the shared part, and alignments with MAPQ below --min_mapq are left
out of the comparison as ambiguous (the component fits two places almost
equally well, as in a region where two haplotypes are identical).

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
  agp_compare.py agentssembler.agp --ref_paf unitigs_vs_ref.paf \\
      --labels reference agentssembler -o out --prefix 13Q010_vs_ref
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
    def __init__(self, path, label, min_components=2, scaffold_regex=None,
                 layout=None, lengths=None):
        self.path, self.label = path, label
        self.entries = defaultdict(list)     # object -> [Comp | Gap] in order
        self.warnings = []
        self.unknown_orient = 0
        if layout is None:
            self._parse()
        else:                                # components placed by alignment
            for obj, es in layout.items():
                self.entries[obj] = sorted(es, key=lambda e: e.obj_beg)
        self.obj_len = {o: max(e.obj_end for e in es)
                        for o, es in self.entries.items()}
        for o in self.obj_len:
            if lengths and o in lengths:
                self.obj_len[o] = max(self.obj_len[o], lengths[o])
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
# reference mode: a layout of the components from their alignments
# --------------------------------------------------------------------------
def open_text(path):
    if path.endswith(".gz"):
        import gzip
        return gzip.open(path, "rt")
    return open(path)


def reference_layout(paf, min_aln, min_mapq):
    """Place component pieces on the reference from a PAF (query = component).

    Returns (layout, reference lengths, counts, ambiguous intervals).  Each
    kept alignment becomes one Comp whose object is the reference sequence.
    Per component, alignments are taken best first (most matching bases) and
    one that overlaps an already taken alignment keeps only its largest free
    part, provided that is still >= min_aln and at least half of it.
    """
    by_query, lengths, stats = defaultdict(list), {}, Counter()
    with open_text(paf) as fh:
        for n, line in enumerate(fh, 1):
            f = line.rstrip("\n").split("\t")
            if len(f) < 12:
                sys.exit(f"{paf}:{n}: expected >=12 tab-separated PAF columns")
            stats["read"] += 1
            lengths[f[5]] = int(f[6])
            if "tp:A:S" in f[12:]:
                stats["secondary"] += 1
                continue
            qs, qe = int(f[2]), int(f[3])
            if qe - qs < min_aln:
                stats["short"] += 1
                continue
            by_query[f[0]].append((int(f[9]), qs, qe, f[4], f[5], int(f[7]),
                                   int(f[8]), int(f[11])))
    layout, ambiguous = defaultdict(list), defaultdict(list)
    for q, alns in by_query.items():
        taken = []
        for _, qs, qe, strand, t, ts, te, mapq in sorted(alns, reverse=True):
            free, pos = [], qs
            for a, b in sorted(taken):
                if b <= pos:
                    continue
                if a >= qe:
                    break
                if a > pos:
                    free.append((pos, a))
                pos = max(pos, b)
            if pos < qe:
                free.append((pos, qe))
            cs, ce = max(free, key=lambda x: x[1] - x[0]) if free else (0, 0)
            if ce - cs < min_aln or (ce - cs) * 2 < qe - qs:
                stats["overlapped"] += 1
                continue
            taken.append((cs, ce))
            if mapq < min_mapq and mapq != 255:
                stats["ambiguous"] += 1
                ambiguous[q].append((cs + 1, ce))
                continue
            scale = (te - ts) / (qe - qs)
            dl, dr = (cs - qs) * scale, (qe - ce) * scale
            if strand == "-":
                dl, dr = dr, dl
            nts, nte = int(round(ts + dl)), int(round(te - dr))
            if nte <= nts:
                continue
            stats["kept"] += 1
            layout[t].append(Comp(t, nts + 1, nte, q, cs + 1, ce,
                                  -1 if strand == "-" else 1))
    return layout, lengths, stats, ambiguous


def covered(intervals):
    """(bp covered at least once, bp of overlap) for 1-based closed intervals."""
    union = total = 0
    end = 0
    for b, e in sorted(intervals):
        total += e - b + 1
        if e > end:
            union += e - max(b, end + 1) + 1
            end = e
    return union, total - union


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
                        # scale is 1 for an AGP; an alignment's span on the
                        # reference can differ from its span on the component
                        scale = (c.obj_end - c.obj_beg + 1) / c.length
                        off = s - c.beg if c.strand > 0 else c.end - e
                        ob = c.obj_beg + int(round(off * scale))
                        oe = ob + max(1, int(round((e - s + 1) * scale))) - 1
                        pl[key] = Place(c.obj, ob, oe, c.strand, c)
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
                "within": pl.comp is pr.comp,
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
    ap.add_argument("agp", nargs="+", metavar="AGP",
                    help="two AGPs: A (supplies scaffold names, e.g. chr1_h1) "
                         "then B (its scaffolds are matched onto A); or one "
                         "AGP with --ref_paf")
    ap.add_argument("--ref_paf", help="reference mode: PAF of the AGP's "
                    "components (query) aligned to a reference genome "
                    "(target). The reference takes the role of A")
    ap.add_argument("--min_aln", type=int, default=10000,
                    help="reference mode: ignore alignments covering fewer "
                         "than this many bp of the component [10000]")
    ap.add_argument("--min_mapq", type=int, default=10,
                    help="reference mode: alignments with a lower MAPQ are "
                         "ambiguous and left out of the comparison [10]")
    ap.add_argument("-o", "--outdir", default=".", help="output directory [.]")
    ap.add_argument("-p", "--prefix", default="agp_compare",
                    help="output file prefix [agp_compare]")
    ap.add_argument("--labels", nargs=2, metavar=("A", "B"),
                    help="short names for A and B [A B; reference assembly "
                         "in reference mode]")
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

    ref_mode = bool(args.ref_paf)
    if len(args.agp) != (1 if ref_mode else 2):
        ap.error("give one AGP with --ref_paf, or two AGPs without it")
    la, lb = args.labels or (["reference", "assembly"] if ref_mode else ["A", "B"])
    if ref_mode:
        layout, ref_len, aln, ambiguous = reference_layout(
            args.ref_paf, args.min_aln, args.min_mapq)
        if not layout:
            sys.exit(f"{args.ref_paf}: no usable alignments")
        # every reference sequence counts as a scaffold unless a regex selects
        A = Agp(args.ref_paf, la, 1, args.scaffold_regex_a, layout, ref_len)
        B = Agp(args.agp[0], lb, args.min_components, args.scaffold_regex_b)
    else:
        A = Agp(args.agp[0], la, args.min_components, args.scaffold_regex_a)
        B = Agp(args.agp[1], lb, args.min_components, args.scaffold_regex_b)
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
        if ref_mode and tag == "A":
            say(f"A = {X.label}: components placed on the reference by "
                f"alignment, {X.path}")
            say(f"    alignments read {aln['read']}: kept {aln['kept']}, "
                f"secondary {aln['secondary']}, shorter than "
                f"{fmt_bp(args.min_aln)} bp {aln['short']},")
            say(f"    overlapped by a better alignment {aln['overlapped']}, "
                f"ambiguous (MAPQ < {args.min_mapq}) {aln['ambiguous']}")
            say(f"    reference sequences with placed components "
                f"{len(X.scaffolds)} ({fmt_bp(sum(X.obj_len[o] for o in X.scaffolds))} bp)")
            continue
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
    say(f"Component IDs: {len(idsA & idsB)} in both A and B, {len(idsA - idsB)} "
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
        what = "placed on the reference" if ref_mode and tag == "A" else "scaffolded"
        say(f"Sequence {what} in {tag} ({fmt_bp(tot)} bp) is, in {other}:")
        for k in ("scaffolded", "unscaffolded", "part not used", "absent"):
            name = "placed" if ref_mode and tag == "B" and k == "scaffolded" else k
            say(f"    {name:<14}{fmt_bp(tally[k]):>16} bp  {pct(tally[k], tot):>8}")
        if ref_mode and tag == "B":
            amb = 0
            for o in B.scaffolds:
                for c in B.comps[o]:
                    amb += sum(max(0, min(c.end, e) - max(c.beg, b) + 1)
                               for b, e in ambiguous.get(c.cid, ()))
            unplaced = tally["part not used"] + tally["absent"]
            say(f"    of the {fmt_bp(unplaced)} bp not placed: {fmt_bp(amb)} bp is "
                f"ambiguous (MAPQ < {args.min_mapq}),")
            say(f"    {fmt_bp(unplaced - amb)} bp has no alignment of "
                f"{fmt_bp(args.min_aln)} bp or more")
    if ref_mode:
        say("  In A: part not used = the component is in B, but not this part of it (e.g.")
        say("  a trimmed overlap); absent = the component is not in B at all.")
        say("  In B: part not used / absent = this part of the component, or all of it,")
        say("  has no unique placement on the reference. Such sequence is not compared.")
    else:
        say("  part not used = the component is in the other AGP, but not this part of it")
        say("  (e.g. a trimmed overlap); absent = the whole component is missing there.")
    markers = [k for k in plA if k in plB
               and A.is_scaffold(plA[k].obj) and B.is_scaffold(plB[k].obj)]
    total_bp = sum(key_len(k) for k in markers)
    say(f"Shared markers scaffolded in both: {len(markers)} "
        f"({fmt_bp(total_bp)} bp). All metrics below use these.")
    if ref_mode:
        spans = defaultdict(list)
        for k in markers:
            spans[plA[k].obj].append((plA[k].beg, plA[k].end))
        cov = [covered(v) for v in spans.values()]
        ref_total = sum(A.obj_len[o] for o in A.scaffolds)
        once, twice = sum(c[0] for c in cov), sum(c[1] for c in cov)
        say(f"Reference covered by B's scaffolds: {fmt_bp(once)} of "
            f"{fmt_bp(ref_total)} bp ({pct(once, ref_total)}); covered more than "
            f"once: {fmt_bp(twice)} bp")
        say("  (reference sequences with no alignment at all are not in the PAF and so")
        say("  are not counted)")
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
        inside = sum(1 for j in js if j["within"] and j["status"] != "adjacent")
        if inside or ref_mode:
            say(f"    of the {len(js) - c['adjacent']} not adjacent, {inside} lie "
                f"inside one {tag} component, not at a join")
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
        if ref_mode and tag == "A":
            say("A: the reference layout has no gaps (runs of N in the reference "
                "are not read)")
            continue
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
                               size="NA", verdict="not_assessable")
                else:
                    o = j["other"]
                    oobjs = sorted({plY[j["left"]].obj, plY[j["right"]].obj},
                                   key=natural_key)
                    row.update(status=j["status"], lm=key_str(j["left"]),
                               rm=key_str(j["right"]), oobj=",".join(oobjs),
                               ogaps=len(o["gaps"]) if o else "NA",
                               ocomps=o["n_comps"] if o else "NA",
                               obp=o["seq_bp"] if o else "NA", size="NA")
                    if j["status"] == "adjacent":
                        row["verdict"] = ("gap_in_both" if o["gaps"]
                                          else "join_gapless_in_other")
                        if not o["gaps"] and len(j["gaps"]) == 1:
                            # sequence the other side has between the markers,
                            # less what this side already has there
                            row["size"] = o["seq_bp"] - j["seq_bp"]
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
            d = sorted(r["size"] for r in filled if r["size"] != "NA")
            if d:
                say(f"    size of the {len(d)} gaps that {other} joins without a gap, "
                    f"implied by {other}: median {fmt_bp(d[len(d) // 2])} bp,")
                say(f"    range {fmt_bp(d[0])} to {fmt_bp(d[-1])} bp (negative = the "
                    f"two sides overlap in {other})")
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
    if ref_mode:      # a reference sequence is measured by its own length
        a_cbp = {a: A.obj_len[a] for a in a_names}
    else:
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
            "source", "scaffold", "left_marker", "right_marker",
            "within_component", "gaps", "components_between", "sequence_bp_between", "status_in_other",
            "other_scaffold_left", "other_scaffold_right", "other_gaps",
            "other_components_between", "other_sequence_bp_between"]) + "\n")
        for js, X, plY in ((jA, A, plB), (jB, B, plA)):
            for j in js:
                o = j["other"]
                fh.write("\t".join(map(str, [
                    X.label, j["obj"], key_str(j["left"]), key_str(j["right"]),
                    "yes" if j["within"] else "no", len(j["gaps"]), j["n_comps"], j["seq_bp"], j["status"],
                    plY[j["left"]].obj, plY[j["right"]].obj,
                    len(o["gaps"]) if o else "NA", o["n_comps"] if o else "NA",
                    o["seq_bp"] if o else "NA"])) + "\n")
    with open(out("gaps.tsv"), "w") as fh:
        fh.write("\t".join([
            "source", "scaffold", "gap_beg", "gap_end", "gap_length", "gap_type",
            "left_component", "right_component", "left_shared_marker",
            "right_shared_marker", "junction_status_in_other", "verdict",
            "other_scaffold", "other_gaps", "other_components_between",
            "other_sequence_bp_between", "implied_gap_size"]) + "\n")
        for r in gap_rows:
            fh.write("\t".join(map(str, [
                r["src"], r["obj"], r["beg"], r["end"], r["len"], r["type"],
                r["left"], r["right"], r["lm"], r["rm"], r["status"],
                r["verdict"], r["oobj"], r["ogaps"], r["ocomps"], r["obp"],
                r["size"]])) + "\n")

    # ---- plot -----------------------------------------------------------
    if not args.no_plot:
        try:
            dotplot(out("dotplot.png"), A, B, a_names, b_names, blocks, best_a,
                    flip, la, lb, ref_mode)
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


def dotplot(path, A, B, a_names, b_names, blocks, best_a, flip, la, lb,
            ref_mode=False):
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
    # name only the sequences long enough for a label to fit
    a_lab = [a for a in a_names if A.obj_len[a] >= 0.004 * a_total]
    b_lab = [b for b in b_order if B.obj_len[b] >= 0.004 * b_total]
    ax.set_xticks([(a_off[a] + A.obj_len[a] / 2) / 1e6 for a in a_lab])
    ax.set_xticklabels(a_lab, rotation=90, fontsize=7, color=ink)
    ax.set_yticks([(b_off[b] + B.obj_len[b] / 2) / 1e6 for b in b_lab])
    ax.set_yticklabels([f"{b} ({'+' if flip.get(b, 1) > 0 else '-'})"
                        for b in b_lab], fontsize=7, color=ink)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_color("#c3c2b7")
    ax.set_xlabel(f"{la} {'sequences' if ref_mode else 'scaffolds'}",
                  color=ink, fontsize=10)
    ax.set_ylabel(f"{lb} scaffolds, ordered and oriented to match "
                  f"(- = reversed)", color=ink, fontsize=10)
    ax.set_title(f"Shared-marker blocks: {la} vs {lb}", color=ink, fontsize=12,
                 loc="left")
    leg = ax.legend(handles=handles, loc="upper left", fontsize=8, frameon=True,
                    facecolor=surface, edgecolor=grid)
    for t in leg.get_texts():
        t.set_color(ink)
    ax.text(1.0, 1.005, "each line = one block of markers with identical "
            "order and orientation in both", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=7, color=muted)
    fig.tight_layout()
    fig.savefig(path, facecolor=surface)
    plt.close(fig)


if __name__ == "__main__":
    main()
