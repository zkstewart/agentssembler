#!/usr/bin/env python3
'''
Re-phase reference-placed unitigs into haplotypes using Hi-C links.

Each unitig is placed on its primary reference chromosome using annotarium's footprint
logic. Within a chromosome, two unitigs "conflict" (are allelic) when their reference
footprints overlap substantially, and conflicting unitigs may not share a haplotype.
Subject to that constraint, unitigs are coloured into `ploidy` haplotypes to maximise
the summed normalised Hi-C links between unitigs of the same haplotype.

The search is a beam search sweeping along the chromosome in footprint order, followed
by single-unitig recolouring until no move improves the score.

Outputs a phasing TSV in annotarium's format (so `annotarium assembly scaffold` can use it)
and a per-chromosome report comparing the Hi-C support of the input and new phasings.
'''

import os, sys, argparse, math
import numpy as np
from collections import defaultdict

# annotarium (github.com/zkstewart/annotarium) supplies the PAF parsing and footprint logic
sys.path.insert(0, os.path.join(os.environ["ANNOTARIUM_DIR"], "modules"))
from minimap2 import parse_paf
from phasing import deduplicate_query, build_footprints, infer_orientations, merge_intervals, intervals_length

def parse_phasing_tsv(fileName):
    rows = {}
    with open(fileName) as fileIn:
        header = fileIn.readline()
        for line in fileIn:
            sl = line.rstrip("\r\n").split("\t")
            rows[sl[0]] = sl
    return header, rows

def parse_links(fileName):
    '''
    Expects lines of "contigA<TAB>contigB<TAB>count" for inter-contig pairs and
    "contig<TAB>.<TAB>count" for each contig's total filtered read count.
    '''
    links, depth = {}, {}
    with open(fileName) as fileIn:
        for line in fileIn:
            a, b, n = line.rstrip("\r\n").split("\t")
            if b == ".":
                depth[a] = int(n)
            else:
                key = (a, b) if a < b else (b, a)
                links[key] = links.get(key, 0) + int(n)
    return links, depth

def overlap_bp(intsA, intsB):
    i = j = total = 0
    while i < len(intsA) and j < len(intsB):
        s = max(intsA[i][0], intsB[j][0]); e = min(intsA[i][1], intsB[j][1])
        if e > s:
            total += e - s
        if intsA[i][1] < intsB[j][1]:
            i += 1
        else:
            j += 1
    return total

class ChromosomePhaser:
    EVICT_MAX_BP = 500000
    EVICT_MAX_FRAC = 0.25

    def __init__(self, nodes, intervals, weights, ploidy, minOverlap, minOverlapFrac, beamWidth, lengths):
        self.nodes = nodes # list of unitig IDs, sorted by footprint start
        self.lengths = lengths
        self.evicted = set()
        self.intervals = intervals
        self.w = weights # dict of node -> {node: weight}
        self.ploidy = ploidy
        self.beamWidth = beamWidth

        # Find conflicting (allelic) pairs
        self.conflicts = defaultdict(set)
        for x in range(len(nodes)):
            a = nodes[x]
            aEnd = intervals[a][-1][1]
            for y in range(x+1, len(nodes)):
                b = nodes[y]
                if intervals[b][0][0] >= aEnd:
                    continue # footprints are not sorted by end, so keep scanning
                ov = overlap_bp(intervals[a], intervals[b])
                shorter = min(intervals_length(intervals[a]), intervals_length(intervals[b]))
                if ov >= minOverlap and ov >= minOverlapFrac * shorter:
                    self.conflicts[a].add(b)
                    self.conflicts[b].add(a)

    def score(self, states):
        total = 0.0
        for a in self.nodes:
            for b, wt in self.w.get(a, {}).items():
                if a < b and b in states and states.get(a) is not None and states[a] == states[b]:
                    total += wt
        return total

    def violations(self, states):
        return sum(1 for a in self.nodes for b in self.conflicts[a]
                   if a < b and states.get(a) is not None and states.get(a) == states.get(b))

    def beam(self):
        index = { n: i for i, n in enumerate(self.nodes) }
        beam = [ (0.0, ()) ] # (score, tuple of colours for nodes[0:k], None = unplaced)
        for k, node in enumerate(self.nodes):
            earlier = [ (index[b], wt) for b, wt in self.w.get(node, {}).items() if index.get(b, len(self.nodes)) < k ]
            earlierConflicts = [ index[b] for b in self.conflicts[node] if index[b] < k ]
            candidates = []
            for score, colours in beam:
                used = { colours[i] for i in earlierConflicts }
                # Break label symmetry: only allow a new colour index one above the max used
                maxColour = max([ c for c in colours if c is not None ], default=-1)
                options = [ c for c in range(min(self.ploidy, maxColour + 2)) if not c in used ]
                for c in options:
                    gain = sum(wt for i, wt in earlier if colours[i] == c)
                    candidates.append((score + gain, colours + (c,)))
                if options == []:
                    candidates.append((score - 1e6, colours + (None,)))
            candidates.sort(key=lambda x: -x[0])
            beam = candidates[:self.beamWidth]
        best = beam[0][1]
        return { n: best[i] for i, n in enumerate(self.nodes) }

    def refine(self, states):
        improved = True
        while improved:
            improved = False
            for node in self.nodes:
                cur = states[node]
                used = { states[b] for b in self.conflicts[node] }
                gains = defaultdict(float)
                for b, wt in self.w.get(node, {}).items():
                    if b in states and states[b] is not None:
                        gains[states[b]] += wt
                if cur is None:
                    # Reinsert unplaced unitigs wherever they fit, even without links
                    options = [ c for c in range(self.ploidy) if not c in used ]
                    if options != []:
                        states[node] = max(options, key=lambda c: gains[c])
                        improved = True
                    continue
                for c in range(self.ploidy):
                    if c != cur and not c in used and gains[c] > gains.get(cur, 0.0) + 1e-12:
                        states[node] = c
                        cur = c
                        improved = True
        return states

    def _matrices(self, rowNodes, colNodes):
        W = np.zeros((len(rowNodes), len(colNodes)))
        K = np.zeros((len(rowNodes), len(colNodes)))
        colIndex = { n: j for j, n in enumerate(colNodes) }
        for i, a in enumerate(rowNodes):
            for b, wt in self.w.get(a, {}).items():
                if b in colIndex:
                    W[i, colIndex[b]] = wt
            for b in self.conflicts[a]:
                if b in colIndex:
                    K[i, colIndex[b]] = 1
        return W, K

    @staticmethod
    def _head_tail(M):
        '''
        Returns HT where HT[i, j] = M[:i, j:].sum() for i in 0..n1 and j in 0..n2.
        '''
        n1, n2 = M.shape
        P = np.zeros((n1 + 1, n2 + 1))
        P[1:, 1:] = M.cumsum(0).cumsum(1) # P[i, j] = M[:i, :j].sum()
        return P[:, n2:n2+1] - P # M[:i, :].sum() - M[:i, :j].sum()

    def tail_swap_refine(self, states):
        '''
        Repairs switches between phase blocks. Each haplotype is a thread of unitigs in
        footprint order; for each pair of haplotypes, every pair of junctions (one per thread)
        is considered for exchanging the two threads' tails. The swap that most raises the
        score while keeping the colouring valid is applied, repeating until none helps.
        Because the junctions may differ between threads, this respects the staggered
        overlaps of unitigs, and lets long-range Hi-C links decide the phase across regions
        that local links cannot bridge.
        '''
        improved = True
        while improved:
            improved = False
            for c1 in range(self.ploidy):
                for c2 in range(c1 + 1, self.ploidy):
                    T1 = [ n for n in self.nodes if states[n] == c1 ]
                    T2 = [ n for n in self.nodes if states[n] == c2 ]
                    if T1 == [] or T2 == []:
                        continue
                    W11, _ = self._matrices(T1, T1)
                    W22, _ = self._matrices(T2, T2)
                    W12, K12 = self._matrices(T1, T2)

                    # Heads are T1[:i] and T2[:j]; tails are T1[i:] and T2[j:]
                    keep1 = np.diag(self._head_tail(W11))[:, None] # W(T1[:i], T1[i:]) as a column
                    keep2 = np.diag(self._head_tail(W22))[None, :] # W(T2[:j], T2[j:]) as a row
                    gain = self._head_tail(W12) + self._head_tail(W12.T).T - keep1 - keep2
                    clash = self._head_tail(K12) + self._head_tail(K12.T).T
                    gain[clash > 0] = -np.inf
                    gain[0, 0] = gain[-1, -1] = -np.inf # whole-thread swaps are relabelling only
                    i, j = np.unravel_index(np.argmax(gain), gain.shape)
                    if gain[i, j] > 1e-9:
                        for n in T1[i:]:
                            states[n] = c2
                        for n in T2[j:]:
                            states[n] = c1
                        improved = True
        return states

    def segment_swap_refine(self, states, maxLen=8):
        '''
        Repairs local misassignments that need several unitigs to move at once: for each
        pair of haplotypes, a run of up to maxLen consecutive unitigs in one thread is
        exchanged with an overlapping run in the other, when the result is valid and
        raises the score. Repeats until no exchange helps.
        '''
        improved = True
        while improved:
            improved = False
            for c1 in range(self.ploidy):
                for c2 in range(c1 + 1, self.ploidy):
                    T1 = [ n for n in self.nodes if states[n] == c1 ]
                    T2 = [ n for n in self.nodes if states[n] == c2 ]
                    span = lambda seg: (self.intervals[seg[0]][0][0], max(self.intervals[n][-1][1] for n in seg))
                    best, bestMove = 1e-9, None
                    for i1 in range(len(T1)):
                        for i2 in range(i1 + 1, min(len(T1), i1 + maxLen) + 1):
                            A = T1[i1:i2]; aStart, aEnd = span(A)
                            for j1 in range(len(T2)):
                                if self.intervals[T2[j1]][0][0] >= aEnd:
                                    break
                                for j2 in range(j1 + 1, min(len(T2), j1 + maxLen) + 1):
                                    B = T2[j1:j2]; bStart, bEnd = span(B)
                                    if bEnd <= aStart:
                                        continue
                                    result = self._swap_gain(A, B, states, c1, c2)
                                    if result is not None and result[0] > best:
                                        best, bestMove = result[0], (A, B, result[1])
                    if bestMove is not None:
                        A, B, evicted = bestMove
                        for n in evicted:
                            states[n] = None
                            self.evicted.add(n)
                        for n in A:
                            states[n] = c2
                        for n in B:
                            states[n] = c1
                        improved = True
        return states

    def _swap_gain(self, A, B, states, c1, c2):
        '''
        Score change from moving A (colour c1) to c2 and B (colour c2) to c1, or None if
        the move would put conflicting unitigs in the same haplotype. Small blocking
        unitigs (e.g. a surplus fourth allele) may be evicted to unplaced, if their total
        length is within EVICT_MAX_BP and EVICT_MAX_FRAC of the moved length.
        Returns (gain, evicted) or None.
        '''
        setA, setB = set(A), set(B)
        blockers = { b for a in A for b in self.conflicts[a] if states[b] == c2 and not b in setB }
        blockers |= { a for b in B for a in self.conflicts[b] if states[a] == c1 and not a in setA }
        if blockers:
            blockLen = sum(self.lengths[x] for x in blockers)
            movedLen = sum(self.lengths[x] for x in A) + sum(self.lengths[x] for x in B)
            if blockLen > self.EVICT_MAX_BP or blockLen > self.EVICT_MAX_FRAC * movedLen:
                return None
        gain = 0.0
        for a in A:
            for x, wt in self.w.get(a, {}).items():
                sx = states.get(x)
                if x in blockers:
                    continue
                if sx == c2 and not x in setB:
                    gain += wt
                elif sx == c1 and not x in setA:
                    gain -= wt
        for b in B:
            for x, wt in self.w.get(b, {}).items():
                sx = states.get(x)
                if x in blockers:
                    continue
                if sx == c1 and not x in setA:
                    gain += wt
                elif sx == c2 and not x in setB:
                    gain -= wt
        # Evicted unitigs lose their links to the rest of their haplotype
        for x in blockers:
            for y, wt in self.w.get(x, {}).items():
                if states.get(y) == states[x] and not y in setA and not y in setB and not y in blockers:
                    gain -= wt
        return gain, blockers

    def support(self, states):
        '''
        For each node, the share of its weight to same-chromosome nodes that goes to
        its own haplotype; ~1/ploidy indicates no Hi-C evidence for the assignment.
        '''
        result = {}
        for node in self.nodes:
            tot = own = 0.0
            for b, wt in self.w.get(node, {}).items():
                if b in states and states[b] is not None:
                    tot += wt
                    if states[b] == states[node]:
                        own += wt
            result[node] = own / tot if tot > 0 else float("nan")
        return result

def core_footprint(intervals, trim=0.05):
    '''
    Clips a unitig's footprint on one chromosome to the span holding the central
    (1 - 2*trim) of its aligned bases, so a few stray alignments elsewhere on the
    chromosome do not distort its position or create false allelic conflicts.
    '''
    total = intervals_length(intervals)
    lo, hi = trim * total, (1 - trim) * total
    cum, start, end = 0, None, None
    for s, e in intervals:
        if start is None and cum + (e - s) > lo:
            start = s + max(0, int(lo - cum))
        if end is None and cum + (e - s) >= hi:
            end = s + int(hi - cum)
        cum += e - s
    core = [ (max(s, start), min(e, end)) for s, e in intervals if e > start and s < end ]
    return core if core != [] else intervals

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--paf", required=True, help="unitig-to-reference PAF")
    p.add_argument("--phasing", required=True, help="annotarium phasing.tsv to start from")
    p.add_argument("--links", required=True, help="Hi-C link counts (see count_links.sh)")
    p.add_argument("--out", required=True, help="output prefix")
    p.add_argument("--ploidy", type=int, default=3)
    p.add_argument("--minQ", type=int, default=1)
    p.add_argument("--minAlignLen", type=int, default=10000)
    p.add_argument("--minOverlap", type=int, default=50000)
    p.add_argument("--minOverlapFrac", type=float, default=0.1)
    p.add_argument("--beam", type=int, default=2000)
    args = p.parse_args()

    header, rows = parse_phasing_tsv(args.phasing)
    alignments = deduplicate_query(parse_paf(args.paf, minQ=args.minQ, minAlignLen=args.minAlignLen))
    footprints = build_footprints(alignments)
    links, depth = parse_links(args.links)

    # Place each assigned unitig on its primary chromosome
    byChrom = defaultdict(list)
    primary = {}
    for qid, sl in rows.items():
        if sl[3] != "assigned" or not qid in footprints:
            continue
        tid = max(footprints[qid], key=lambda t: intervals_length(footprints[qid][t]))
        primary[qid] = tid
        byChrom[tid].append(qid)

    # Normalised weights: links / sqrt(depth_a * depth_b), restricted to same chromosome
    weights = defaultdict(dict)
    for (a, b), n in links.items():
        if a in primary and b in primary and primary[a] == primary[b]:
            da, db = depth.get(a, 0), depth.get(b, 0)
            if da > 0 and db > 0:
                wt = n / math.sqrt(da * db)
                weights[a][b] = wt
                weights[b][a] = wt

    newStates, supportAll, coreSpan, evictedAll = {}, {}, {}, set()
    report = open(args.out + ".report.tsv", "w")
    report.write("#chrom\tunitigs\tconflict_pairs\tinput_score\tinput_violations\tinput_mean_support\tnew_score\tnew_violations\tnew_mean_support\tchanged_unitigs\n")
    for tid in sorted(byChrom):
        nodes = byChrom[tid]
        intervals = { q: core_footprint(footprints[q][tid]) for q in nodes }
        coreSpan.update({ q: (iv[0][0], iv[-1][1]) for q, iv in intervals.items() })
        nodes.sort(key=lambda q: (intervals[q][0][0], q))
        lengths = { q: int(rows[q][1]) for q in nodes }
        phaser = ChromosomePhaser(nodes, intervals, weights, args.ploidy, args.minOverlap, args.minOverlapFrac, args.beam, lengths)

        inStates = { q: int(rows[q][2]) - 1 for q in nodes }
        states = phaser.beam()
        prev = None
        for _ in range(20):
            if prev == states:
                break
            prev = dict(states)
            states = phaser.refine(phaser.segment_swap_refine(phaser.tail_swap_refine(states)))

        # Align new labels to the input labels by maximum length agreement, so "changed"
        # counts reflect real haplotype changes rather than label permutations
        from itertools import permutations
        bestPerm = max(permutations(range(args.ploidy)),
                       key=lambda perm: sum(lengths[q] for q in nodes if states[q] is not None and perm[states[q]] == inStates[q]))
        states = { q: (bestPerm[s] if s is not None else None) for q, s in states.items() }

        inSupport = phaser.support(inStates)
        newSupport = phaser.support(states)
        mean = lambda d: sum(v for v in d.values() if v == v) / max(1, sum(1 for v in d.values() if v == v))
        numConflicts = sum(len(v) for v in phaser.conflicts.values()) // 2
        changed = sum(1 for q in nodes if states[q] != inStates[q])
        report.write(f"{tid}\t{len(nodes)}\t{numConflicts}\t{phaser.score(inStates):.4f}\t{phaser.violations(inStates)}\t{mean(inSupport):.3f}\t" +
                     f"{phaser.score(states):.4f}\t{phaser.violations(states)}\t{mean(newSupport):.3f}\t{changed}\n")
        newStates.update(states)
        evictedAll |= { q for q in phaser.evicted if states[q] is None }
        supportAll.update(newSupport)
    report.close()

    with open(args.out + ".phasing.tsv", "w") as fileOut:
        fileOut.write(header)
        for qid, sl in rows.items():
            sl = list(sl)
            if qid in newStates:
                if newStates[qid] is None:
                    sl[2], sl[3] = ".", ("evicted" if qid in evictedAll else "unplaceable")
                else:
                    sl[2] = str(newStates[qid] + 1)
            fileOut.write("\t".join(sl) + "\n")

    with open(args.out + ".support.tsv", "w") as fileOut:
        fileOut.write("#contig\tchrom\tcore_start\tcore_end\tlength\thaplotype\tsupport\tdepth\n")
        for qid in newStates:
            fileOut.write(f"{qid}\t{primary[qid]}\t{coreSpan[qid][0]}\t{coreSpan[qid][1]}\t{rows[qid][1]}\t{'.' if newStates[qid] is None else newStates[qid]+1}\t{supportAll[qid]:.3f}\t{depth.get(qid, 0)}\n")

if __name__ == "__main__":
    main()
