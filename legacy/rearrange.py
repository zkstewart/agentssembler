#!/usr/bin/env python3
'''
Rebuild pseudochromosomes from segments of existing ones, choosing segment order and
orientation by Hi-C.

A spec file lists one derived pseudochromosome per line:
    name<TAB>scaffold:firstUnitig..lastUnitig[,scaffold:u..v ...]
Each segment is the run of AGP components from firstUnitig to lastUnitig (inclusive) of
that scaffold, in AGP order. A segment may be a single unitig (scaffold:u..u). For each
derived pseudochromosome every order (the first segment's orientation is fixed) and
orientation of its segments is tried, scoring the Hi-C pairs that link the facing
--window bp at each junction; the best arrangement is written out.

Scaffolds named in the spec are replaced; all others are copied unchanged. Writes
<out>.agp and a report of junction scores to stdout.
'''

import gzip, argparse
from itertools import permutations, product
from collections import defaultdict

def read_agp(fileName):
    comps = defaultdict(list)
    for line in open(fileName):
        if line.startswith("#"): continue
        sl = line.rstrip("\n").split("\t")
        if sl[4] == "W":
            comps[sl[0]].append((sl[5], int(sl[7]) - int(sl[6]) + 1, sl[8]))
    return comps

def segment(comps, spec):
    scaf, rng = spec.split(":")
    first, last = rng.split("..")
    names = [ c[0] for c in comps[scaf] ]
    i, j = names.index(first), names.index(last)
    return comps[scaf][i:j+1]

def flip(seg):
    return [ (u, L, "-" if o == "+" else "+") for u, L, o in reversed(seg) ]

def ends(seg, window):
    '''Unitig intervals within window bp of each end: {unitig: [(lo, hi), ...]} in unitig
    coordinates (1-based), for the 'left' and 'right' ends of the oriented segment.'''
    def take(items):
        acc, out = 0, defaultdict(list)
        for u, L, o in items:
            need = window - acc
            if need <= 0: break
            span = min(L, need)
            # the part of u nearest this end of the segment
            out[u].append((1, span) if o == "+" else (L - span + 1, L))
            acc += L
        return out
    left = take(seg)
    right = take([ (u, L, "-" if o == "+" else "+") for u, L, o in reversed(seg) ])
    return left, right

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--agp", required=True)
    p.add_argument("--pairs", required=True)
    p.add_argument("--spec", required=True)
    p.add_argument("--window", type=int, default=1000000)
    p.add_argument("--gap", type=int, default=100)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    comps = read_agp(args.agp)
    derived, replaced = [], set()
    for line in open(args.spec):
        if not line.strip() or line.startswith("#"): continue
        name, segs = line.rstrip("\n").split("\t")
        segs = segs.split(",")
        replaced |= { s.split(":")[0] for s in segs }
        derived.append((name, [ segment(comps, s) for s in segs ], segs))

    # Every component of a replaced scaffold must be used exactly once, or sequence is lost
    used = [ u for _, segs, _ in derived for seg in segs for u, _, _ in seg ]
    dups = sorted({ u for u in used if used.count(u) > 1 })
    missing = sorted({ c[0] for s in replaced for c in comps[s] } - set(used))
    if dups or missing:
        raise ValueError(f"spec leaves out {missing} and/or repeats {dups} from the replaced scaffolds")

    # Every unitig that sits in any segment, so the pairs file is read only once
    wanted = { u for _, segs, _ in derived for seg in segs for u, _, _ in seg }
    pairs = defaultdict(list)
    with gzip.open(args.pairs, "rt") as fileIn:
        for line in fileIn:
            a, pa, b, pb = line.split("\t")
            if a in wanted and b in wanted and a != b:
                pairs[(a, b)].append((int(pa), int(pb)))

    def links(endX, endY):
        n = 0
        for a, ia in endX.items():
            for b, ib in endY.items():
                for key, swap in (((a, b), False), ((b, a), True)):
                    for pa, pb in pairs.get(key, ()):
                        x, y = (pb, pa) if swap else (pa, pb)
                        if any(lo <= x <= hi for lo, hi in ia) and any(lo <= y <= hi for lo, hi in ib):
                            n += 1
        return n

    out = open(args.out + ".agp", "w")
    out.write("##agp-version 2.1\n")
    def write(name, items, joins):
        pos, part = 1, 1
        for k, (u, L, o) in enumerate(items):
            if k > 0:
                ev = "proximity_ligation" if k in joins else "align_genus"
                out.write(f"{name}\t{pos}\t{pos+args.gap-1}\t{part}\tU\t{args.gap}\tscaffold\tyes\t{ev}\n")
                pos += args.gap; part += 1
            out.write(f"{name}\t{pos}\t{pos+L-1}\t{part}\tW\t{u}\t1\t{L}\t{o}\n")
            pos += L; part += 1

    for scaf in comps:
        if not scaf in replaced:
            write(scaf, comps[scaf], set())
    for name, segs, labels in derived:
        best = None
        n = len(segs)
        for order in permutations(range(n)):
            for orient in product((False, True), repeat=n):
                if orient[0]:
                    continue # fix segment 0's orientation; the mirror image scores the same
                oriented = [ flip(segs[i]) if orient[i] else segs[i] for i in order ]
                score, junction = 0, []
                for x, y in zip(oriented, oriented[1:]):
                    s = links(ends(x, args.window)[1], ends(y, args.window)[0])
                    score += s; junction.append(s)
                if best is None or score > best[0]:
                    best = (score, order, orient, oriented, junction)
        score, order, orient, oriented, junction = best
        items, joins, k = [], set(), 0
        for seg in oriented:
            if items:
                joins.add(len(items))
            items.extend(seg)
        write(name, items, joins)
        desc = " + ".join(f"{labels[i]}{'(rev)' if orient[i] else ''}" for i in order)
        print(f"{name}\t{sum(L for _, L, _ in items)/1e6:.1f} Mb\t{desc}\tjunction_links={junction}")
    out.close()

if __name__ == "__main__":
    main()
