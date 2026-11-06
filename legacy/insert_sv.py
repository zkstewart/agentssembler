#!/usr/bin/env python3
'''
Move individual unitigs to another pseudochromosome, placing each at the gap (and in
the orientation) where its ends have the most Hi-C links to the flanking --window bp.

Moves are given as a TSV of unitig<TAB>target_scaffold. Writes <out>.agp and reports the
chosen position and its score, alongside the second-best position for comparison.
'''

import gzip, argparse
from collections import defaultdict
from rearrange import read_agp, ends, flip

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--agp", required=True)
    p.add_argument("--pairs", required=True)
    p.add_argument("--moves", required=True)
    p.add_argument("--window", type=int, default=1000000)
    p.add_argument("--gap", type=int, default=100)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    comps = read_agp(args.agp)
    moves = [ line.rstrip("\n").split("\t") for line in open(args.moves) if line.strip() and not line.startswith("#") ]
    unknown = [ t for _, t in moves if not t in comps ]
    if unknown:
        raise ValueError(f"target scaffold(s) not in the AGP: {sorted(set(unknown))}")
    moving = { u for u, _ in moves }
    targets = { t for _, t in moves }
    wanted = moving | { u for t in targets for u, _, _ in comps[t] }

    pairs = defaultdict(list)
    with gzip.open(args.pairs, "rt") as fileIn:
        for line in fileIn:
            a, pa, b, pb = line.split("\t")
            if a != b and (a in moving or b in moving) and a in wanted and b in wanted:
                pairs[(a, b)].append((int(pa), int(pb)))

    def links(endX, endY):
        n = 0
        for a, ia in endX.items():
            for b, ib in endY.items():
                for (x0, y0), swap in (((a, b), False), ((b, a), True)):
                    for pa, pb in pairs.get((x0, y0), ()):
                        x, y = (pb, pa) if swap else (pa, pb)
                        if any(lo <= x <= hi for lo, hi in ia) and any(lo <= y <= hi for lo, hi in ib):
                            n += 1
        return n

    # Remove moving unitigs from wherever they are now
    source = {}
    for scaf in comps:
        for c in comps[scaf]:
            if c[0] in moving:
                source[c[0]] = c
        comps[scaf] = [ c for c in comps[scaf] if not c[0] in moving ]

    for u, target in moves:
        _, L, o = source[u]
        half = max(1, L // 2)
        unit = [(u, L, "+")]
        scored = []
        items = comps[target]
        for k in range(len(items) + 1):
            left, right = items[:k], items[k:]
            for orient in ("+", "-"):
                seg = unit if orient == "+" else flip(unit)
                uLeft, uRight = ends(seg, half)
                s = 0
                if left: s += links(ends(left, args.window)[1], uLeft)
                if right: s += links(uRight, ends(right, args.window)[0])
                scored.append((s, k, orient))
        scored.sort(reverse=True)
        s, k, orient = scored[0]
        second = next((x for x in scored[1:] if x[1] != k), (0, None, None))
        comps[target] = items[:k] + [(u, L, orient)] + items[k:]
        where = "start" if k == 0 else ("end" if k == len(items) else f"after {items[k-1][0]}")
        print(f"{u}\t-> {target}\t{where} ({orient})\tlinks={s}\tnext_best={second[0]}")

    with open(args.out + ".agp", "w") as out:
        out.write("##agp-version 2.1\n")
        for scaf, items in comps.items():
            pos, part = 1, 1
            for k, (u, L, o) in enumerate(items):
                if k > 0:
                    ev = "proximity_ligation" if u in moving or items[k-1][0] in moving else "align_genus"
                    out.write(f"{scaf}\t{pos}\t{pos+args.gap-1}\t{part}\tU\t{args.gap}\tscaffold\tyes\t{ev}\n")
                    pos += args.gap; part += 1
                out.write(f"{scaf}\t{pos}\t{pos+L-1}\t{part}\tW\t{u}\t1\t{L}\t{o}\n")
                pos += L; part += 1

if __name__ == "__main__":
    main()
