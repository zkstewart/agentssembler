#!/usr/bin/env python3
'''
Finish an AGP produced by rearrange.py / insert_sv.py into a deliverable:

1. Telomere capping: a moved unitig (from the moves file) that carries a terminal telomere
   is moved to a pseudochromosome end lacking a telomere, telomere outward (the same rule
   build_final.py applies to Hi-C-placed unitigs).
2. Pseudochromosomes are written in natural order (chr1_h1, chr1_h2, ..., chr10_h1).
3. Gap evidence: a join that already existed in the main assembly keeps that AGP's evidence;
   every new join (translocation junctions, moved unitigs) is 'proximity_ligation'.
'''
import argparse
from rearrange import read_agp
from build_final import cap_telomeres
from common import natural_key

def joins(agp):
    out, prev, gap = {}, None, None
    for l in open(agp):
        if l.startswith("#"): continue
        x = l.rstrip("\n").split("\t")
        if x[4] == "W":
            if prev and prev[0] == x[0]: out[frozenset((prev[1], x[5]))] = gap
            prev = (x[0], x[5])
        else:
            gap = x[8]
    return out

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--agp", required=True, help="AGP after rearrange.py/insert_sv.py")
    p.add_argument("--main", required=True, help="main assembly AGP (source of join evidence)")
    p.add_argument("--moves", default=None, help="insert_sv.py moves TSV")
    p.add_argument("--telomeric", required=True, help="telomeres.py calls for the unitigs")
    p.add_argument("--gap", type=int, default=100)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    comps = read_agp(args.agp)
    moved = set()
    if args.moves:
        moved = { l.split("\t")[0] for l in open(args.moves) if l.strip() and not l.startswith("#") }
    telomeric = {}
    for l in open(args.telomeric):
        if l.startswith("#"): continue
        sl = l.rstrip("\n").split("\t")
        telomeric[sl[0]] = sl[4]

    for scaf, items in comps.items():
        merged = [ (u, o, "hic" if u in moved else "ref") for u, L, o in items ]
        merged, moves = cap_telomeres(merged, telomeric)
        for u, end in moves:
            print(f"{scaf}\t{u}\tcapped {end}")
        lengths = { u: L for u, L, _ in items }
        comps[scaf] = [ (u, lengths[u], o) for u, o, _ in merged ]

    evidence = joins(args.main)
    new = 0
    with open(args.out, "w") as out:
        out.write("##agp-version 2.1\n")
        for scaf in sorted(comps, key=natural_key):
            pos, part = 1, 1
            for k, (u, L, o) in enumerate(comps[scaf]):
                if k > 0:
                    ev = evidence.get(frozenset((comps[scaf][k-1][0], u)))
                    if ev is None:
                        ev, new = "proximity_ligation", new + 1
                    out.write(f"{scaf}\t{pos}\t{pos+args.gap-1}\t{part}\tU\t{args.gap}\tscaffold\tyes\t{ev}\n")
                    pos += args.gap; part += 1
                out.write(f"{scaf}\t{pos}\t{pos+L-1}\t{part}\tW\t{u}\t1\t{L}\t{o}\n")
                pos += L; part += 1
    print(f"new joins (not in the main assembly): {new}")

if __name__ == "__main__":
    main()
