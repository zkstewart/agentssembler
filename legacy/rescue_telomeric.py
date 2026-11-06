#!/usr/bin/env python3
'''
Rescue "redundant" unitigs that carry a terminal telomere and have a clear Hi-C home.

annotarium marks unitigs whose reference footprint is already covered as "redundant".
Some of these hold real chromosome ends (repeat-rich subtelomeric regions align ambiguously
to several reference chromosome ends). Each redundant telomeric unitig not yet in a group is
added to the (chromosome, haplotype) group that receives at least --minShare of its Hi-C
links, if it has at least --minLinks links to grouped unitigs. build_final.py then places
it at the free telomeric end of that pseudochromosome.

Writes the input membership plus the rescued rows (how = hic_rescued_telomeric) to --out,
and reports every candidate to stderr.
'''
import sys, argparse
from collections import defaultdict

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--membership", required=True, help="make_groups.py membership.tsv")
    p.add_argument("--phasing", required=True, help="hic_phase.py phasing TSV (for 'redundant' status)")
    p.add_argument("--telomeric", required=True, help="telomeres.py calls for telomeric unitigs")
    p.add_argument("--links", required=True)
    p.add_argument("--minLinks", type=int, default=500)
    p.add_argument("--minShare", type=float, default=0.2)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    lines = open(args.membership).readlines()
    lab = {}
    for l in lines:
        if l.startswith("#"): continue
        q, g, _ = l.rstrip("\n").split("\t"); lab[q] = g
    status = { l.split("\t")[0]: l.split("\t")[3] for l in open(args.phasing) if not l.startswith("#") }
    tel = {}
    for l in open(args.telomeric):
        if l.startswith("#"): continue
        sl = l.rstrip("\n").split("\t")
        if sl[4] != "--": tel[sl[0]] = sl[4]
    cand = { q for q in tel if status.get(q) == "redundant" and q not in lab }

    cnt = defaultdict(lambda: defaultdict(int))
    for l in open(args.links):
        a, b, n = l.rstrip("\n").split("\t")
        if b == ".": continue
        for x, y in ((a, b), (b, a)):
            if x in cand and y in lab: cnt[x][lab[y]] += int(n)

    with open(args.out, "w") as out:
        out.writelines(lines)
        for q in sorted(cand):
            tot = sum(cnt[q].values())
            if tot < args.minLinks:
                print(f"{q}\ttoo few links ({tot})", file=sys.stderr); continue
            g, n = max(cnt[q].items(), key=lambda x: x[1])
            if n / tot >= args.minShare:
                out.write(f"{q}\t{g}\thic_rescued_telomeric\n")
                print(f"{q}\t-> {g} ({n/tot:.0%} of {tot} links)", file=sys.stderr)
            else:
                print(f"{q}\tno clear group (best {g} {n/tot:.0%})", file=sys.stderr)

if __name__ == "__main__":
    main()
