#!/usr/bin/env python3
'''
Make a null Hi-C link set for calibrating hic_phase.py's "mean support".

Within each reference chromosome, the identities of the assigned unitigs are randomly
permuted in the link file, so each chromosome keeps its link count distribution but the
links no longer carry haplotype information. Running hic_phase.py on the result shows how
much support the search achieves on noise; real support should be clearly higher.
Unitig read depths (lines with '.' as the partner) are kept unchanged.

Usage: null_links.py <phasing.tsv> <links.tsv> <out.tsv> [--seed N]
'''
import random, argparse
from collections import defaultdict

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("phasing"); p.add_argument("links"); p.add_argument("out")
    p.add_argument("--seed", type=int, default=11)
    args = p.parse_args()
    random.seed(args.seed)

    chrom = {}
    for l in open(args.phasing):
        s = l.rstrip("\n").split("\t")
        if l.startswith("#") or s[3] != "assigned": continue
        regs, bps = s[6].split(","), [ int(x) for x in s[7].split(",") ]
        chrom[s[0]] = max(zip(bps, regs))[1].split(":")[0] if len(regs) > 1 else regs[0].split(":")[0]
    by = defaultdict(list)
    for q, c in chrom.items(): by[c].append(q)
    perm = {}
    for c, qs in by.items():
        sh = qs[:]; random.shuffle(sh); perm.update(zip(qs, sh))
    with open(args.out, "w") as out:
        for l in open(args.links):
            a, b, n = l.rstrip("\n").split("\t")
            if b == ".": out.write(l); continue
            out.write(f"{perm.get(a, a)}\t{perm.get(b, b)}\t{n}\n")

if __name__ == "__main__":
    main()
