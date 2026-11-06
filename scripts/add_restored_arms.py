#!/usr/bin/env python3
'''Restored arms become chains of their own in the broken-chain AGP, named chain9001 onward. They have no Hi-C pairs, being near-identical to a sibling.
Usage: add_restored_arms.py restored_arms.tsv unitigs.fasta.fai chains.agp restored_chains.tsv'''
import sys
arms, fai, agp, out = sys.argv[1:5]
L = {l.split()[0]: int(l.split()[1]) for l in open(fai)}
U = [l.split("\t")[0] for l in open(arms) if l[0] != "#"]
with open(agp, "a") as o, open(out, "w") as t:
    for i, u in enumerate(U, 1):
        c = f"chain9{i:03d}"; o.write(f"{c}\t1\t{L[u]}\t1\tW\t{u}\t1\t{L[u]}\t+\n"); t.write(f"{c}\t{u}\n")
print(len(U), "restored arms")
