#!/usr/bin/env python3
'''
Names for groups: <chrom>_h<k>. Names are conventions and change no layout. With a reference, <chrom> is the reference sequence that most of the chromosome's
aligned sequence hits, and a scaffold is flagged to be written in the reference's direction. Without one, chromosomes are named chr1.. by size. Haplotype numbers
are arbitrary and mean nothing across chromosomes. Check for duplicates in the output: a reference name claimed by two chromosomes is an alarm, not a choice.
Usage: auto_names.py --groups groups.tsv --chroms groups.chroms.tsv --out names.tsv [--order order.tsv --agp chains.agp --paf unitigs_vs_ref.paf.gz --prefix chr]
'''
import argparse, collections, gzip, re
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in ("groups", "chroms", "out"): ap.add_argument("--" + x, required=True)
ap.add_argument("--order"); ap.add_argument("--agp"); ap.add_argument("--paf"); ap.add_argument("--prefix", default="chr")
a = ap.parse_args()
size = collections.Counter(); members = collections.defaultdict(list)
for l in open(a.groups):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if re.match(r"g\d+$", s[2]): size[s[2]] += int(s[1]); members[s[2]].append(s[0])
chromOf = {l.split("\t")[0]: l.split("\t")[1].strip() for l in open(a.chroms) if l[0] != "#"}
byChrom = collections.defaultdict(list)
for g in size: byChrom[chromOf[g]].append(g)
refOf = {}; hits = collections.defaultdict(list)
if a.paf and a.agp:
    comp = collections.defaultdict(set)
    for l in open(a.agp):
        if l[0] == "#": continue
        s = l.rstrip("\n").split("\t")
        if s[4] == "W": comp[s[5]].add(s[0].split(".")[0])
    gOf = {c.split(".")[0]: g for g, cs in members.items() for c in cs}
    chromHits = collections.defaultdict(collections.Counter)
    with gzip.open(a.paf, "rt") as fh:
        for l in fh:
            s = l.split("\t", 12)
            if "tp:A:P" not in s[12] or int(s[11]) < 20: continue
            for c in comp.get(s[0], ()):
                g = gOf.get(c)
                if g: chromHits[chromOf[g]][s[5]] += int(s[3]) - int(s[2]); hits[c].append((s[5], (int(s[7]) + int(s[8])) // 2))
    for k, c in chromHits.items(): refOf[k] = c.most_common(1)[0][0]
flip = {}
if a.order and refOf:
    O = collections.defaultdict(list)
    for l in open(a.order):
        if l[0] == "#": continue
        s = l.rstrip("\n").split("\t"); O[s[0]].append(s)
    for g, v in O.items():
        main = refOf.get(chromOf.get(g)); pos = []
        for s in v:
            h = sorted(p for r, p in hits.get(s[2].split(".")[0], ()) if r == main)
            if h: pos.append(h[len(h) // 2])
        up = sum(y > x for x, y in zip(pos, pos[1:])); flip[g] = up < len(pos) - 1 - up
ranks = sorted(byChrom, key=lambda k: -sum(size[g] for g in byChrom[k]))
seen = collections.Counter()
with open(a.out, "w") as o:
    o.write("#group\tname\tflip\n")
    for i, k in enumerate(ranks, 1):
        base = refOf.get(k, f"{a.prefix}{i}")
        seen[base] += 1
        if seen[base] > 1: print(f"ALARM: {base} claimed by more than one chromosome ({k})")
        tag = base if seen[base] == 1 else f"{base}.{seen[base]}"
        for h, g in enumerate(sorted(byChrom[k], key=lambda g: -size[g]), 1):
            o.write(f"{g}\t{tag}_h{h}\t{'flip' if flip.get(g) else '.'}\n")
