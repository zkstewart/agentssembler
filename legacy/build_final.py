#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Build the final 27 pseudochromosomes (chr1-9 x h1-3) by merging two orderings:

- Skeleton: RagTag's reference-guided order and orientation of each group's unitigs
  (partition K of the reference-guided Hi-C phasing is haplotype hK).
- Insertions: group members RagTag did not place on that chromosome (Hi-C relocated or
  rescued unitigs, or ones RagTag left unplaced) are placed by HapHiC's Hi-C tour of the
  group: each run of such unitigs is attached beside its nearest skeleton neighbour in the
  tour, flipped as needed to follow the skeleton's orientation of that neighbour.

Writes <out>.fasta and <out>.agp (AGP 2.1; evidence 'align_genus' for skeleton joins and
'proximity_ligation' for Hi-C insertions), plus <out>.unplaced.fasta holding the
non-redundant unitigs that are in no group.
'''

import argparse
from collections import defaultdict
from common import natural_key

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def read_fasta(fileName):
    seqs, name, chunks = {}, None, []
    with open(fileName) as fileIn:
        for line in fileIn:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks)
    return seqs

def read_ragtag_agp(fileName):
    order = defaultdict(list) # scaffold -> [(unitig, orient)]
    for line in open(fileName):
        if line.startswith("#"): continue
        sl = line.rstrip("\n").split("\t")
        if sl[4] == "W":
            order[sl[0]].append((sl[5], sl[8]))
    return order

def read_tour(fileName):
    last = None
    for line in open(fileName):
        if line.strip() and not line.startswith(">"):
            last = line.split()
    return [ (x[:-1], x[-1]) for x in last ]

def flip(o):
    return "-" if o == "+" else "+"

def merge(skeleton, tour, members):
    skelSet = { u for u, _ in skeleton }
    skelOrient = dict(skeleton)
    tour = [ (u, o) for u, o in tour if u in members ]
    toursSet = { u for u, _ in tour }
    # Members missing from the tour (e.g. split by HapHiC's correction) go at the end
    tail = [ (u, "+") for u in members if not u in skelSet and not u in toursSet ]

    # Group non-skeleton runs with the anchor before them (or after, for a leading run)
    after, before = defaultdict(list), defaultdict(list)
    run, prevAnchor = [], None
    for u, o in tour:
        if u in skelSet:
            if run:
                if prevAnchor is None:
                    before[(u, o)].extend(run)
                else:
                    after[prevAnchor].extend(run)
            run, prevAnchor = [], (u, o)
        else:
            run.append((u, o))
    if run:
        if prevAnchor is None:
            tail = run + tail # no skeleton in the tour at all
        else:
            after[prevAnchor].extend(run)

    result = []
    tourOrient = { u: o for u, o in tour }
    for u, o in skeleton:
        anchor = (u, tourOrient.get(u))
        same = tourOrient.get(u) == o
        pre, post = before.get(anchor, []), after.get(anchor, [])
        if not same:
            # Anchor is reversed relative to the tour: mirror the attached runs
            pre, post = [ (x, flip(y)) for x, y in reversed(post) ], [ (x, flip(y)) for x, y in reversed(pre) ]
        result.extend([ (x, y, "hic") for x, y in pre ])
        result.append((u, o, "ref"))
        result.extend([ (x, y, "hic") for x, y in post ])
    result.extend([ (x, y, "hic") for x, y in tail ])
    return result

def cap_telomeres(merged, telomeric):
    '''
    Moves Hi-C-placed unitigs carrying a terminal telomere (telomeric[u] is 'T-' for a
    telomere at the unitig's start, '-T' at its end) to a pseudochromosome end that lacks a
    telomere, oriented telomere-outward. If both ends lack one, the end nearer the unitig's
    current position is used. Reference-placed unitigs are never moved.
    Returns (merged, moves) where moves lists (unitig, 'start'|'end').
    '''
    def outward(u, o, end):
        # True if u in orientation o presents its telomere at the given pseudochromosome end
        atStart = (telomeric[u] == "T-") == (o == "+")
        return atStart if end == "start" else not atStart
    def capped(items, end):
        u, o, _ = items[0] if end == "start" else items[-1]
        return u in telomeric and outward(u, o, end)
    moves = []
    for u in [ x for x, _, ev in merged if ev == "hic" and x in telomeric and telomeric[x] in ("T-", "-T") ]:
        index = [ x for x, _, _ in merged ].index(u)
        rest = merged[:index] + merged[index+1:]
        free = [ end for end in ("start", "end") if not capped(rest, end) ] if rest else ["start"]
        if free == []:
            continue
        end = free[0] if len(free) == 1 else ("start" if index < len(merged) / 2 else "end")
        o = "+" if outward(u, "+", end) else "-"
        merged = [(u, o, "hic")] + rest if end == "start" else rest + [(u, o, "hic")]
        moves.append((u, end))
    return merged, moves

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--membership", required=True)
    p.add_argument("--ragtag", nargs=3, required=True, help="RagTag AGPs of partitions 1, 2 and 3")
    p.add_argument("--tours", required=True, help="HapHiC final_tours directory")
    p.add_argument("--fasta", required=True, help="unitig FASTA")
    p.add_argument("--phasing", required=True, help="phasing TSV, to find redundant unitigs")
    p.add_argument("--telomeric", default=None, help="telomeres.py output for the unitigs; enables end-capping")
    p.add_argument("--gap", type=int, default=100)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    members = defaultdict(set)
    for line in open(args.membership):
        if line.startswith("#"): continue
        q, g, _ = line.rstrip("\n").split("\t")
        members[g].add(q)
    ragtag = { str(k + 1): read_ragtag_agp(f) for k, f in enumerate(args.ragtag) }
    telomeric = {}
    if args.telomeric:
        for line in open(args.telomeric):
            if line.startswith("#"): continue
            sl = line.rstrip("\n").split("\t")
            telomeric[sl[0]] = sl[4]
    capLog = open(args.out + ".telomere_moves.tsv", "w")
    capLog.write("#pseudochromosome\tunitig\tmoved_to\n")
    seqs = read_fasta(args.fasta)

    agp = open(args.out + ".agp", "w")
    agp.write("##agp-version 2.1\n")
    fasta = open(args.out + ".fasta", "w")
    used = set()
    for g in sorted(members, key=natural_key):
        chrom, hap = g.rsplit("h", 1)
        skeleton = [ (u, o) for u, o in ragtag[hap].get(f"{chrom}_RagTag", []) if u in members[g] ]
        tour = read_tour(f"{args.tours}/{g}.tour")
        tour = [ (u.split(":")[0], o) for u, o in tour ] # collapse HapHiC-corrected fragments
        seen, dedup = set(), []
        for u, o in tour:
            if not u in seen:
                seen.add(u); dedup.append((u, o))
        merged = merge(skeleton, dedup, members[g])
        if telomeric:
            merged, moves = cap_telomeres(merged, telomeric)
            for u, end in moves:
                capLog.write(f"{chrom}_h{hap}\t{u}\t{end}\n")

        name = f"{chrom}_h{hap}"
        pos, part, chunks = 1, 1, []
        for k, (u, o, ev) in enumerate(merged):
            if k > 0:
                evidence = "proximity_ligation" if ev == "hic" or merged[k-1][2] == "hic" else "align_genus"
                agp.write(f"{name}\t{pos}\t{pos+args.gap-1}\t{part}\tU\t{args.gap}\tscaffold\tyes\t{evidence}\n")
                chunks.append("N" * args.gap); pos += args.gap; part += 1
            s = seqs[u] if o == "+" else seqs[u].translate(COMP)[::-1]
            agp.write(f"{name}\t{pos}\t{pos+len(s)-1}\t{part}\tW\t{u}\t1\t{len(s)}\t{o}\n")
            chunks.append(s); pos += len(s); part += 1
            used.add(u)
        seq = "".join(chunks)
        fasta.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            fasta.write(seq[i:i+80] + "\n")
        print(f"{name}\t{len(seq)}\t{len(merged)} unitigs\t{sum(1 for x in merged if x[2] == 'hic')} Hi-C placed")
    agp.close(); fasta.close(); capLog.close()

    redundant = set()
    for line in open(args.phasing):
        sl = line.split("\t")
        if len(sl) > 3 and sl[3] == "redundant":
            redundant.add(sl[0])
    with open(args.out + ".unplaced.fasta", "w") as fileOut:
        for u, s in seqs.items():
            if not u in used and not u in redundant:
                fileOut.write(f">{u}\n")
                for i in range(0, len(s), 80):
                    fileOut.write(s[i:i+80] + "\n")

if __name__ == "__main__":
    main()
