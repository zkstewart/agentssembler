#!/usr/bin/env python3
'''
Build the final 27 (chromosome x haplotype) groups for Hi-C ordering with HapHiC sort.

Groups start from the reference-guided Hi-C phasing (hic_phase.py output). Unitigs flagged
by relocate_scan.py whose reference alignment is unreliable (aln_reliable_frac below
--minReliable) are moved to the group their Hi-C links point to; this includes previously
"redundant" unitigs with a clear Hi-C home. Flagged unitigs with reliable alignments are
left in place and listed as candidate structural differences from the reference.

Writes, under the output directory:
    groups/<chrom>h<hap>.txt  -- HapHiC group files (#Contig RECounts Length)
    split_clms/<chrom>h<hap>.clm -- the HapHiC CLM records within each group
    membership.tsv            -- every unitig's final group and how it got there
'''

import os, argparse
from collections import defaultdict
from common import natural_key, re_site_regex

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

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--phasing", required=True)
    p.add_argument("--scan", required=True, help="relocate_scan.py output")
    p.add_argument("--fasta", required=True, help="HapHiC corrected_asm.fa")
    p.add_argument("--clm", required=True, help="HapHiC paired_links.clm")
    p.add_argument("--minReliable", type=float, default=0.5)
    p.add_argument("--RE", default="GATC,GANTC", help="restriction sites, N for any base (default Arima)")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    group, how = {}, {}
    for line in open(args.phasing):
        if line.startswith("#"): continue
        sl = line.rstrip("\n").split("\t")
        if sl[3] != "assigned": continue
        regs, bps = sl[6].split(","), [int(x) for x in sl[7].split(",")]
        group[sl[0]] = f"{regs[bps.index(max(bps))].split(':')[0]}h{sl[2]}"
        how[sl[0]] = "reference+hic_phase"

    candidates = []
    for line in open(args.scan):
        if line.startswith("#"): continue
        sl = line.rstrip("\n").split("\t")
        q, status, best, reliable = sl[0], sl[1], sl[6], float(sl[11])
        if reliable < args.minReliable:
            group[q] = best
            how[q] = "hic_relocated" if status == "assigned" else f"hic_rescued_from_{status}"
        else:
            candidates.append(sl)

    # Map original unitig names onto HapHiC's corrected fragment names
    seqs = read_fasta(args.fasta)
    fragments = defaultdict(list)
    for name in seqs:
        fragments[name.split(":")[0]].append(name)
    reSite = re_site_regex(args.RE)

    os.makedirs(f"{args.out}/groups", exist_ok=True)
    os.makedirs(f"{args.out}/split_clms", exist_ok=True)
    members = defaultdict(list)
    fragGroup = {}
    for q, g in group.items():
        for frag in fragments.get(q, []):
            members[g].append(frag)
            fragGroup[frag] = g
    for g, frags in members.items():
        with open(f"{args.out}/groups/{g}.txt", "w") as fileOut:
            fileOut.write("#Contig\tRECounts\tLength\n")
            for frag in sorted(frags, key=lambda f: -len(seqs[f])):
                fileOut.write(f"{frag}\t{len(reSite.findall(seqs[frag].upper()))}\t{len(seqs[frag])}\n")

    clmOuts = { g: open(f"{args.out}/split_clms/{g}.clm", "w") for g in members }
    with open(args.clm) as fileIn:
        for line in fileIn:
            a, b = line.split("\t", 1)[0].split(" ")
            ga, gb = fragGroup.get(a[:-1]), fragGroup.get(b[:-1])
            if ga is not None and ga == gb:
                clmOuts[ga].write(line)
    for f in clmOuts.values():
        f.close()

    with open(f"{args.out}/membership.tsv", "w") as fileOut:
        fileOut.write("#contig\tgroup\thow\n")
        for q in sorted(group):
            fileOut.write(f"{q}\t{group[q]}\t{how[q]}\n")
    with open(f"{args.out}/sv_candidates.tsv", "w") as fileOut:
        fileOut.write("#contig\tstatus\tlength\tref_chrom\tlabel\tlinks\tbest_group\tbest_enrich\tbest_own_chrom_group\town_enrich\tratio\taln_reliable_frac\n")
        for sl in candidates:
            fileOut.write("\t".join(sl) + "\n")

    for g in sorted(members, key=natural_key):
        print(g, len(members[g]), sum(len(seqs[f]) for f in members[g]))

if __name__ == "__main__":
    main()
