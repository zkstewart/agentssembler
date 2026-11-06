#!/usr/bin/env python3
'''
Infer the Hi-C restriction enzyme(s) from ligation junctions in read 1.

When a Hi-C fragment is cut, filled in and ligated, the junction duplicates part of the
site, leaving a characteristic motif in reads that cross it (DpnII GATC -> GATCGATC). This
counts each enzyme's junction motifs in the first --reads reads and compares the share of
reads containing them with the share expected by chance (from base composition and read
length, or better, from how often the motif occurs in the genome itself with --genome). An
enrichment well above 1 indicates that enzyme. Omni-C/Micro-C (nuclease-based) libraries
show no enriched junction.

Usage: detect_enzyme.py reads_R1.fq.gz --genome assembly.fasta

The suggested --RE value (for HapHiC and make_groups.py) is printed on the last line.
'''
import re, gzip, argparse

# name: (RE sites for HapHiC, junction motifs; N = any base)
ENZYMES = {
    "DpnII/MboI (GATC)":     ("GATC",       ["GATCGATC"]),
    "Arima (GATC + GANTC)":  ("GATC,GANTC", ["GANTGATC", "GATCANTC", "GANTANTC"]), # plus GATCGATC
    "HindIII (AAGCTT)":      ("AAGCTT",     ["AAGCTAGCTT"]),
    "NlaIII (CATG)":         ("CATG",       ["CATGCATG"]),
    "DdeI (CTNAG)":          ("CTNAG",      ["CTNATNAG"]),
    "MseI (TTAA)":           ("TTAA",       ["TTATAA"]),
}

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fastq", help="read 1 FASTQ (gzipped or not)")
    p.add_argument("--reads", type=int, default=500000)
    p.add_argument("--genome", default=None, help="assembly FASTA: motif rates in the genome give the expected background")
    p.add_argument("--genomeBp", type=int, default=50000000, help="bp of --genome to scan (default 50 Mb)")
    args = p.parse_args()

    pats = { name: [ re.compile(m.replace("N", "[ACGT]")) for m in motifs ] for name, (_, motifs) in ENZYMES.items() }
    hits = { name: 0 for name in ENZYMES }
    base = { b: 0 for b in "ACGT" }
    lengths = 0
    opener = gzip.open if args.fastq.endswith(".gz") else open
    n = 0
    with opener(args.fastq, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 != 1: continue
            n += 1; read = line.strip(); lengths += len(read)
            for b in base: base[b] += read.count(b)
            for name, ps in pats.items():
                if any(pt.search(read) for pt in ps): hits[name] += 1
            if n >= args.reads: break

    # Expected share of reads containing a motif by chance: from the motif's rate per bp in the
    # genome (both strands; the motifs are palindromes or come as reverse-complement pairs), or
    # else from base composition and read length
    total = sum(base.values()); freq = { b: c / total for b, c in base.items() }
    readLen = lengths / n
    genomeRate = None
    if args.genome:
        genomeRate = { name: 0 for name in ENZYMES }; scanned = 0
        with (gzip.open if args.genome.endswith(".gz") else open)(args.genome, "rt") as f:
            chunk = []
            for line in f:
                if line.startswith(">"):
                    seq = "".join(chunk).upper(); chunk = []
                    for name, ps in pats.items():
                        genomeRate[name] += sum(len(pt.findall(seq)) for pt in ps)
                    scanned += len(seq)
                    if scanned >= args.genomeBp: break
                else:
                    chunk.append(line.strip())
        genomeRate = { k: v / scanned for k, v in genomeRate.items() }
    def expected(motifs, name):
        if genomeRate is not None:
            return 1 - (1 - genomeRate[name]) ** readLen
        miss = 1.0
        for m in motifs:
            prob = 1.0
            for ch in m: prob *= 1.0 if ch == "N" else freq[ch]
            miss *= (1 - prob) ** max(1, readLen - len(m) + 1)
        return 1 - miss

    print(f"# {n} reads scanned (mean length {readLen:.0f} bp); expected share from "
          f"{'the genome' if genomeRate is not None else 'base composition (use --genome for a better background)'}")
    print(f"{'enzyme':24s} {'observed':>9s} {'expected':>9s} {'enrichment':>10s}")
    enrich = {}
    for name, (_, motifs) in ENZYMES.items():
        exp = expected(motifs, name); obs = hits[name] / n
        enrich[name] = obs / exp if exp > 0 else 0
        print(f"{name:24s} {100*obs:8.2f}% {100*exp:8.2f}% {enrich[name]:9.1f}x")
    # Arima is identified by enriched GANTC junctions; otherwise take the most enriched enzyme
    # (its N-containing motifs have a high chance rate, so a lower threshold applies)
    if enrich["Arima (GATC + GANTC)"] >= 2 and enrich["DpnII/MboI (GATC)"] >= 3:
        best = "Arima (GATC + GANTC)"
    else:
        best = max((k for k in ENZYMES if k != "Arima (GATC + GANTC)"), key=lambda k: enrich[k])
    if enrich[best] < 3:
        print("# No junction motif is enriched >= 3x: possibly a nuclease-based (Omni-C/Micro-C) library;")
        print("# HapHiC's --RE then has no real meaning (it only affects normalisation).")
    print(f"RE={ENZYMES[best][0]}")

if __name__ == "__main__":
    main()
