#!/usr/bin/env python3
'''
compleasm_copies.py -- ploidy-aware gene completeness from one pooled compleasm run per assembly,
using the labels of hic_crosstalk.py so every tool is judged the same way (no haplotype labels).

Run compleasm on all of an assembly's own sequences (the hic_crosstalk mapping target without its
decoys, so sequence names are the labels' target IDs):
    awk '/^>/{keep = ($1 !~ /^>decoy_/)} keep' xt/<tool>/prepare/mapping_target.fa > <tool>.fa
    compleasm run -a <tool>.fa -o compleasm/<tool> -l eudicotyledons_odb12 -t 12

Then, per gene, count its complete copies (Single/Duplicated rows) and compare with --ploidy:
    C   complete: exactly ploidy copies        <- the expected state of a phased assembly
        [chromosomal: every copy on a backbone (chromosome-scale scaffold), whatever its
                      chromosome label, so translocations are not penalised
         unscaffolded: at least one copy on a fragment or unplaced sequence]
    D-  fewer copies (1 .. ploidy-1)           collapsed or lost haplotype copy
    D+  more copies than ploidy                false duplication / unpurged haplotig
    F   fragmented only,  M  missing
C + D- + D+ + F + M = 100%, and chromosomal + unscaffolded = C.

Usage
    compleasm_copies.py --ploidy 2 --run agentssembler compleasm/agentssembler/eudicotyledons_odb12/full_table.tsv \\
        xt/agentssembler/prepare/labels.v2.tsv --run yahs <full_table.tsv> <labels.v2.tsv> [--detail out.tsv]
'''
import argparse, collections, sys

def read_labels(path):
    out = {}
    with open(path) as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            out[r["target_id"]] = r
    return out

def tally(table, labels, ploidy):
    copies, frag, genes = collections.defaultdict(list), set(), set()
    with open(table) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            g = r["Gene"]
            genes.add(g)
            if r["Status"] in ("Single", "Duplicated"):
                lab = labels.get(r["Sequence"])
                if lab is None:
                    sys.exit(f"{r['Sequence']} (in {table}) is not a target_id of the labels file; "
                             f"was compleasm run on the mapping target without decoys?")
                copies[g].append((r["Sequence"], lab["class"], lab["chrom"]))
            elif r["Status"] in ("Fragmented", "Interspersed"):
                frag.add(g)
    cat, detail = collections.Counter(), []
    for g in sorted(genes):
        c = copies.get(g, [])
        n = len(c)
        if n == ploidy:
            k = "C"
            if all(cls == "backbone" for s, cls, ch in c):
                cat["C_chromosomal"] += 1
            else:
                cat["C_unscaffolded"] += 1
        elif 0 < n < ploidy:
            k = "D-"
        elif n > ploidy:
            k = "D+"
        else:
            k = "F" if g in frag else "M"
        cat[k] += 1
        detail.append((g, k, n, ",".join(f"{s}:{cls}:{ch}" for s, cls, ch in c) or "."))
    return len(genes), cat, detail

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ploidy", type=int, required=True)
    ap.add_argument("--run", nargs=3, action="append", required=True,
                    metavar=("NAME", "FULL_TABLE", "LABELS"), help="repeat once per assembly")
    ap.add_argument("--detail", help="per-gene TSV of every assembly")
    a = ap.parse_args()
    pct = lambda x, n: f"{100 * x / n:.1f}%" if n else "nan"
    det = open(a.detail, "w") if a.detail else None
    if det:
        det.write("#assembly\tgene\tcategory\tcomplete_copies\tcopies(sequence:class:chrom)\n")
    print("#assembly\tC\tC_chromosomal\tC_unscaffolded\tD-\tD+\tF\tM\tn\tline")
    for name, table, labels in a.run:
        n, c, detail = tally(table, read_labels(labels), a.ploidy)
        vals = [pct(c[k], n) for k in ("C", "C_chromosomal", "C_unscaffolded", "D-", "D+", "F", "M")]
        line = (f"C:{vals[0]}[chromosomal:{vals[1]},unscaffolded:{vals[2]}],D-:{vals[3]},D+:{vals[4]},"
                f"F:{vals[5]},M:{vals[6]},n:{n}")
        print("\t".join([name] + vals + [str(n), line]))
        if det:
            for row in detail:
                det.write("\t".join(map(str, (name,) + row)) + "\n")
    if det:
        det.close()

if __name__ == "__main__":
    main()
