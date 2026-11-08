#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Call telomeres at the ends of each sequence by counting the telomeric repeat unit (either
strand) in the terminal --window bp. An end with >= --minCopies copies is called telomeric.

Output columns: seq, length, start_repeats, end_repeats, call, where call is two characters
(start, end), each 'T' (telomeric) or '-'. For unitigs, 'T-' means a telomere at the
unitig's start and '-T' at its end.

The default motif is the Arabidopsis-type plant repeat (CCCTAAA at chromosome starts,
TTTAGGG at ends); use CCCTAA for the vertebrate/TTAGGG-type repeat, for example.
'''
import argparse

def revcomp(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fasta")
    p.add_argument("--motif", default="CCCTAAA", help="telomeric repeat unit (default: plant CCCTAAA)")
    p.add_argument("--window", type=int, default=20000)
    p.add_argument("--minCopies", type=int, default=50)
    args = p.parse_args()
    fwd, rev = args.motif.upper(), revcomp(args.motif.upper())

    def ends(name, seq):
        s = seq.upper()
        head, tail = s[:args.window], s[-args.window:]
        h = head.count(fwd) + head.count(rev)
        t = tail.count(rev) + tail.count(fwd)
        call = ("T" if h >= args.minCopies else "-") + ("T" if t >= args.minCopies else "-")
        print(f"{name}\t{len(seq)}\t{h}\t{t}\t{call}")

    print("#seq\tlength\tstart_repeats\tend_repeats\tcall")
    name, chunks = None, []
    for line in open(args.fasta):
        if line.startswith(">"):
            if name: ends(name, "".join(chunks))
            name, chunks = line[1:].split()[0], []
        else:
            chunks.append(line.strip())
    if name: ends(name, "".join(chunks))

if __name__ == "__main__":
    main()
