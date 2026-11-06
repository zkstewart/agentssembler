#!/usr/bin/env python3
'''
Write the sequences described by an AGP (W components from a contig FASTA, U/N gaps as Ns).
Usage: agp_to_fasta.py <in.agp> <contigs.fasta>[,<more.fasta>...] <out.fasta>
'''
import sys
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def read_fasta(fileName):
    seqs, name, chunks = {}, None, []
    with open(fileName) as fileIn:
        for line in fileIn:
            if line.startswith(">"):
                if name is not None: seqs[name] = "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else: chunks.append(line.strip())
    if name is not None: seqs[name] = "".join(chunks)
    return seqs

agpFile, fastaFile, outFile = sys.argv[1:4]
seqs = {}
for f in fastaFile.split(","): seqs.update(read_fasta(f))
order, parts = [], {}
for line in open(agpFile):
    if line.startswith("#"): continue
    sl = line.rstrip("\n").split("\t")
    if not sl[0] in parts:
        order.append(sl[0]); parts[sl[0]] = []
    if sl[4] in ("U", "N"):
        parts[sl[0]].append("N" * int(sl[5]))
    else:
        s = seqs[sl[5]][int(sl[6]) - 1:int(sl[7])]
        parts[sl[0]].append(s if sl[8] != "-" else s.translate(COMP)[::-1])
with open(outFile, "w") as out:
    for name in order:
        seq = "".join(parts[name])
        out.write(f">{name}\n")
        for i in range(0, len(seq), 80):
            out.write(seq[i:i+80] + "\n")
