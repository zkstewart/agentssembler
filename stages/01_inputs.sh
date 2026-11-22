#!/bin/bash
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Link inputs, index them, and check they fit each other. With compleasm set, genes on the unitigs (01_genes/, about 1 h).
source "$AS_REPO/stages/lib.sh"
mkdir -p 00_inputs
ln -sfn "$AS_UNITIGS" 00_inputs/unitigs.fasta; ln -sfn "$AS_GFA" 00_inputs/unitigs.gfa
ln -sfn "$AS_READS" "$AS_READS_LINK"
ln -sfn "$AS_HIC_R1" "$AS_HIC_R1_LINK"; ln -sfn "$AS_HIC_R2" "$AS_HIC_R2_LINK"
has_ref && ln -sfn "$AS_REFERENCE" 00_inputs/reference.fasta
[ -s 00_inputs/unitigs.fasta.fai ] || "$AS_SAMTOOLS" faidx 00_inputs/unitigs.fasta -o 00_inputs/unitigs.fasta.fai
if has_ref; then [ -s 00_inputs/reference.fasta.fai ] || "$AS_SAMTOOLS" faidx 00_inputs/reference.fasta -o 00_inputs/reference.fasta.fai; fi
"$AS_REPO/stages/check_tools.sh"
# every S line of the GFA should be in the FASTA, with the same length where one is given
python3 - <<'PY'
L = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open("00_inputs/unitigs.fasta.fai")}
n = bad = 0
for l in open("00_inputs/unitigs.gfa"):
    if l[0] != "S": continue
    s = l.split("\t"); n += 1
    if s[1] not in L: bad += 1
print(f"GFA S lines {n}, FASTA records {len(L)}, GFA records not in FASTA {bad}")
if bad or n != len(L): raise SystemExit("unitigs.gfa does not match unitigs.fasta")
PY
# Gene copies on the unitigs: what the assembly could hold. Later stages lift these hits onto chains to say how many genes each
# decision gains or loses (D04/D05 rescue, D16 tails, D13 unloc) without running compleasm again. Here, so that a wrong
# compleasm path or lineage stops the run before hours of mapping. About 1 h at 16 threads.
if [ -n "${AS_COMPLEASM:-}" ] && [ ! -s 01_genes/summary.txt ]; then
  "$AS_COMPLEASM" run -a 00_inputs/unitigs.fasta -o 01_genes -l "$AS_LINEAGE" -L "$AS_BUSCO_LIB" -t "$AS_THREADS" > 01_genes.log 2>&1 || { echo "ERROR: compleasm on the unitigs failed: see 01_genes.log" >&2; exit 1; }
fi
