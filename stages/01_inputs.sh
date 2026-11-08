#!/bin/bash
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Link inputs, index them, and check they fit each other.
source "$AS_REPO/stages/lib.sh"
mkdir -p 00_inputs
ln -sfn "$AS_UNITIGS" 00_inputs/unitigs.fasta; ln -sfn "$AS_GFA" 00_inputs/unitigs.gfa
ln -sfn "$AS_READS" 00_inputs/ont_corrected_reads.fasta
ln -sfn "$AS_HIC_R1" 00_inputs/hic_R1.fq.gz; ln -sfn "$AS_HIC_R2" 00_inputs/hic_R2.fq.gz
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
