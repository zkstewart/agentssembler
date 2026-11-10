#!/bin/bash
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Non-redundant unitigs, then chains: first round, redundancy and dovetails, join verification, rule-based breaks. Writes 04_copynumber/copy_number.tsv, 05_nonredundant, 06_chains, 07_chains2, 08_final_chains. Stops for D03.
source "$AS_REPO/stages/lib.sh"
mkdir -p 05_nonredundant 06_chains 07_chains2 08_final_chains
U=00_inputs/unitigs.fasta; F=00_inputs/unitigs.fasta.fai; G=00_inputs/unitigs.gfa; R=00_inputs/ont_corrected_reads.fasta
step 04_copynumber/copy_number.tsv python3 "$S/copy_number.py" 04_copynumber/class_depth.tsv 03_unitig_align/unitigs_self.paf.gz 04_copynumber/copy_number --fai "$F" --hap "$AS_HAP" --repeat "$AS_REPEAT"
[ -s 05_nonredundant/kept.txt ] || awk -F'\t' '!/^#/ && $10=="keep"{print $1}' 04_copynumber/copy_number.tsv > 05_nonredundant/kept.txt
step 06_chains/chains.agp python3 "$S/build_chains.py" --keep 05_nonredundant/kept.txt --gfa "$G" --fai "$F" --out 06_chains/chains
[ -s 06_chains/chains.fasta ] || agp2fa 06_chains/chains.agp "$U" 06_chains/chains.fasta
[ -s 06_chains/chains.fasta.fai ] || "$AS_SAMTOOLS" faidx 06_chains/chains.fasta
map_paf 06_chains/chains_self.paf.gz -x asm5 -c -t "$AS_THREADS" -K 200M -D --dual=no --secondary=yes -N 20 -p 0.3 06_chains/chains.fasta 06_chains/chains.fasta
map_paf 06_chains/ont_vs_chains.paf.gz -x lr:hq -t "$AS_THREADS" --secondary=no 06_chains/chains.fasta "$R"
step 06_chains/depth.windows.tsv python3 "$S/unitig_depth.py" 06_chains/ont_vs_chains.paf.gz 06_chains/chains.fasta.fai 06_chains/depth
step 07_chains2/kept2.txt python3 "$S/round2.py" 06_chains/chains.agp 06_chains/chains.fasta.fai 06_chains/depth.windows.tsv 06_chains/chains_self.paf.gz 05_nonredundant/kept.txt 07_chains2 --hap "$AS_HAP" --repeat "$AS_REPEAT"
[ -s 07_chains2/links.gfa ] || cat <(grep '^L' "$G") 07_chains2/dovetails.gfa > 07_chains2/links.gfa
step 07_chains2/chains.agp python3 "$S/build_chains.py" --keep 07_chains2/kept2.txt --gfa 07_chains2/links.gfa --fai "$F" --out 07_chains2/chains
[ -s 07_chains2/chains.fasta ] || agp2fa 07_chains2/chains.agp "$U" 07_chains2/chains.fasta
[ -s 07_chains2/chains.fasta.fai ] || "$AS_SAMTOOLS" faidx 07_chains2/chains.fasta
map_paf 07_chains2/ont_vs_chains.paf.gz -x lr:hq -c -t "$AS_THREADS" --secondary=no 07_chains2/chains.fasta "$R"
step 07_chains2/joins.verify.tsv python3 "$S/verify_joins.py" 07_chains2/ont_vs_chains.paf.gz 07_chains2/chains.joins.tsv 07_chains2/joins.verify.tsv
if collapsed; then step 07_chains2/restored_arms.tsv python3 "$S/restore_arms.py" 04_copynumber/copy_number.tsv 07_chains2/kept2.txt 07_chains2/repeat_chains.txt "$G" 07_chains2/restored_arms.tsv --lo "$AS_RESTORE_LO" --hi "$AS_RESTORE_HI"; fi
# rule: under 3 spanning reads breaks the join. Further breaks go in edits/join_breaks.txt.
[ -e 08_final_chains/breaks.txt ] || { awk -F'\t' '!/^#/ && $6<3 {print $3, ($4 ~ /\+$/ ? substr($4,1,length($4)-1)"-" : substr($4,1,length($4)-1)"+")}' 07_chains2/joins.verify.tsv > 08_final_chains/breaks.txt; has_rows edits/join_breaks.txt && grep -v '^#' edits/join_breaks.txt | cut -f1 >> 08_final_chains/breaks.txt || true; }
step 08_final_chains/chains.agp python3 "$S/build_chains.py" --keep 07_chains2/kept2.txt --gfa 07_chains2/links.gfa --fai "$F" --out 08_final_chains/chains --breaks 08_final_chains/breaks.txt
[ -s 08_final_chains/chains.fasta ] || agp2fa 08_final_chains/chains.agp "$U" 08_final_chains/chains.fasta
[ -s 08_final_chains/chains.fasta.fai ] || "$AS_SAMTOOLS" faidx 08_final_chains/chains.fasta
