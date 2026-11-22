#!/bin/bash
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Hi-C to the final chains (09_hic_map), and the chimera table for D04 (10_groups/chain_breaks.tsv). This is the only Hi-C mapping: there is no 01_ directory. About 1-2 h. Filters: MAPQ >= 1, NM < 3.
source "$AS_REPO/stages/lib.sh"
mkdir -p 09_hic_map 10_groups
cd 09_hic_map
ln -sfn ../08_final_chains/chains.fasta chains.fasta; ln -sfn ../08_final_chains/chains.fasta.fai chains.fasta.fai
# Decoys: unitigs of the high-copy chains that stage 04 set aside (07_chains2/repeat_chains.txt), organelles among them.
# Without them, chloroplast Hi-C reads land on nuclear plastid insertions as if unique: on 13Q010, 21,394 pairs at one
# point joined a chr9 and a chr5 chain and changed the grouping. Pairs touching a decoy are dropped.
if [ ! -s pairs.tsv.gz ]; then
  guard
  [ -s targets.fasta ] || {
    awk '!/^#/{n=split($3,u,","); for(i=1;i<=n;i++){x=u[i]; sub(/[+-]$/,"",x); print x}}' ../07_chains2/repeat_chains.txt 2>/dev/null | sort -u > decoys.all.txt
    awk -F'\t' '!/^#/ && $5=="W"{print $6}' ../08_final_chains/chains.agp | sort -u | comm -23 decoys.all.txt - > decoys.txt
    cp chains.fasta targets.fasta.part
    if [ -s decoys.txt ]; then "$AS_SAMTOOLS" faidx ../00_inputs/unitigs.fasta -r decoys.txt | sed 's/^>/>decoy_/' >> targets.fasta.part; fi
    mv targets.fasta.part targets.fasta; rm -f decoys.all.txt
  }
  "$AS_BWA" index targets.fasta > bwa_index.log 2>&1
  "$AS_BWA" mem -5SP -t "$AS_THREADS" targets.fasta "../$AS_HIC_R1_LINK" "../$AS_HIC_R2_LINK" 2> bwa.log | "$AS_SAMBLASTER" 2> samblaster.log | "$AS_SAMTOOLS" view - -@ 4 -S -h -b -F 3340 -o HiC.bam.part
  mv HiC.bam.part HiC.bam
  "$AS_FILTER_BAM" HiC.bam 1 --nm 3 --threads "$AS_THREADS" | "$AS_SAMTOOLS" view - -b -@ 4 -o HiC.filtered.bam.part
  mv HiC.filtered.bam.part HiC.filtered.bam
  "$AS_SAMTOOLS" view -@ 4 -f 64 HiC.filtered.bam | awk -F'\t' -v OFS='\t' '$3 !~ /^decoy_/ && $7 !~ /^decoy_/ {print $3, $4, ($7 == "=" ? $3 : $7), $8}' | gzip -1 > pairs.tsv.gz.part
  mv pairs.tsv.gz.part pairs.tsv.gz
fi
# Pairs with an allele-ambiguous end (MAPQ 0, every other hit on an allele partner), for settings.hic_allele_rescue. Made whatever the
# setting, at about 10 min for 50 M pairs, so that D04 can see whether grouping will need them (rescue_signal.txt) and turning the
# rescue on later costs no Hi-C pass. Alleles of these chains first.
if [ ! -s pairs.rescued.tsv.gz ]; then
  [ -s alleles08.pairs.tsv ] || python3 "$S/allele_table.py" ../08_final_chains/chains.agp ../03_unitig_align/unitigs_self.paf.gz alleles08
  "$AS_SAMTOOLS" view -@ 2 HiC.bam | python3 "$S/rescue/hic_rescue_pairs.py" alleles08.pairs.tsv pairs.rescued.tsv.gz.part 2> rescue.log
  mv pairs.rescued.tsv.gz.part pairs.rescued.tsv.gz
fi
cd ..
G=(); has_genes && G=(--genes "$GENES" --agp 08_final_chains/chains.agp)
python3 "$S/rescue/rescue_signal.py" 09_hic_map/pairs.tsv.gz 09_hic_map/pairs.rescued.tsv.gz 08_final_chains/chains.fasta.fai 09_hic_map/alleles08.pairs.tsv 09_hic_map/rescue_signal.tsv "${G[@]}" > 09_hic_map/rescue_signal.txt
step 10_groups/chain_breaks.tsv python3 "$S/chain_breaks.py" 09_hic_map/pairs.tsv.gz 08_final_chains/chains.agp 10_groups/chain_breaks.tsv --dist "$AS_BREAK_DIST" --flank "$AS_BREAK_FLANK" --minSide "$AS_BREAK_MIN_SIDE"
