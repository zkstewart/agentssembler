#!/bin/bash
# Hi-C to the final chains, and the chimera table for D04. About 1-2 h. Filters: MAPQ >= 1, NM < 3.
source "$AS_REPO/stages/lib.sh"
mkdir -p 09_hic_map 10_groups
cd 09_hic_map
ln -sfn ../08_final_chains/chains.fasta chains.fasta; ln -sfn ../08_final_chains/chains.fasta.fai chains.fasta.fai
if [ ! -s pairs.tsv.gz ]; then
  guard
  "$AS_BWA" index chains.fasta > bwa_index.log 2>&1
  "$AS_BWA" mem -5SP -t "$AS_THREADS" chains.fasta ../00_inputs/hic_R1.fq.gz ../00_inputs/hic_R2.fq.gz 2> bwa.log | "$AS_SAMBLASTER" 2> samblaster.log | "$AS_SAMTOOLS" view - -@ 4 -S -h -b -F 3340 -o HiC.bam.part
  mv HiC.bam.part HiC.bam
  "$AS_FILTER_BAM" HiC.bam 1 --nm 3 --threads "$AS_THREADS" | "$AS_SAMTOOLS" view - -b -@ 4 -o HiC.filtered.bam.part
  mv HiC.filtered.bam.part HiC.filtered.bam
  "$AS_SAMTOOLS" view -@ 4 -f 64 HiC.filtered.bam | awk -F'\t' -v OFS='\t' '{print $3, $4, ($7 == "=" ? $3 : $7), $8}' | gzip -1 > pairs.tsv.gz.part
  mv pairs.tsv.gz.part pairs.tsv.gz
fi
cd ..
step 10_groups/chain_breaks.tsv python3 "$S/chain_breaks.py" 09_hic_map/pairs.tsv.gz 08_final_chains/chains.agp 10_groups/chain_breaks.tsv
