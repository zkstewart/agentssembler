#!/bin/bash
# Apply D04 breaks, restore arms if collapsed, then gather per-chain evidence: reads, alleles, Hi-C density, repeats, telomeres, end links, depth class.
source "$AS_REPO/stages/lib.sh"
mkdir -p 11_ont_chains 12_broken 10_groups 16_validate
U=00_inputs/unitigs.fasta; F=00_inputs/unitigs.fasta.fai
[ -e 12_broken/breaks.txt ] || { if has_rows edits/chain_breaks.txt; then cp edits/chain_breaks.txt 12_broken/breaks.txt; else : > 12_broken/breaks.txt; fi; }
if [ ! -s 12_broken/chains.agp ]; then
  python3 "$S/break_chain.py" 12_broken/breaks.txt 08_final_chains/chains.agp 09_hic_map/pairs.tsv.gz 12_broken "$U"
  if collapsed && [ -s 07_chains2/restored_arms.tsv ]; then
    python3 "$S/add_restored_arms.py" 07_chains2/restored_arms.tsv "$F" 12_broken/chains.agp 12_broken/restored_chains.tsv
    agp2fa 12_broken/chains.agp "$U" 12_broken/chains.fasta
  fi
fi
[ -s 12_broken/chains.fasta.fai ] || "$AS_SAMTOOLS" faidx 12_broken/chains.fasta
map_paf 11_ont_chains/ont_vs_chains.paf.gz -x lr:hq -t "$AS_THREADS" --secondary=no 12_broken/chains.fasta 00_inputs/ont_corrected_reads.fasta
step 11_ont_chains/depth.windows.tsv python3 "$S/unitig_depth.py" 11_ont_chains/ont_vs_chains.paf.gz 12_broken/chains.fasta.fai 11_ont_chains/depth
step 10_groups/alleles.pairs.tsv python3 "$S/allele_table.py" 12_broken/chains.agp 03_unitig_align/unitigs_self.paf.gz 10_groups/alleles
step 10_groups/chain_stats.tsv python3 "$S/chain_stats.py" 12_broken/pairs.tsv.gz 12_broken/chains.fasta.fai 10_groups/alleles.pairs.tsv - 10_groups/chain_stats.tsv
step 10_groups/chain_repeat.tsv python3 "$S/chain_repeat.py" 12_broken/chains.fasta 10_groups/chain_repeat.tsv 20000
[ -s 16_validate/telomeres.chains.tsv ] || python3 "$S/telomeres.py" 12_broken/chains.fasta --motif "$AS_MOTIF" > 16_validate/telomeres.chains.tsv
step 10_groups/end_links.tsv python3 "$S/chain_end_links.py" 12_broken/chains.agp 00_inputs/unitigs.gfa 11_ont_chains/ont_vs_chains.paf.gz 12_broken/chains.fasta.fai 10_groups/end_links.tsv
step 11_ont_chains/chain_class.tsv python3 "$S/chain_class.py" 11_ont_chains/depth.windows.tsv 12_broken/chains.fasta.fai 10_groups/chain_repeat.tsv 11_ont_chains/chain_class.tsv --oneLo "$AS_ONE_LO" --sharedLo "$AS_SHARED_LO" --sharedHi "$AS_SHARED_HI" --sharedOut 10_groups/shared.candidates.txt
