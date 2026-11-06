#!/bin/bash
# Apply group edits, write shared copies when collapsed, find graph bridges (8 hops) and rank them. Stops for D07.
source "$AS_REPO/stages/lib.sh"
D=12_shared; mkdir -p "$D" 13_scaffold/bridges
if has_rows edits/group_edits.tsv; then python3 "$S/apply_group_edits.py" 10_groups/groups edits/group_edits.tsv 10_groups/final; else for x in tsv chroms.tsv; do cp 10_groups/groups.$x 10_groups/final.$x; done; [ -s 10_groups/groups.shared.tsv ] && cp 10_groups/groups.shared.tsv 10_groups/final.shared.tsv || true; fi
if collapsed; then
  [ -s $D/groups.tsv ] || python3 "$S/make_shared_copies.py" 10_groups/final 12_broken/chains.agp 12_broken/chains.fasta.fai 12_broken/pairs.tsv.gz 10_groups/end_links.tsv 10_groups/chain_stats.tsv 16_validate/telomeres.chains.tsv $D > $D/copies.log 2>&1
else
  cp 10_groups/final.tsv $D/groups.tsv
  ln -sfn ../12_broken/chains.agp $D/chains.agp; ln -sfn ../12_broken/chains.fasta.fai $D/chains.fasta.fai; ln -sfn ../12_broken/pairs.tsv.gz $D/pairs.tsv.gz
  ln -sfn ../10_groups/end_links.tsv $D/end_links.tsv; ln -sfn ../10_groups/chain_stats.tsv $D/chain_stats.tsv; ln -sfn ../16_validate/telomeres.chains.tsv $D/telomeres.tsv
fi
step $D/graph.bridges.tsv python3 "$S/chain_bridges.py" $D/end_links.tsv $D/groups.tsv $D/graph --hops "$AS_HOPS"
if collapsed; then step $D/base.bridges.tsv python3 "$S/assign_arms.py" prepare $D/graph.bridges.tsv $D/base.bridges.tsv; B=$D/base.bridges.tsv; else B=$D/graph.bridges.tsv; fi
step 13_scaffold/bridges/choices.tsv python3 "$S/choose_bridges.py" "$B" 10_groups/end_links.tsv 12_broken/chains.fasta 11_ont_chains/ont_vs_chains.paf.gz 00_inputs/ont_corrected_reads.fasta 13_scaffold/bridges
rm -f 13_scaffold/bridges/reads.fasta
