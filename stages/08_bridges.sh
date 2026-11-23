#!/bin/bash
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Apply group edits, write shared copies when collapsed, find graph bridges (up to `hops` steps) and rank them. Stops for D07.
source "$AS_REPO/stages/lib.sh"
D=12_shared; mkdir -p "$D" 13_scaffold/bridges
if has_rows edits/group_edits.tsv; then python3 "$S/apply_group_edits.py" 10_groups/groups edits/group_edits.tsv 10_groups/final; else for x in tsv chroms.tsv; do cp 10_groups/groups.$x 10_groups/final.$x; done; for x in shared dosage; do [ -s 10_groups/groups.$x.tsv ] && cp 10_groups/groups.$x.tsv 10_groups/final.$x.tsv || true; done; fi
if collapsed; then
  [ -s $D/groups.tsv ] || python3 "$S/make_shared_copies.py" 10_groups/final 12_broken/chains.agp 12_broken/chains.fasta.fai 12_broken/pairs.tsv.gz 10_groups/end_links.tsv 10_groups/chain_stats.tsv 16_validate/telomeres.chains.tsv $D > $D/copies.log 2>&1
else
  cp 10_groups/final.tsv $D/groups.tsv
  ln -sfn ../12_broken/chains.agp $D/chains.agp; ln -sfn ../12_broken/chains.fasta.fai $D/chains.fasta.fai; ln -sfn ../12_broken/pairs.tsv.gz $D/pairs.tsv.gz
  ln -sfn ../10_groups/end_links.tsv $D/end_links.tsv; ln -sfn ../10_groups/chain_stats.tsv $D/chain_stats.tsv; ln -sfn ../16_validate/telomeres.chains.tsv $D/telomeres.tsv
fi
# settings.hic_allele_rescue: rescued pairs for ordering the chains rescue_groups.py placed. Kept apart from pairs.tsv.gz, which holds
# unique pairs only: stage 09 orders on both, and stage 15 validates on unique pairs alone.
if hic_rescue && [ -s 10_groups/groups.rescued.tsv ] && [ ! -s $D/pairs.rescue_order.tsv.gz ]; then
  python3 "$S/rescue/project_rescued.py" $D/groups.tsv 12_broken/pairs.rescued.tsv.gz 10_groups/alleles.blocks.tsv 10_groups/groups.rescued.tsv $D/pairs.rescue_order.tsv.gz 2> $D/projected.log
fi
step $D/graph.bridges.tsv python3 "$S/chain_bridges.py" $D/end_links.tsv $D/groups.tsv $D/graph --hops "$AS_HOPS"
if collapsed; then step $D/base.bridges.tsv python3 "$S/assign_arms.py" prepare $D/graph.bridges.tsv $D/base.bridges.tsv; B=$D/base.bridges.tsv; else B=$D/graph.bridges.tsv; fi
step 13_scaffold/bridges/choices.tsv python3 "$S/choose_bridges.py" "$B" 10_groups/end_links.tsv 12_broken/chains.fasta 11_ont_chains/ont_vs_chains.paf.gz "$AS_READS_LINK" 13_scaffold/bridges --preset "$AS_READ_PRESET" --minimap2 "$AS_MINIMAP2" --seqkit "$AS_SEQKIT" --threads "$AS_THREADS"
rm -f 13_scaffold/bridges/reads.fasta
