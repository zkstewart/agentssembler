#!/bin/bash
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Names, arms at bridges, a pre-build, one-sided tails and their tests, and the lost-route fork scan for D16's tangle check. Writes into 12_shared and 13_scaffold. Stops for D17, D15 (both only if collapsed) and D16: nothing is placed without edits/accepted_tails.tsv.
source "$AS_REPO/stages/lib.sh"
D=12_shared; T=13_scaffold; mkdir -p $T $D/pre
O=$D/order.tsv; GRP=$D/groups.tsv; [ -s $D/groups.edited.tsv ] && GRP=$D/groups.edited.tsv
# Hand names must end in _h<k> and be unique: stages 12-15 find each scaffold's haplotype and its homologs by that suffix
if has_rows edits/names.tsv; then
  awk -F'\t' '$1!~/^#/ && NF{ if ($2 !~ /_h[0-9]+$/) { print "ERROR: edits/names.tsv: " $2 " does not end in _h<k> (a haplotype number)"; e=1 } if (n[$2]++) { print "ERROR: edits/names.tsv: " $2 " is used twice"; e=1 } } END{exit e}' edits/names.tsv >&2
  cp edits/names.tsv $T/names.tsv
elif [ ! -s $T/names.tsv ]; then
  A=(--groups "$GRP" --chroms 10_groups/final.chroms.tsv --out $T/names.tsv --order "$O")
  has_ref && A+=(--agp 12_broken/chains.agp --paf 03_unitig_align/unitigs_vs_ref.paf.gz)
  python3 "$S/auto_names.py" "${A[@]}"
fi
step $T/repeat_chains.fasta python3 "$S/extract_repeat_chains.py" 07_chains2/repeat_chains.txt 06_chains/chains.fasta $T/repeat_chains.fasta
if collapsed; then
  step $D/choices.tsv python3 "$S/assign_arms2.py" --order "$O" --bridges $D/graph.bridges.tsv --ranked $T/bridges/choices.tsv --groups "$GRP" --links 12_broken/links.tsv --agp 12_broken/chains.agp --gfa 00_inputs/unitigs.gfa --paf 11_ont_chains/ont_vs_chains.paf.gz --out $D/choices.tsv --report $D/arm_report.tsv --depth $D/arm_depth.tsv --unused $D/arm_usage.tsv --twoX "$AS_TWO_X" --hap "$AS_HAP" --ploidy "$AS_PLOIDY" --minAllele "$AS_MIN_ALLELE" --maxPair "$AS_MAX_PAIR"
else
  [ -s $D/choices.tsv ] || cp $T/bridges/choices.tsv $D/choices.tsv
fi
BA=(--order "$O" --agp $D/chains.agp --links $D/end_links.tsv --choices $D/choices.tsv --stats $D/chain_stats.tsv --names $T/names.tsv --telomeres $D/telomeres.tsv --unitigs 00_inputs/unitigs.fasta --extraUnplaced $T/repeat_chains.fasta --minWithin "$AS_MIN_WITHIN")
[ -s 12_broken/connector_chains.tsv ] && BA+=(--connectors 12_broken/connector_chains.tsv)
printf '%s\n' "${BA[@]}" > $D/build_args.txt
step $D/pre/$AS_SAMPLE.junctions.tsv python3 "$S/build_scaffolds.py" "${BA[@]}" --out $D/pre/$AS_SAMPLE
# Hi-C gaps a swap of two chains between homologs would close, on the pre-build: a swap is a group edit and `reset 08_bridges`, cheapest now (D16, D13)
RS=(); [ -s 10_groups/groups.rescued.tsv ] && RS=(--rescued 10_groups/groups.rescued.tsv)
# Homologs written in opposite directions, on the pre-build: listed in edits/flip.txt at D16, stage 12 applies the flip before it maps reads (D12)
step $D/homolog_direction.tsv python3 "$S/review/homolog_direction.py" $D/pre/$AS_SAMPLE.agp 03_unitig_align/unitigs_self.paf.gz $D/homolog_direction.tsv
step $D/rescue_swaps.tsv python3 "$S/review/rescue_swaps.py" $D/pre/$AS_SAMPLE.junctions.tsv $D/end_links.tsv $D/graph.bridges.tsv $D/rescue_swaps.tsv "${RS[@]}"
step $D/tails.tsv python3 "$S/make_tails.py" --order "$O" --junctions $D/pre/$AS_SAMPLE.junctions.tsv --placement $D/pre/$AS_SAMPLE.placement.tsv --links $D/end_links.tsv --groups "$GRP" --chroms 10_groups/final.chroms.tsv --choices $D/choices.tsv --telomeres $D/telomeres.tsv --hic 12_broken/links.tsv --class 11_ont_chains/chain_class.tsv --agp 12_broken/chains.agp --gfa 00_inputs/unitigs.gfa --paf 11_ont_chains/ont_vs_chains.paf.gz --exclude 12_broken/connector_chains.tsv --out $D/tails.tsv --report $D/tail_report.tsv --twoX "$AS_TWO_X" --hap "$AS_HAP" --ploidy "$AS_PLOIDY" --maxDepth "$AS_MAX_DEPTH"
# Genes on each tail (D16): what a tail would add, from compleasm on the unitigs (stage 01)
if has_genes && [ -s $D/tails.tsv ] && [ ! -s $D/tail_genes.tsv ]; then python3 "$S/review/tail_genes.py" --genes "$GENES" --agp 12_broken/chains.agp --placement $D/pre/$AS_SAMPLE.placement.tsv --tails $D/tails.tsv --ploidy "$AS_PLOIDY" --out $D/tail_genes.tsv > $D/tail_genes.txt; fi
[ -s $D/pre/target.fasta ] || cat $D/pre/$AS_SAMPLE.fasta $D/pre/$AS_SAMPLE.unplaced.fasta > $D/pre/target.fasta
# D16 tangle check: placed or tail chains running through a 2x fork whose other exit reaches an unplaced 1x chain (a lost haplotype route).
if collapsed && [ ! -s $D/chain_forks.tsv ]; then
  python3 "$S/review/chain_forks.py" 12_broken/chains.agp 00_inputs/unitigs.gfa 04_copynumber/copy_number.tsv 10_groups/final.tsv 10_groups/final.chroms.tsv --hap "$AS_HAP" --windows 11_ont_chains/depth.windows.tsv --placement $D/pre/$AS_SAMPLE.placement.tsv --chainClass 11_ont_chains/chain_class.tsv --oneLo "$AS_ONE_LO" --oneHi "$AS_SHARED_LO" --tails $D/tails.tsv > $D/chain_forks.tsv.tmp && mv $D/chain_forks.tsv.tmp $D/chain_forks.tsv
fi
# settings.raw_phase_early: do raw reads phase the arbitrary arms (D15)? The stage 13 test, run on the pre-build so that D15 can use it:
# raw reads >= 5 kb mapped to each haplotype's pre-build (shared sequence masked), then deleted. Costs a raw mapping (about 2.5 h on 292).
if collapsed && has_raw && [ "${AS_RAW_PHASE_EARLY:-false}" = "true" ] && [ ! -s $D/phase_reads.pre.tsv ] \
   && awk -F'\t' '$5=="bridge" && $11 ~ /^arbitrary/{f=1} END{exit !f}' $D/pre/$AS_SAMPLE.junctions.tsv; then
  python3 "$S/mask_shared.py" $D/pre/$AS_SAMPLE.agp $D/pre/target.fasta $D/pre/rawmask
  PP=(); SF=()
  for h in $(haps_of $D/pre/$AS_SAMPLE.agp); do map_stream $D/pre/raw.$h.paf.gz $D/pre/rawmask.$h.fasta -; PP+=($D/pre/raw.$h.paf.gz); SF+=(_$h); done
  python3 "$S/review/phase_reads.py" $D/pre/$AS_SAMPLE.agp 12_shared/chains.agp $D/pre/$AS_SAMPLE.junctions.tsv 12_shared/groups.tsv $D/phase_reads.pre.tsv.tmp --paf "${PP[@]}" --suffix "${SF[@]}"
  mv $D/phase_reads.pre.tsv.tmp $D/phase_reads.pre.tsv; rm -f $D/pre/raw.h*.paf.gz* $D/pre/rawmask.*
fi
# Genes on chains left unplaced and off every tail, with the lost routes leading to them (D16): what leaving them out costs
if has_genes && [ ! -s $D/unplaced_genes.tsv ]; then
  python3 "$S/review/tail_genes.py" --genes "$GENES" --agp 12_broken/chains.agp --placement $D/pre/$AS_SAMPLE.placement.tsv --ploidy "$AS_PLOIDY" --unplacedOut $D/unplaced_genes.tsv \
    $([ -s $D/tails.tsv ] && echo "--tails $D/tails.tsv") $([ -s $D/chain_forks.tsv ] && echo "--forks $D/chain_forks.tsv") > $D/unplaced_genes.txt
fi
if [ -s $D/tails.tsv ] && [ "$(grep -vc '^#' $D/tails.tsv)" -gt 0 ]; then
  step $D/tail_checks.tsv python3 "$S/rawreads/tail_checks.py" --tails $D/tails.tsv --chains 12_broken/chains.fasta --target $D/pre/target.fasta --names $T/names.tsv --class 11_ont_chains/chain_class.tsv --report $D/tail_report.tsv --hapDepth "$AS_HAP" --ploidy "$AS_PLOIDY" --minimap2 "$AS_MINIMAP2" --threads "$AS_THREADS" --out $D/tail_checks.tsv
fi
