#!/bin/bash
# Names, arms at bridges, a pre-build, one-sided tails and their tests. Stops for D15, D16: nothing is placed without edits/accepted_tails.tsv.
source "$AS_REPO/stages/lib.sh"
D=12_shared; T=13_scaffold; mkdir -p $T $D/pre
O=$D/order.tsv; GRP=$D/groups.tsv; [ -s $D/groups.edited.tsv ] && GRP=$D/groups.edited.tsv
if has_rows edits/names.tsv; then cp edits/names.tsv $T/names.tsv; elif [ ! -s $T/names.tsv ]; then
  A=(--groups "$GRP" --chroms 10_groups/final.chroms.tsv --out $T/names.tsv --order "$O")
  has_ref && A+=(--agp 12_broken/chains.agp --paf 03_unitig_align/unitigs_vs_ref.paf.gz)
  python3 "$S/auto_names.py" "${A[@]}"
fi
step $T/repeat_chains.fasta python3 "$S/extract_repeat_chains.py" 07_chains2/repeat_chains.txt 06_chains/chains.fasta $T/repeat_chains.fasta
if collapsed; then
  step $D/choices.tsv python3 "$S/assign_arms2.py" --order "$O" --bridges $D/graph.bridges.tsv --ranked $T/bridges/choices.tsv --groups "$GRP" --links 12_broken/links.tsv --agp 12_broken/chains.agp --gfa 00_inputs/unitigs.gfa --paf 11_ont_chains/ont_vs_chains.paf.gz --out $D/choices.tsv --report $D/arm_report.tsv --depth $D/arm_depth.tsv --unused $D/arm_usage.tsv --twoX "$AS_TWO_X" --minAllele "$AS_MIN_ALLELE" --maxPair "$AS_MAX_PAIR"
else
  [ -s $D/choices.tsv ] || cp $T/bridges/choices.tsv $D/choices.tsv
fi
BA=(--order "$O" --agp $D/chains.agp --links $D/end_links.tsv --choices $D/choices.tsv --stats $D/chain_stats.tsv --names $T/names.tsv --telomeres $D/telomeres.tsv --unitigs 00_inputs/unitigs.fasta --extraUnplaced $T/repeat_chains.fasta --minWithin "$AS_MIN_WITHIN")
[ -s 12_broken/connector_chains.tsv ] && BA+=(--connectors 12_broken/connector_chains.tsv)
printf '%s\n' "${BA[@]}" > $D/build_args.txt
step $D/pre/$AS_SAMPLE.junctions.tsv python3 "$S/build_scaffolds.py" "${BA[@]}" --out $D/pre/$AS_SAMPLE
step $D/tails.tsv python3 "$S/make_tails.py" --order "$O" --junctions $D/pre/$AS_SAMPLE.junctions.tsv --placement $D/pre/$AS_SAMPLE.placement.tsv --links $D/end_links.tsv --groups "$GRP" --chroms 10_groups/final.chroms.tsv --choices $D/choices.tsv --telomeres $D/telomeres.tsv --hic 12_broken/links.tsv --class 11_ont_chains/chain_class.tsv --agp 12_broken/chains.agp --gfa 00_inputs/unitigs.gfa --paf 11_ont_chains/ont_vs_chains.paf.gz --exclude 12_broken/connector_chains.tsv --out $D/tails.tsv --report $D/tail_report.tsv --twoX "$AS_TWO_X" --maxDepth "$AS_MAX_DEPTH"
[ -s $D/pre/target.fasta ] || cat $D/pre/$AS_SAMPLE.fasta $D/pre/$AS_SAMPLE.unplaced.fasta > $D/pre/target.fasta
if [ -s $D/tails.tsv ] && [ "$(grep -vc '^#' $D/tails.tsv)" -gt 0 ]; then
  step $D/tail_checks.tsv python3 "$S/rawreads/tail_checks.py" --tails $D/tails.tsv --chains 12_broken/chains.fasta --target $D/pre/target.fasta --names $T/names.tsv --class 11_ont_chains/chain_class.tsv --report $D/tail_report.tsv --hapDepth "$AS_HAP" --ploidy "$AS_PLOIDY" --minimap2 "$AS_MINIMAP2" --threads "$AS_THREADS" --out $D/tail_checks.tsv
fi
