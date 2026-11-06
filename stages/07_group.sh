#!/bin/bash
# Group chains by Hi-C under the allele conflict constraint. Reference check is an alarm only. Stops for D05, D06, D14.
source "$AS_REPO/stages/lib.sh"
ARGS=(--chroms "$AS_CHROMS" --ploidy "$AS_PLOIDY" --stats 10_groups/chain_stats.tsv --repeat 10_groups/chain_repeat.tsv)
[ "$AS_HAP_FIRST" = "true" ] && ARGS+=(--hapFirst)
if collapsed; then
  [ -s 10_groups/shared.txt ] || { if has_rows edits/shared.txt; then cp edits/shared.txt 10_groups/shared.txt; else cp 10_groups/shared.candidates.txt 10_groups/shared.txt; fi; }
  ARGS+=(--shared 10_groups/shared.txt)
fi
if [ ! -s 10_groups/groups.tsv ]; then python3 "$S/hic_group.py" 12_broken/links.tsv 10_groups/alleles.pairs.tsv 12_broken/chains.fasta.fai 10_groups/groups "${ARGS[@]}" | tee 10_groups/hic_group.log; fi
if has_ref; then
  B="${AS_BASELINE_AGP:-/dev/null}"
  python3 "$S/eval_groups.py" 10_groups/groups.tsv 12_broken/chains.agp 03_unitig_align/unitigs_vs_ref.paf.gz "$B" > 10_groups/eval_groups.txt || true
fi
