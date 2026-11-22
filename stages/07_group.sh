#!/bin/bash
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Group chains by Hi-C under the allele conflict constraint. Reference check is an alarm only. Stops for D05, D06, D14.
source "$AS_REPO/stages/lib.sh"
ARGS=(--chroms "$AS_CHROMS" --ploidy "$AS_PLOIDY" --seedLen "$AS_SEED_LEN" --minWithin "$AS_MIN_WITHIN" --stats 10_groups/chain_stats.tsv --repeat 10_groups/chain_repeat.tsv)
[ "$AS_HAP_FIRST" = "true" ] && ARGS+=(--hapFirst)
if collapsed; then
  [ -s 10_groups/shared.txt ] || { if has_rows edits/shared.txt; then cp edits/shared.txt 10_groups/shared.txt; else cp 10_groups/shared.candidates.txt 10_groups/shared.txt; fi; }
  ARGS+=(--shared 10_groups/shared.txt --classDepth 11_ont_chains/chain_class.tsv --hap "$AS_HAP" --endLinks 10_groups/end_links.tsv)
fi
if [ ! -s 10_groups/groups.tsv ]; then python3 "$S/hic_group.py" 12_broken/links.tsv 10_groups/alleles.pairs.tsv 12_broken/chains.fasta.fai 10_groups/groups "${ARGS[@]}" | tee 10_groups/hic_group.log; fi
# settings.hic_allele_rescue: chains left unassigned, placed by allele-ambiguous Hi-C; calls in 10_groups/groups.rescued.tsv (D04, D05).
# When it is off, the same step runs on a copy (10_groups/rescue_preview/) so that D05 sees what turning it on would place:
# that costs `reset 07_group`, since stages 05-06 made the pairs either way.
RP=(12_broken/pairs.rescued.tsv.gz 12_broken/pairs.tsv.gz 10_groups/alleles.pairs.tsv 10_groups/chain_repeat.tsv); RT=""
if hic_rescue; then
  [ -s 10_groups/groups.rescued.tsv ] || python3 "$S/rescue/rescue_groups.py" 10_groups/groups "${RP[@]}" | tee -a 10_groups/hic_group.log
  RT=10_groups/groups.rescued.tsv; GB=10_groups/groups.hic.tsv
elif [ -s 12_broken/pairs.rescued.tsv.gz ]; then
  if [ ! -s 10_groups/rescue_preview/groups.rescued.tsv ]; then
    mkdir -p 10_groups/rescue_preview; cp 10_groups/groups.tsv 10_groups/groups.chroms.tsv 10_groups/rescue_preview/
    python3 "$S/rescue/rescue_groups.py" 10_groups/rescue_preview/groups "${RP[@]}" > 10_groups/rescue_preview.txt
  fi
  RT=10_groups/rescue_preview/groups.rescued.tsv; GB=10_groups/groups.tsv
fi
# Genes short of a copy in the groups, and how many the rescue places (D05); needs compleasm on the unitigs (stage 01)
if has_genes && [ ! -s 10_groups/genes_placed.txt ]; then
  python3 "$S/review/genes_placed.py" --genes "$GENES" --agp 12_broken/chains.agp --groups "${GB:-10_groups/groups.tsv}" --ploidy "$AS_PLOIDY" --out 10_groups/genes_placed.tsv ${RT:+--rescued "$RT"} > 10_groups/genes_placed.txt
fi
if has_ref; then
  B="${AS_BASELINE_AGP:-/dev/null}"
  python3 "$S/eval_groups.py" 10_groups/groups.tsv 12_broken/chains.agp 03_unitig_align/unitigs_vs_ref.paf.gz "$B" > 10_groups/eval_groups.txt || true
fi
# Chains running through a collapsed fork cannot be seen until groups exist. A break costs a read mapping, so read this before spending time on order (D06).
if collapsed && [ -s 10_groups/groups.chroms.tsv ]; then python3 "$S/review/chain_forks.py" 12_broken/chains.agp 00_inputs/unitigs.gfa 04_copynumber/copy_number.tsv 10_groups/groups.tsv 10_groups/groups.chroms.tsv --hap "$AS_HAP" > 10_groups/chain_forks.txt; fi
