#!/bin/bash
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Final scaffolds with the accepted tails, if any.
source "$AS_REPO/stages/lib.sh"
D=12_shared; T=13_scaffold
mapfile -t BA < $D/build_args.txt
if has_rows edits/accepted_tails.tsv; then BA+=(--tails edits/accepted_tails.tsv); fi
python3 "$S/build_scaffolds.py" "${BA[@]}" --out $T/$AS_SAMPLE.graph > $T/build.log
cat $T/$AS_SAMPLE.graph.fasta $T/$AS_SAMPLE.graph.unplaced.fasta > $T/target.fasta
"$AS_SAMTOOLS" faidx $T/target.fasta; "$AS_SAMTOOLS" faidx $T/$AS_SAMPLE.graph.fasta
# Homologs written in opposite directions, again with the tails (stage 10 checked the pre-build): a flip is cheap before stage 12 maps reads, costly after (D12)
python3 "$S/review/homolog_direction.py" $T/$AS_SAMPLE.graph.agp 03_unitig_align/unitigs_self.paf.gz $T/homolog_direction.tsv >> $T/build.log
# Hi-C gaps a swap of two chains between homologs would close: usually the allele rescue calling against the graph (D10, D13)
RS=(); [ -s 10_groups/groups.rescued.tsv ] && RS=(--rescued 10_groups/groups.rescued.tsv)
python3 "$S/review/rescue_swaps.py" $T/$AS_SAMPLE.graph.junctions.tsv $D/end_links.tsv $D/graph.bridges.tsv $T/rescue_swaps.tsv "${RS[@]}" >> $T/build.log
