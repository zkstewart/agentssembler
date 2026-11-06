#!/bin/bash
# Final scaffolds with the accepted tails, if any.
source "$AS_REPO/stages/lib.sh"
D=12_shared; T=13_scaffold
mapfile -t BA < $D/build_args.txt
if has_rows edits/accepted_tails.tsv; then BA+=(--tails edits/accepted_tails.tsv); fi
python3 "$S/build_scaffolds.py" "${BA[@]}" --out $T/$AS_SAMPLE.graph > $T/build.log
cat $T/$AS_SAMPLE.graph.fasta $T/$AS_SAMPLE.graph.unplaced.fasta > $T/target.fasta
"$AS_SAMTOOLS" faidx $T/target.fasta; "$AS_SAMTOOLS" faidx $T/$AS_SAMPLE.graph.fasta
