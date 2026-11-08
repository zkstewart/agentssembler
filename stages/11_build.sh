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
