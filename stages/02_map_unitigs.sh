#!/bin/bash
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Corrected reads to all unitigs, unitig self-alignment, unitigs to the reference (alarm only). Takes hours; can run elsewhere and be put in place.
source "$AS_REPO/stages/lib.sh"
mkdir -p 02_ont_map 03_unitig_align
U=00_inputs/unitigs.fasta
map_paf 02_ont_map/ont_vs_unitigs.paf.gz -x lr:hq -c -t "$AS_THREADS" --secondary=yes -N 5 -p 0.9 "$U" 00_inputs/ont_corrected_reads.fasta
map_paf 03_unitig_align/unitigs_self.paf.gz -x asm20 -c -t "$AS_THREADS" -D --dual=no --secondary=yes -N 20 -p 0.3 "$U" "$U"
if has_ref; then map_paf 03_unitig_align/unitigs_vs_ref.paf.gz -x asm20 -c -t "$AS_THREADS" --secondary=yes -N 5 00_inputs/reference.fasta "$U"; fi
