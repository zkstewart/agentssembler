#!/bin/bash
# Depth per unitig, for decision D01.
source "$AS_REPO/stages/lib.sh"
mkdir -p 04_copynumber
step 04_copynumber/depth.windows.tsv python3 "$S/unitig_depth.py" 02_ont_map/ont_vs_unitigs.paf.gz 00_inputs/unitigs.fasta.fai 04_copynumber/depth
step 04_copynumber/class_depth.tsv python3 "$S/class_depth.py" 04_copynumber/depth.windows.tsv 03_unitig_align/unitigs_self.paf.gz 04_copynumber/class_depth.tsv
