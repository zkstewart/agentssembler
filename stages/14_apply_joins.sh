#!/bin/bash
# Apply D18: the raw candidate, less every join in edits/join_decisions.tsv. Without raw reads, the 15_final files stand.
source "$AS_REPO/stages/lib.sh"
has_raw || { log "no raw BAMs: final files stand"; exit 0; }
A=17_raw_reads; F=15_final; O=15_final_raw; mkdir -p $O
if has_rows edits/join_decisions.tsv; then
  # a join may be reversed against the pre-fill AGP when closed in an earlier step, or the current final AGP
  { grep -v '^#' 13_scaffold/$AS_SAMPLE.graph.agp; } > $O/pre.agp
  python3 "$S/apply_join_decisions.py" $A/fill/candidate.agp $O/pre.agp edits/join_decisions.tsv $O/$AS_SAMPLE.final.agp; rm -f $O/pre.agp
else cp $A/fill/candidate.agp $O/$AS_SAMPLE.final.agp; fi
agp2fa $O/$AS_SAMPLE.final.agp 00_inputs/unitigs.fasta,14_readfill/fills.fasta,$A/fill/fills.fasta $O/$AS_SAMPLE.final.fasta
cp $F/$AS_SAMPLE.final.unplaced.fasta $F/$AS_SAMPLE.final.unplaced.fasta.fai $O/
python3 "$S/agp_stats.py" $O/$AS_SAMPLE.final.agp | tee $O/stats.txt
