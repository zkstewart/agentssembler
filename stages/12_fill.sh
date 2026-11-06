#!/bin/bash
# Reads across gaps, then the first final files. Collapsed samples map per haplotype to masked targets and leave --target off, since excision would treat true shared copies as false duplication. Stops for D11.
source "$AS_REPO/stages/lib.sh"
T=13_scaffold; R=14_readfill; F=15_final; mkdir -p $R $F; AGP=$T/$AS_SAMPLE.graph.agp; READS=00_inputs/ont_corrected_reads.fasta
RF=(--seqkit "$AS_SEQKIT" --minimap2 "$AS_MINIMAP2")
if collapsed; then
  [ -s $R/target.masked.tsv ] || python3 "$S/mask_shared.py" $AGP $T/target.fasta $R/target
  HAPS=$(haps_of $AGP); agps=(); fills=()
  for h in $HAPS; do
    "$AS_SAMTOOLS" faidx $R/target.$h.fasta
    map_paf $R/ont_vs_target.$h.paf.gz -x lr:hq -t "$AS_THREADS" --secondary=no $R/target.$h.fasta $READS
  done
  n=$(echo $HAPS | wc -w); t=$(( AS_THREADS / n )); [ $t -ge 1 ] || t=1
  for h in $HAPS; do
   (
    d=$R/$h; mkdir -p $d
    awk -F'\t' -v h=_$h '$1~/^#/||$1 ~ h"$"' $AGP > $d/in.agp
    other=$(echo $HAPS | tr ' ' '\n' | grep -vx "$h" | sed 's/^/_/;s/$/$/' | paste -sd'|')
    grep '^>' $T/target.fasta | sed 's/^>//' | awk '{print $1}' | grep -Ev "$other" > $d/names.txt
    "$AS_SEQKIT" grep -f $d/names.txt $T/target.fasta -w 0 > $d/in.fasta 2>/dev/null; "$AS_SAMTOOLS" faidx $d/in.fasta
    zcat $R/ont_vs_target.$h.paf.gz | awk -F'\t' -v o="$other" '$6 !~ o' | gzip -1 > $d/ont.paf.gz
    [ -s $d/fill.$h.agp ] || python3 "$S/readfill.py" --agp $d/in.agp --fasta $d/in.fasta --paf $d/ont.paf.gz --reads $READS --out $d/fill.$h "${RF[@]}" --threads $t > $d/fill.log 2>&1
   ) &
  done
  wait
  for h in $HAPS; do grep -v '^#' $R/$h/fill.$h.agp; done > $R/merged.agp
  cat $R/h*/fill.h*.fills.fasta > $R/fills.fasta
  # D11: re-fill listed scaffolds at another MAPQ (a shared flank can leave anchored reads at MAPQ 0)
  if has_rows edits/fill_retry.tsv; then
    grep -v '^#' edits/fill_retry.tsv | while IFS=$'\t' read -r sc mq why; do
      h=${sc##*_}; d=$R/retry_$sc; mkdir -p $d; awk -F'\t' -v s="$sc" '$1~/^#/||$1==s' $R/$h/in.agp > $d/in.agp
      python3 "$S/readfill.py" --agp $d/in.agp --fasta $R/$h/in.fasta --paf $R/$h/ont.paf.gz --reads $READS --out $d/fill "${RF[@]}" --threads "$AS_THREADS" --minMapq "$mq" > $d/fill.log 2>&1
      awk -F'\t' -v s="$sc" '$1!=s' $R/merged.agp > $R/merged.tmp; grep -v '^#' $d/fill.agp >> $R/merged.tmp; mv $R/merged.tmp $R/merged.agp
      cat $d/fill.fills.fasta >> $R/fills.fasta
    done
  fi
  sort -t$'\t' -k1,1 -k2,2n $R/merged.agp > $F/$AS_SAMPLE.final.agp
else
  map_paf $R/ont_vs_target.paf.gz -x lr:hq -t "$AS_THREADS" --secondary=no $T/target.fasta $READS
  [ -s $R/$AS_SAMPLE.graph.agp ] || python3 "$S/readfill.py" --agp $AGP --fasta $T/$AS_SAMPLE.graph.fasta --paf $R/ont_vs_target.paf.gz --reads $READS --out $R/$AS_SAMPLE.graph "${RF[@]}" --threads "$AS_THREADS" --target $T/target.fasta > $R/fill.log 2>&1
  cp $R/$AS_SAMPLE.graph.agp $F/$AS_SAMPLE.final.agp; cp $R/$AS_SAMPLE.graph.fills.fasta $R/fills.fasta
fi
agp2fa $F/$AS_SAMPLE.final.agp 00_inputs/unitigs.fasta,$R/fills.fasta $F/$AS_SAMPLE.final.fasta
# D-names trap: a stale direction flag leaves a scaffold reversed against its homolog; list it in edits/flip.txt
if has_rows edits/flip.txt; then python3 "$S/flip_scaffold.py" $F/$AS_SAMPLE.final.agp $F/$AS_SAMPLE.final.fasta $F/flip.agp $F/flip.fasta $(grep -v '^#' edits/flip.txt | cut -f1); mv $F/flip.agp $F/$AS_SAMPLE.final.agp; mv $F/flip.fasta $F/$AS_SAMPLE.final.fasta; "$AS_SAMTOOLS" faidx $F/$AS_SAMPLE.final.fasta; fi
cp $T/$AS_SAMPLE.graph.unplaced.fasta $F/$AS_SAMPLE.final.unplaced.fasta; "$AS_SAMTOOLS" faidx $F/$AS_SAMPLE.final.unplaced.fasta
python3 "$S/agp_stats.py" $F/$AS_SAMPLE.final.agp | tee $F/stats.txt
