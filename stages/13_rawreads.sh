#!/bin/bash
# Optional, with raw basecalled BAMs. Streams reads >= 5 kb from the BAMs to minimap2: never write a raw FASTA. Compares raw with corrected reads at gaps and control sites, re-fills from raw reads near gaps, audits every read-derived join. Stops for D18.
source "$AS_REPO/stages/lib.sh"
has_raw || { log "no raw BAMs: nothing to do"; exit 0; }
A=17_raw_reads; T=13_scaffold; R=14_readfill; F=15_final; mkdir -p $A/fill; AGP=$T/$AS_SAMPLE.graph.agp; MM="$AS_MINIMAP2"
IFS=: read -r -a BAMS <<< "$AS_RAW_BAMS"
stream() { # stream TAB|- : FASTA of reads >= 5 kb to stdout
  local tab="$1"; for b in "${BAMS[@]}"; do "$AS_SAMTOOLS" view -@ 2 "$b"; done | awk -F'\t' -v tab="$tab" 'BEGIN{OFS="\t"; if(tab!="-") cmd="gzip -1 > " tab}
    {l=length($10); if(tab!="-"){q="NA";for(i=12;i<=NF;i++)if($i~/^qs:f:/){q=substr($i,6);break}; print $1,l,q | cmd} if(l>=5000) print ">"$1"\n"$10}'
}
map_stream() { # OUT target tab
  local o="$1" t="$2" tab="$3"; [ -s "$o" ] && { log "skip, exists: $o"; return 0; }
  guard; log "stream map: $o"
  stream "$tab" | "$MM" -x lr:hq -t "$AS_THREADS" --secondary=no "$t" - 2> "$o.log" | gzip -1 > "$o.part" &
  local p=$!
  while kill -0 $p 2>/dev/null; do f=$(df -B1G "$AS_GUARD_PATH" | awk 'NR==2{print $4}'); if [ "$f" -lt "$AS_MIN_FREE_GB" ]; then kill $p; pkill -P $$ || true; echo "ERROR: disk guard" >&2; exit 3; fi; sleep 30; done
  wait $p; gzip -t "$o.part"; mv "$o.part" "$o"
}
[ -s $A/corrected.lens.tsv.gz ] || "$AS_SEQKIT" fx2tab -n -i -l 00_inputs/ont_corrected_reads.fasta | cut -f1,2 | gzip -1 > $A/corrected.lens.tsv.gz
SETS=(); tab=$A/reads.tsv.gz; [ -s $tab ] && tab=-
if collapsed; then
  for h in $(haps_of $AGP); do
    map_stream $A/raw_vs_target.$h.paf.gz $R/target.$h.fasta "$tab"; tab=-
    "$AS_SAMTOOLS" faidx $R/$h/in.fasta
    SETS+=("$R/$h/in.agp,$R/$h/in.fasta.fai,$A/raw_vs_target.$h.paf.gz,$R/$h/ont.paf.gz")
  done
  for h in $(haps_of $AGP); do [ -s $A/raw.$h.paf.gz ] || { other=$(haps_of $AGP | grep -vx "$h" | sed 's/^/_/;s/$/$/' | paste -sd'|'); zcat $A/raw_vs_target.$h.paf.gz | awk -F'\t' -v o="$other" '$6 !~ o' | gzip -1 > $A/raw.$h.paf.gz; }; done
else
  map_stream $A/raw_vs_target.paf.gz $T/target.fasta "$tab"
  SETS+=("$AGP,$T/$AS_SAMPLE.graph.fasta.fai,$A/raw_vs_target.paf.gz,$R/ont_vs_target.paf.gz")
fi
[ -s $A/cmp.by_length.tsv ] || python3 "$S/rawreads/compare_reads.py" $A/cmp $A/reads.tsv.gz $A/corrected.lens.tsv.gz "${SETS[@]}" > $A/compare.log 2>&1
# reads near gaps only, in a small FASTA
if [ ! -s $A/neargap.fasta ]; then
  python3 - "$AGP" $A/neargap.ids $A/raw_vs_target*.paf.gz <<'PY'
import sys, gzip, collections
agp, out, *pafs = sys.argv[1:]; G = collections.defaultdict(list)
for l in open(agp):
    if l[0] == "#": continue
    s = l.split("\t")
    if s[4] in "NU": G[s[0]].append((int(s[1]) - 1, int(s[2])))
ids = set()
for p in pafs:
    for l in gzip.open(p, "rt"):
        s = l.split("\t", 10)
        if s[5] in G and any(int(s[7]) < b + 60000 and int(s[8]) > a - 60000 for a, b in G[s[5]]): ids.add(s[0])
open(out, "w").write("\n".join(sorted(ids)) + "\n")
PY
  for b in "${BAMS[@]}"; do "$AS_SAMTOOLS" view -@ 4 -N $A/neargap.ids "$b" | awk -F'\t' '{print ">"$1"\n"$10}'; done > $A/neargap.part; mv $A/neargap.part $A/neargap.fasta
fi
# re-fill the current final AGP from raw reads, then map to the candidate with CIGARs for the audit
cat $F/$AS_SAMPLE.final.fasta $F/$AS_SAMPLE.final.unplaced.fasta > $A/fill/in.fasta; "$AS_SAMTOOLS" faidx $A/fill/in.fasta
RF=(--seqkit "$AS_SEQKIT" --minimap2 "$MM"); AUD=(); REPORTS=()
if collapsed; then
  python3 "$S/mask_shared.py" $F/$AS_SAMPLE.final.agp $A/fill/in.fasta $A/fill/target
  HAPS=$(haps_of $AGP)
  for h in $HAPS; do
    "$AS_SAMTOOLS" faidx $A/fill/target.$h.fasta
    map_paf $A/fill/raw.$h.paf.gz -x lr:hq -t "$AS_THREADS" --secondary=no $A/fill/target.$h.fasta $A/neargap.fasta
    other=$(echo $HAPS | tr ' ' '\n' | grep -vx "$h" | sed 's/^/_/;s/$/$/' | paste -sd'|')
    d=$A/fill/$h; mkdir -p $d; awk -F'\t' -v h=_$h '$1~/^#/||$1 ~ h"$"' $F/$AS_SAMPLE.final.agp > $d/in.agp
    grep '^>' $A/fill/in.fasta | sed 's/^>//' | awk '{print $1}' | grep -Ev "$other" > $d/names.txt
    "$AS_SEQKIT" grep -f $d/names.txt $A/fill/in.fasta -w 0 > $d/in.fasta 2>/dev/null; "$AS_SAMTOOLS" faidx $d/in.fasta
    zcat $A/fill/raw.$h.paf.gz | awk -F'\t' -v o="$other" '$6 !~ o' | gzip -1 > $d/raw.paf.gz
    python3 "$S/readfill.py" --agp $d/in.agp --fasta $d/in.fasta --paf $d/raw.paf.gz --reads $A/neargap.fasta --out $d/fill "${RF[@]}" --threads "$AS_THREADS" > $d/fill.log 2>&1
    REPORTS+=(--report raw=$d/fill.report.tsv,$d/in.fasta --report corrected=$R/$h/fill.$h.report.tsv,$R/$h/in.fasta)
  done
  for h in $HAPS; do grep -v '^#' $A/fill/$h/fill.agp; done > $A/fill/candidate.agp; cat $A/fill/h*/fill.fills.fasta > $A/fill/fills.fasta
else
  map_paf $A/fill/raw.paf.gz -x lr:hq -t "$AS_THREADS" --secondary=no $A/fill/in.fasta $A/neargap.fasta
  python3 "$S/readfill.py" --agp $F/$AS_SAMPLE.final.agp --fasta $A/fill/in.fasta --paf $A/fill/raw.paf.gz --reads $A/neargap.fasta --out $A/fill/fill --target $A/fill/in.fasta "${RF[@]}" --threads "$AS_THREADS" > $A/fill/fill.log 2>&1
  cp $A/fill/fill.agp $A/fill/candidate.agp; cp $A/fill/fill.fills.fasta $A/fill/fills.fasta
  REPORTS+=(--report raw=$A/fill/fill.report.tsv,$A/fill/in.fasta --report corrected=$R/$AS_SAMPLE.graph.report.tsv,$T/$AS_SAMPLE.graph.fasta)
fi
agp2fa $A/fill/candidate.agp 00_inputs/unitigs.fasta,$R/fills.fasta,$A/fill/fills.fasta $A/fill/candidate.fasta
cat $A/fill/candidate.fasta $F/$AS_SAMPLE.final.unplaced.fasta > $A/fill/audit.fasta
if collapsed; then
  python3 "$S/mask_shared.py" $A/fill/candidate.agp $A/fill/audit.fasta $A/fill/audit
  for h in $(haps_of $AGP); do map_paf $A/fill/audit.$h.paf.gz -x lr:hq -c --eqx -t "$AS_THREADS" --secondary=no $A/fill/audit.$h.fasta $A/neargap.fasta; AUD+=($A/fill/audit.$h.paf.gz); done
else
  map_paf $A/fill/audit.paf.gz -x lr:hq -c --eqx -t "$AS_THREADS" --secondary=no $A/fill/audit.fasta $A/neargap.fasta; AUD+=($A/fill/audit.paf.gz)
fi
python3 "$S/rawreads/audit_joins.py" --agp $A/fill/candidate.agp --paf "${AUD[@]}" "${REPORTS[@]}" --minimap2 "$MM" --hapDepth "$AS_HAP" --out $A/audit.tsv
rm -f $A/fill/target.h*.fasta* $A/fill/audit.h*.fasta* $A/fill/audit.fasta
