#!/bin/bash
# Everything measured on the final files. Corrected reads are remapped to final coordinates. Stops for D12, D13: the human Juicebox review and SV calls.
source "$AS_REPO/stages/lib.sh"
V=16_validate; F=15_final; [ -s 15_final_raw/$AS_SAMPLE.final.agp ] && F=15_final_raw
AGP=$F/$AS_SAMPLE.final.agp; FA=$F/$AS_SAMPLE.final.fasta; UN=$F/$AS_SAMPLE.final.unplaced.fasta; D=12_shared
mkdir -p $V/maps/specific $V/maps/scaffold $V/juicebox $V/sv/inv $V/reads $V/compleasm
python3 "$S/agp_stats.py" $AGP | tee $V/stats.txt
python3 "$S/telomeres.py" $FA --motif "$AS_MOTIF" > $V/telomeres.final.tsv
has_ref && python3 "$S/dosage.py" $AGP 03_unitig_align/unitigs_vs_ref.paf.gz 00_inputs/reference.fasta.fai --ploidy "$AS_PLOIDY" > $V/dosage.txt 2>&1 || true
has_ref && python3 "$S/sv_candidates.py" $AGP 03_unitig_align/unitigs_vs_ref.paf.gz $V/sv_candidates.tsv > /dev/null 2>&1 || true
awk 'BEGIN{OFS="\t"}{print $1,1,$2,1,"W",$1,1,$2,"+"}' $FA.fai > $V/identity.agp
python3 "$S/lift_pairs.py" $D/pairs.tsv.gz $D/chains.agp $AGP $D/groups.tsv 13_scaffold/names.tsv $V/pairs.scaffold.tsv.gz $V/pairs.specific.tsv.gz > $V/lift.log 2>&1
cut -f1 $FA.fai | awk '{c=$1; sub(/_h[0-9]+$/,"",c); a[c]=a[c]","$1; if(!(c in s)){s[c]=1; o[++n]=c}} END{for(i=1;i<=n;i++) print o[i]"\t"substr(a[o[i]],2)}' > $V/homolog_groups.tsv
if [ -n "${AS_PLOT_PYTHON:-}" ]; then for t in specific scaffold; do "$AS_PLOT_PYTHON" "$S/contact_map.py" --pairs $V/pairs.$t.tsv.gz --agp $V/identity.agp --groups $V/homolog_groups.tsv --out $V/maps/$t > $V/maps/$t.log 2>&1 || echo "contact_map $t failed: see log"; done; fi
python3 "$S/validate/phase_scan.py" $V/pairs.specific.tsv.gz $AGP > $V/phase_scan.txt 2>&1
# corrected reads to the final coordinates
cat $FA $UN > $V/reads/all.fasta; PAFS=()
if collapsed; then
  python3 "$S/mask_shared.py" $AGP $V/reads/all.fasta $V/reads/target > $V/reads/mask.log 2>&1; HAPS=$(haps_of $AGP)
  for h in $HAPS; do other=$(echo $HAPS | tr ' ' '\n' | grep -vx "$h" | sed 's/^/_/;s/$/$/' | paste -sd'|')
    map_paf $V/reads/all.$h.paf.gz -x lr:hq -t "$AS_THREADS" --secondary=no $V/reads/target.$h.fasta 00_inputs/ont_corrected_reads.fasta
    [ -s $V/reads/ont.$h.paf.gz ] || zcat $V/reads/all.$h.paf.gz | awk -F'\t' -v o="$other" '$6 !~ o' | gzip -1 > $V/reads/ont.$h.paf.gz; PAFS+=($V/reads/ont.$h.paf.gz); done
  rm -f $V/reads/target.h*.fasta*
  python3 "$S/validate/depth_classes.py" $AGP $V/reads/target.masked.tsv "${PAFS[@]}" $V/depth.windows.tsv > $V/depth_classes.log 2>&1
else
  map_paf $V/reads/ont.paf.gz -x lr:hq -t "$AS_THREADS" --secondary=no $V/reads/all.fasta 00_inputs/ont_corrected_reads.fasta; PAFS+=($V/reads/ont.paf.gz)
  "$AS_SAMTOOLS" faidx $V/reads/all.fasta
  python3 "$S/depth_flat.py" $V/reads/ont.paf.gz $V/reads/all.fasta.fai > $V/depth_flat.txt 2>&1
fi
rm -f $V/reads/all.fasta $V/reads/all.fasta.fai
# homologs against each other: inverted blocks are SV candidates, and a wholly inverted scaffold is a direction flag to fix (edits/flip.txt)
python3 - "$FA.fai" "$V" "$FA" "$AS_MINIMAP2" "$AS_THREADS" "$AS_SAMTOOLS" <<'PY'
import sys, collections, subprocess, os
fai, V, fa, mm, t, st = sys.argv[1:7]; H = collections.defaultdict(list)
for l in open(fai):
    n = l.split("\t")[0]; c, _, h = n.rpartition("_h")
    if c: H[c].append(n)
for c, ns in H.items():
    ns.sort(); ref = ns[0]; r = f"{V}/sv/inv/{c}.ref.fa"
    open(r, "w").write(subprocess.run([st, "faidx", fa, ref], capture_output=True, text=True, check=True).stdout)
    for q in ns[1:]:
        qf = f"{V}/sv/inv/{c}.q.fa"; open(qf, "w").write(subprocess.run([st, "faidx", fa, q], capture_output=True, text=True, check=True).stdout)
        with open(f"{V}/sv/inv/{q}_vs_{ref}.paf", "w") as o: subprocess.run([mm, "-x", "asm20", "-t", t, "--secondary=no", r, qf], stdout=o, stderr=subprocess.DEVNULL, check=True)
        os.remove(qf)
    os.remove(r)
PY
python3 "$S/merge_inverted.py" $V/sv/inv/*_vs_*.paf > $V/sv/inv/merged.tsv
if [ -s $V/sv/inv/merged.tsv ]; then
  B=$(cut -f1,2 $V/sv/inv/merged.tsv | tr '\t' ' ' | tr '\n' ' ')
  python3 "$S/inversion_test.py" $V/identity.agp $V/pairs.scaffold.tsv.gz $B > $V/sv/inversion_test.txt 2>&1 || true
  SP=$(IFS=,; echo "${PAFS[*]}")
  python3 "$S/bp_reads.py" $AGP "$SP" $B > $V/sv/bp_reads.tsv 2>&1 || true
fi
if has_ref; then "$AS_MINIMAP2" -x asm20 -K 100M -t "$AS_THREADS" --secondary=no 00_inputs/reference.fasta $FA 2>/dev/null > $V/sv/final_vs_ref.paf || echo "reference alignment failed; only an alarm, so continuing"; fi
if [ -s $D/tails.tsv ] && has_rows edits/accepted_tails.tsv; then python3 "$S/rawreads/tail_joins.py" --tails edits/accepted_tails.tsv --chainsAgp $D/chains.agp --agp $AGP --names 13_scaffold/names.tsv --paf "${PAFS[@]}" --out $V/tail_joins.tsv > $V/tail_joins.log 2>&1 || true; fi
if [ -n "${AS_JUICER_TOOLS:-}" ] && [ ! -e $V/juicebox/done.ok ]; then
 ( cd $V/juicebox; python3 "$S/pairs_to_jbat.py" "$AS_WORK/$AGP" ../pairs.scaffold.tsv.gz out_JBAT && java -Djava.awt.headless=true -Xmx8g -jar "$AS_JUICER_TOOLS" pre out_JBAT.txt out_JBAT.hic.part out_JBAT.sizes > juicer_tools.log 2>&1 && mv out_JBAT.hic.part out_JBAT.hic && rm -f out_JBAT.txt && touch done.ok )
fi
if [ -n "${AS_COMPLEASM:-}" ]; then
 ( cd $V/compleasm
   for h in $(haps_of "$AS_WORK/$AGP"); do [ -e $h/summary.txt ] && continue; "$AS_SEQKIT" grep -r -p "_$h\$" "$AS_WORK/$FA" > $h.fasta 2>/dev/null; "$AS_COMPLEASM" run -a $h.fasta -o $h -l "$AS_LINEAGE" -L "$AS_BUSCO_LIB" -t "$AS_THREADS" > $h.log 2>&1 || echo "compleasm failed for $h"; rm -f $h.fasta; done
   [ -e unplaced/summary.txt ] || "$AS_COMPLEASM" run -a "$AS_WORK/$UN" -o unplaced -l "$AS_LINEAGE" -L "$AS_BUSCO_LIB" -t "$AS_THREADS" > unplaced.log 2>&1 || echo "compleasm failed for unplaced" )
fi
