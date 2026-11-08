#!/bin/bash
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Haplotype-aware gap filling from the hifiasm unitig graph, with Hi-C arbitration and checks.
# Runs after run_pipeline.sh (and apply_sv.sh for the SV version). See REPLICATION.md.
#
# Usage: gapfill.sh <config.sh> [sv|main]      (default: sv if OUTDIR/13_sv exists, else main)
# Needs GFA (hifiasm *.p_utg.gfa, the graph of UNITIGS) in the config.
# Output: OUTDIR/15_gapfill/<sv|main>/
#   <name>.gapfilled.{agp,fasta}       the gap-filled assembly (components are still unitigs)
#   <name>.gapfilled.unplaced.fasta    unplaced unitigs minus those used as fills
#   gapfill.report.tsv                 one row per gap: method, path, Hi-C scores, notes
#   gapfill.validate.tsv               Hi-C own-haplotype share of each fill vs its flanks
#   gapfill.unphased.bed               fills without independent Hi-C phase support
#   crosstalk_{before,after}.tsv       homolog crosstalk (must not rise)
#   telomeres, confidence, juicebox/
set -eo pipefail

[[ $# -ge 1 && -f $1 ]] || { echo "Usage: $0 <config.sh> [sv|main]" >&2; exit 1; }
PIPE=$(cd "$(dirname "$0")" && pwd)
source "${PIPE}/lib.sh"
load_config "$(readlink -f "$1")"
cd "${OUTDIR}"
TARGET=${2:-$([[ -d 13_sv ]] && echo sv || echo main)}
case ${TARGET} in
    sv)   AGP=13_sv/${SAMPLE}.hic_phased.sv.agp; NAME=${SAMPLE}.hic_phased.sv ;;
    main) AGP=10_final/${SAMPLE}.hic_phased.agp; NAME=${SAMPLE}.hic_phased ;;
    *) echo "Target must be sv or main" >&2; exit 1 ;;
esac
[[ -f ${AGP} ]] || { echo "ERROR: ${AGP} not found (run run_pipeline.sh / apply_sv.sh first)" >&2; exit 1; }
[[ -f ${GFA} ]] || { echo "ERROR: GFA '${GFA}' not found; set GFA to hifiasm's *.p_utg.gfa in the config" >&2; exit 1; }
OUT=15_gapfill/${TARGET}; mkdir -p ${OUT}

# The GFA must describe the same unitigs as the FASTA
use_env ${ENV_PYTHON}
python3 - "${GFA}" 00_ref/utg.fasta.fai <<'EOF'
import sys
gfa = {}
for l in open(sys.argv[1]):
    if l[0] == "S":
        s = l.split("\t"); ln = [ x for x in s[3:] if x.startswith("LN:i:") ]
        gfa[s[1]] = int(ln[0][5:].strip()) if ln else len(s[2])
fai = { l.split("\t")[0]: int(l.split("\t")[1]) for l in open(sys.argv[2]) }
bad = [ u for u in fai if gfa.get(u) != fai[u] ]
if bad: sys.exit(f"ERROR: {len(bad)} unitigs differ between the GFA and UNITIGS, e.g. {bad[:3]}")
EOF

echo "# Filling gaps of ${AGP}"
python3 ${PIPE}/gapfill.py --agp ${AGP} --gfa "${GFA}" --fasta 00_ref/utg.fasta --links 03_links/links.tsv \
    --footprints 02_refplace/phasing.tsv --minimap2 ${MINIMAP2} --threads ${THREADS} \
    --window ${GAPFILL_WINDOW:-1000000} --minLinks ${GAPFILL_MIN_LINKS:-10} --ratio ${GAPFILL_RATIO:-2.0} \
    --minEquivIdentity ${GAPFILL_EQUIV_IDENTITY:-0.999} --out ${OUT}/gapfill | tee ${OUT}/gapfill.summary.txt
mv ${OUT}/gapfill.agp ${OUT}/${NAME}.gapfilled.agp
mv ${OUT}/gapfill.unphased.bed ${OUT}/gapfill.unphased.bed.tmp && sort -k1,1V -k2,2n ${OUT}/gapfill.unphased.bed.tmp > ${OUT}/gapfill.unphased.bed && rm ${OUT}/gapfill.unphased.bed.tmp
python3 ${PIPE}/agp_to_fasta.py ${OUT}/${NAME}.gapfilled.agp 00_ref/utg.fasta ${OUT}/${NAME}.gapfilled.fasta

# Unitig bookkeeping: every original unitig is still placed, and new unitigs appear once
# (except homozygous copies the report marks 'shared')
lost=$(comm -23 <(awk -F'\t' '$5=="W"{print $6}' ${AGP} | sort -u) <(awk -F'\t' '$5=="W"{print $6}' ${OUT}/${NAME}.gapfilled.agp | sort -u) | wc -l)
extra=$(comm -13 <(awk -F'\t' '$5=="W"{print $6}' ${AGP} | sort -u) <(awk -F'\t' '$5=="W"{print $6}' ${OUT}/${NAME}.gapfilled.agp | sort -u) \
    | awk 'NR==FNR{f[$1]; next} !($1 in f)' ${OUT}/gapfill.membership.tsv - | wc -l)
[[ ${lost} -eq 0 && ${extra} -eq 0 ]] || { echo "ERROR: ${lost} original unitigs lost, ${extra} unexplained new unitigs" >&2; exit 1; }
dups=$(awk -F'\t' '$5=="W"{print $6}' ${OUT}/${NAME}.gapfilled.agp | sort | uniq -d | wc -l)
shared=$(grep -c '_shared$' ${OUT}/gapfill.membership.tsv || true)
echo "# unitigs placed twice: ${dups} (homozygous copies marked shared: ${shared})"
use_env ${ENV_ASSEMBLY}
cut -f1 ${OUT}/gapfill.membership.tsv | grep -v '^#' > ${OUT}/fill_unitigs.txt || true
seqkit grep -v -f ${OUT}/fill_unitigs.txt 10_final/${SAMPLE}.hic_phased.unplaced.fasta > ${OUT}/${NAME}.gapfilled.unplaced.fasta 2> /dev/null
rm -f ${OUT}/fill_unitigs.txt

# Hi-C checks: each fill against its flanks, and homolog crosstalk before vs after
use_env ${ENV_PYTHON}
python3 ${PIPE}/gapfill_validate.py ${OUT}/${NAME}.gapfilled.agp ${OUT}/gapfill.report.tsv 03_links/pairs.tsv.gz > ${OUT}/gapfill.validate.tsv
awk -F'\t' '!/^#/{n[$13]++} END{for(k in n) printf "# fills, Hi-C %s: %d\n", k, n[k]}' ${OUT}/gapfill.validate.tsv
awk -F'\t' '$5=="W"{print $1}' ${AGP} | sort -uV | awk '{k=$1; sub(/_[^_]*$/, "", k); g[k]=g[k] (g[k]?",":"") $1} END{for(k in g) print k"\t"g[k]}' \
    | sort -V > ${OUT}/homolog_groups.tsv
use_env ${ENV_HAPHIC}
for x in before:${AGP} after:${OUT}/${NAME}.gapfilled.agp; do
    mkdir -p ${OUT}/maps_${x%%:*}
    python3 ${PIPE}/contact_map.py --pairs 03_links/pairs.tsv.gz --agp ${x#*:} --groups ${OUT}/homolog_groups.tsv \
        --out ${OUT}/maps_${x%%:*} > /dev/null
    cp ${OUT}/maps_${x%%:*}/crosstalk.tsv ${OUT}/crosstalk_${x%%:*}.tsv
    awk -F'\t' -v w=${x%%:*} '!/^#/{i+=$3; h+=$4} END{printf "# homolog crosstalk %s: %.4f\n", w, h/(i+h)}' ${OUT}/crosstalk_${x%%:*}.tsv
done

# Telomeres, confidence (fills count as weak unless Hi-C chose them), Juicebox
use_env ${ENV_PYTHON}
python3 ${PIPE}/telomeres.py ${OUT}/${NAME}.gapfilled.fasta --motif ${TELOMERE_MOTIF} > ${OUT}/${NAME}.gapfilled.telomeres.tsv
awk 'NR>1{n+=gsub(/T/,"",$5); e+=2} END{print "# telomeric ends:", n, "of", e}' ${OUT}/${NAME}.gapfilled.telomeres.tsv
cat ${OUT}/gapfill.membership.tsv 10_final/membership.tsv > ${OUT}/membership.tsv # original rows win for shared unitigs
python3 ${PIPE}/confidence.py ${OUT}/${NAME}.gapfilled.agp ${OUT}/membership.tsv 05_phase/phase.support.tsv \
    03_links/links.tsv > ${OUT}/${NAME}.gapfilled.confidence.tsv
tail -1 ${OUT}/${NAME}.gapfilled.confidence.tsv
awk -F'\t' '!/^#/ && $5!="W"{n++} END{print "# gaps remaining:", n}' ${OUT}/${NAME}.gapfilled.agp
make_juicebox ${OUT}/${NAME}.gapfilled.agp ${OUT}/juicebox
echo "# Done: ${OUTDIR}/${OUT}"
