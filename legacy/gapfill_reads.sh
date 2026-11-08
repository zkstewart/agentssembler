#!/bin/bash
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Haplotype-aware gap filling from long reads, after the graph-based gapfill.sh. See REPLICATION.md.
# The reads are mapped once to the whole gap-filled assembly (all haplotypes + unplaced), so each
# read chooses its haplotype; readfill.py fills a gap only from reads anchored uniquely (MAPQ)
# in both flanks that agree with each other, and excises segments that the reads skip and that
# are duplicated elsewhere (false duplications).
#
# Usage: gapfill_reads.sh <config.sh> [sv|main]      (default: sv if OUTDIR/15_gapfill/sv exists)
# Needs READS (long reads, FASTA/FASTQ) in the config; READS_PRESET is minimap2's -x for them.
# Output: OUTDIR/16_readfill/<sv|main>/
#   <name>.readfilled.{agp,fasta}       the assembly (components: unitigs + rf_* read fills)
#   <name>.readfilled.components.fasta  the rf_* fill sequences (with the unitigs, the AGP's components)
#   <name>.readfilled.unplaced.fasta    unplaced unitigs (unchanged by this step)
#   <name>.readfilled.excised.fasta     unitigs excised as false duplications and no longer placed
#   readfill.report.tsv                 one row per gap: method, read counts, trims, fill length, notes
#   readfill.excised.tsv, .skipped.tsv  segments reads skip: excised (duplicated elsewhere) or only reported
#   readfill.links.tsv                  where uniquely anchored reads go instead at unfilled gap edges
#   reads.paf.gz                        the competitive read mapping (reused if present)
#   crosstalk_after.tsv, telomeres, confidence, juicebox/
set -eo pipefail

[[ $# -ge 1 && -f $1 ]] || { echo "Usage: $0 <config.sh> [sv|main]" >&2; exit 1; }
PIPE=$(cd "$(dirname "$0")" && pwd)
source "${PIPE}/lib.sh"
load_config "$(readlink -f "$1")"
cd "${OUTDIR}"
TARGET=${2:-$([[ -d 15_gapfill/sv ]] && echo sv || echo main)}
case ${TARGET} in
    sv)   NAME=${SAMPLE}.hic_phased.sv ;;
    main) NAME=${SAMPLE}.hic_phased ;;
    *) echo "Target must be sv or main" >&2; exit 1 ;;
esac
IN=15_gapfill/${TARGET}
[[ -f ${IN}/${NAME}.gapfilled.agp ]] || { echo "ERROR: ${IN}/${NAME}.gapfilled.agp not found (run gapfill.sh first)" >&2; exit 1; }
[[ -f ${READS} ]] || { echo "ERROR: READS '${READS}' not found; set READS to the long reads in the config" >&2; exit 1; }
OUT=16_readfill/${TARGET}; mkdir -p ${OUT}

# Competitive mapping: every haplotype plus the unplaced unitigs, no secondary alignments
cat ${IN}/${NAME}.gapfilled.fasta ${IN}/${NAME}.gapfilled.unplaced.fasta > ${OUT}/target.fasta
if [[ ! -s ${OUT}/reads.paf.gz ]]; then
    echo "# Mapping ${READS}"
    ${MINIMAP2} -x ${READS_PRESET:-lr:hq} -t ${THREADS} --secondary=no ${OUT}/target.fasta "${READS}" 2> ${OUT}/reads.minimap2.log \
        | gzip -1 > ${OUT}/reads.paf.gz.part
    mv ${OUT}/reads.paf.gz.part ${OUT}/reads.paf.gz
fi

echo "# Filling gaps of ${IN}/${NAME}.gapfilled.agp from reads"
use_env ${ENV_ASSEMBLY}; SEQKIT=$(command -v seqkit)
use_env ${ENV_PYTHON}
[[ -f ${IN}/${NAME}.gapfilled.fasta.fai ]] || ${SEQKIT} faidx ${IN}/${NAME}.gapfilled.fasta
python3 ${PIPE}/readfill.py --agp ${IN}/${NAME}.gapfilled.agp --fasta ${IN}/${NAME}.gapfilled.fasta \
    --paf ${OUT}/reads.paf.gz --reads "${READS}" --target ${OUT}/target.fasta --gapfillReport ${IN}/gapfill.report.tsv \
    --minimap2 ${MINIMAP2} --seqkit ${SEQKIT} --preset ${READS_PRESET:-lr:hq} --threads ${THREADS} \
    --minMapq ${READFILL_MIN_MAPQ:-20} --minAnchor ${READFILL_MIN_ANCHOR:-5000} --minReads ${READFILL_MIN_READS:-2} \
    --minIdentity ${READFILL_MIN_IDENTITY:-0.98} --maxTrim ${READFILL_MAX_TRIM:-2000} \
    --maxExcise ${READFILL_MAX_EXCISE:-2000000} --minDupCover ${READFILL_MIN_DUP_COVER:-0.9} \
    --out ${OUT}/readfill 2>&1 | grep -v '^\[INFO\]' | tee ${OUT}/readfill.summary.txt
mv ${OUT}/readfill.agp ${OUT}/${NAME}.readfilled.agp
mv ${OUT}/readfill.fills.fasta ${OUT}/${NAME}.readfilled.components.fasta
${SEQKIT} faidx ${OUT}/${NAME}.readfilled.components.fasta
python3 ${PIPE}/agp_to_fasta.py ${OUT}/${NAME}.readfilled.agp 00_ref/utg.fasta,${OUT}/${NAME}.readfilled.components.fasta \
    ${OUT}/${NAME}.readfilled.fasta
cp ${IN}/${NAME}.gapfilled.unplaced.fasta ${OUT}/${NAME}.readfilled.unplaced.fasta
rm -f ${OUT}/target.fasta

# Bookkeeping: every unitig of the input is still placed (some flanks trimmed) unless excised as a
# false duplication; nothing new but read fills. Excised unitigs that are no longer placed anywhere
# are near-identical copies of placed sequence: kept aside, like hifiasm's redundant unitigs.
lost=$(comm -23 <(awk -F'\t' '$5=="W"{print $6}' ${IN}/${NAME}.gapfilled.agp | sort -u) <(awk -F'\t' '$5=="W"{print $6}' ${OUT}/${NAME}.readfilled.agp | sort -u))
bad=$(comm -23 <(echo "${lost}" | grep . | sort -u) <(grep -v '^#' ${OUT}/readfill.excised.tsv | cut -f5 | sort -u) | wc -l)
extra=$(comm -13 <(awk -F'\t' '$5=="W"{print $6}' ${IN}/${NAME}.gapfilled.agp | sort -u) <(awk -F'\t' '$5=="W"{print $6}' ${OUT}/${NAME}.readfilled.agp | sort -u) | grep -vc '^rf_' || true)
[[ ${bad} -eq 0 && ${extra} -eq 0 ]] || { echo "ERROR: ${bad} unitigs lost without excision, ${extra} unexplained new components" >&2; exit 1; }
echo "${lost}" | grep . > ${OUT}/excised_unitigs.txt || true
use_env ${ENV_ASSEMBLY}
${SEQKIT} grep -f ${OUT}/excised_unitigs.txt 00_ref/utg.fasta > ${OUT}/${NAME}.readfilled.excised.fasta 2> /dev/null || true
echo "# excised segments: $(grep -vc '^#' ${OUT}/readfill.excised.tsv || true) components; unitigs no longer placed: $(grep -c . ${OUT}/excised_unitigs.txt || true)"
rm -f ${OUT}/excised_unitigs.txt

# Hi-C homolog crosstalk (the read fills carry no Hi-C; this checks the trims and joins), telomeres,
# confidence (read fills count as no-signal), Juicebox
awk -F'\t' '$5=="W"{print $1}' ${OUT}/${NAME}.readfilled.agp | sort -uV | awk '{k=$1; sub(/_[^_]*$/, "", k); g[k]=g[k] (g[k]?",":"") $1} END{for(k in g) print k"\t"g[k]}' \
    | sort -V > ${OUT}/homolog_groups.tsv
use_env ${ENV_HAPHIC}
mkdir -p ${OUT}/maps_after
python3 ${PIPE}/contact_map.py --pairs 03_links/pairs.tsv.gz --agp ${OUT}/${NAME}.readfilled.agp --groups ${OUT}/homolog_groups.tsv \
    --out ${OUT}/maps_after > /dev/null
cp ${OUT}/maps_after/crosstalk.tsv ${OUT}/crosstalk_after.tsv
for x in before:${IN}/crosstalk_after.tsv after:${OUT}/crosstalk_after.tsv; do
    awk -F'\t' -v w=${x%%:*} '!/^#/{i+=$3; h+=$4} END{printf "# homolog crosstalk %s: %.4f\n", w, h/(i+h)}' ${x#*:}
done
use_env ${ENV_PYTHON}
python3 ${PIPE}/telomeres.py ${OUT}/${NAME}.readfilled.fasta --motif ${TELOMERE_MOTIF} > ${OUT}/${NAME}.readfilled.telomeres.tsv
awk 'NR>1{n+=gsub(/T/,"",$5); e+=2} END{print "# telomeric ends:", n, "of", e}' ${OUT}/${NAME}.readfilled.telomeres.tsv
cat ${OUT}/readfill.membership.tsv ${IN}/membership.tsv > ${OUT}/membership.tsv
python3 ${PIPE}/confidence.py ${OUT}/${NAME}.readfilled.agp ${OUT}/membership.tsv 05_phase/phase.support.tsv \
    03_links/links.tsv > ${OUT}/${NAME}.readfilled.confidence.tsv
tail -1 ${OUT}/${NAME}.readfilled.confidence.tsv
awk -F'\t' '!/^#/ && $5!="W"{n++; h[gensub(/.*_/, "", 1, $1)]++} END{printf "# gaps remaining: %d (", n; for(k in h) printf " %s %d", k, h[k]; print " )"}' ${OUT}/${NAME}.readfilled.agp
cat 00_ref/utg.fasta.fai ${OUT}/${NAME}.readfilled.components.fasta.fai > ${OUT}/components.fai
make_juicebox ${OUT}/${NAME}.readfilled.agp ${OUT}/juicebox ${OUT}/components.fai
echo "# Done: ${OUTDIR}/${OUT}"
