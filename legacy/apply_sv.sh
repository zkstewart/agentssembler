#!/bin/bash
# Apply structural-variant models to the main assembly from run_pipeline.sh.
#
# Usage: apply_sv.sh <config.sh>
# SV_DIR (set in the config) holds the decisions, made by inspecting contact maps and
# sv_breakpoints.py output (see REPLICATION.md):
#   *.spec     rearrange.py specs (e.g. reciprocal translocations), applied in file-name order
#   moves.tsv  optional insert_sv.py moves: unitig<TAB>target pseudochromosome
# Output: OUTDIR/13_sv/<SAMPLE>.hic_phased.sv.{agp,fasta}, telomeres, confidence, juicebox/
set -eo pipefail

[[ $# -ge 1 && -f $1 ]] || { echo "Usage: $0 <config.sh>" >&2; exit 1; }
PIPE=$(cd "$(dirname "$0")" && pwd)
source "${PIPE}/lib.sh"
load_config "$(readlink -f "$1")"
cd "${OUTDIR}"
OUT=13_sv; MAIN=10_final/${SAMPLE}.hic_phased.agp
mkdir -p ${OUT}; use_env ${ENV_PYTHON}

agp=${MAIN}; k=0
for spec in $(ls "${SV_DIR}"/*.spec | sort); do
    k=$((k + 1))
    echo "# Applying $(basename ${spec})"
    python3 ${PIPE}/rearrange.py --agp ${agp} --pairs 03_links/pairs.tsv.gz --spec ${spec} --out ${OUT}/step${k}
    agp=${OUT}/step${k}.agp
done
MOVES=""
if [[ -f ${SV_DIR}/moves.tsv ]]; then
    echo "# Moving single unitigs"
    python3 ${PIPE}/insert_sv.py --agp ${agp} --pairs 03_links/pairs.tsv.gz --moves ${SV_DIR}/moves.tsv --out ${OUT}/moved
    agp=${OUT}/moved.agp; MOVES="--moves ${SV_DIR}/moves.tsv"
fi
python3 ${PIPE}/finalize_sv.py --agp ${agp} --main ${MAIN} ${MOVES} --telomeric 10_final/telomeric_unitigs.tsv \
    --out ${OUT}/${SAMPLE}.hic_phased.sv.agp
python3 ${PIPE}/agp_to_fasta.py ${OUT}/${SAMPLE}.hic_phased.sv.agp 00_ref/utg.fasta ${OUT}/${SAMPLE}.hic_phased.sv.fasta

# The SV version must hold exactly the main assembly's unitigs, each once
cmp -s <(awk -F'\t' '$5=="W"{print $6}' ${MAIN} | sort) <(awk -F'\t' '$5=="W"{print $6}' ${OUT}/${SAMPLE}.hic_phased.sv.agp | sort) \
    || { echo "ERROR: the SV AGP does not contain the same unitigs as the main AGP" >&2; exit 1; }
python3 ${PIPE}/telomeres.py ${OUT}/${SAMPLE}.hic_phased.sv.fasta --motif ${TELOMERE_MOTIF} > ${OUT}/${SAMPLE}.hic_phased.sv.telomeres.tsv
python3 ${PIPE}/confidence.py ${OUT}/${SAMPLE}.hic_phased.sv.agp 10_final/membership.tsv 05_phase/phase.support.tsv \
    03_links/links.tsv > ${OUT}/${SAMPLE}.hic_phased.sv.confidence.tsv
awk 'NR>1{n+=gsub(/T/,"",$5); e+=2} END{print "# telomeric ends:", n, "of", e}' ${OUT}/${SAMPLE}.hic_phased.sv.telomeres.tsv
make_juicebox ${OUT}/${SAMPLE}.hic_phased.sv.agp ${OUT}/juicebox
echo "# Done: ${OUTDIR}/${OUT}"
