# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Shared functions for run_pipeline.sh, apply_sv.sh and plot_contact_maps.sh (source, don't run).

# Load a config file, resolve its relative paths against the config's directory, and derive
# CHROMS (natural order), NGROUPS (chromosomes x ploidy) and the environment the scripts need
load_config() {
    local config=$1 dir
    dir=$(dirname "$(readlink -f "${config}")")
    source "${config}"
    for var in UNITIGS GFA READS HIC_R1 HIC_R2 REFERENCE OUTDIR HIC_BAM PHASE_DIR BASELINE_DIR SV_DIR; do
        local val=${!var}
        if [[ -n ${val} && ${val} != /* ]]; then
            printf -v ${var} '%s' "$(realpath -m "${dir}/${val}")"
        fi
    done
    source "${CONDA_SH}"
    if [[ -z ${CHROMS} ]]; then
        [[ -f ${REFERENCE}.fai ]] || { use_env ${ENV_ASSEMBLY}; samtools faidx "${REFERENCE}"; }
        CHROMS=$(cut -f1 "${REFERENCE}.fai" | sort -V | tr '\n' ' ')
    fi
    NGROUPS=$(( $(echo ${CHROMS} | wc -w) * PLOIDY ))
    export ANNOTARIUM_DIR
}

use_env() { conda activate "$1"; }

# Juicebox (JBAT) files for an AGP over the unitigs: $1 = AGP, $2 = output dir,
# [$3 = .fai of the AGP's components when it has more than the unitigs, e.g. read fills].
# Load out_JBAT.hic together with out_JBAT.assembly in Juicebox for manual curation.
make_juicebox() {
    local agp=$(readlink -f "$1") out=$2 U=${HAPHIC_DIR}/utils fai=$(readlink -f "${3:-${OUTDIR}/00_ref/utg.fasta.fai}")
    mkdir -p "${out}"; use_env ${ENV_ASSEMBLY}
    (cd "${out}" \
        && ${U}/juicer pre -a -q 1 -o out_JBAT "${OUTDIR}/01_map/HiC.filtered.bam" "${agp}" "${fai}" > out_JBAT.log 2>&1 \
        && java -Djava.awt.headless=true -Xmx${JAVA_MEM} -jar ${U}/juicer_tools.1.9.9_jcuda.0.8.jar pre out_JBAT.txt out_JBAT.hic.part \
            <(grep PRE_C_SIZE out_JBAT.log | awk '{print $2" "$3}') > juicer_tools.log 2>&1 \
        && mv out_JBAT.hic.part out_JBAT.hic && rm -f out_JBAT.txt)
}

# Dotplot of each haplotype set (names ending _h<N>) against the reference: $1 = FASTA, $2 = output dir
make_dotplots() {
    local fasta=$1 out=$2
    mkdir -p "${out}"; use_env ${ENV_ASSEMBLY}
    for h in $(seq 1 ${PLOIDY}); do
        seqkit grep -r -p "_h${h}\$" "${fasta}" > "${out}/h${h}.fasta" 2> /dev/null
        ${MINIMAP2} -x asm20 -t $(( THREADS / 2 )) "${REFERENCE}" "${out}/h${h}.fasta" > "${out}/h${h}.paf" 2> "${out}/h${h}.minimap2.log"
        ${PAF2DOTPLOT} -b --plot-size=15 --min-ref-len=1000000 --min-query-length=400000 --min-alignment-length=10000 \
            -o "${out}/h${h}.vs_reference" "${out}/h${h}.paf" > /dev/null 2>&1
        rm -f "${out}/h${h}.fasta"
    done
}
