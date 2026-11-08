#!/bin/bash
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Hi-C contact maps for the pipeline's assemblies (see contact_map.py for the method).
#
# Usage: plot_contact_maps.sh <config.sh> [section ...]
#   sections: baseline phasing haphic_order main sv pairs genome summary   (default: all;
#             baseline needs BASELINE_DIR, sv/pairs need apply_sv.sh output, pairs needs SV_PAIRS)
# Output: OUTDIR/14_contact_maps/<section>/. Each per-chromosome run writes one map per
# chromosome (its homologs side by side), all_chroms.png (grid) and genome.png (all
# pseudochromosomes at 1 Mb, for contacts between chromosomes), plus crosstalk.tsv files.
set -eo pipefail

[[ $# -ge 1 && -f $1 ]] || { echo "Usage: $0 <config.sh> [section ...]" >&2; exit 1; }
PIPE=$(cd "$(dirname "$0")" && pwd)
source "${PIPE}/lib.sh"
load_config "$(readlink -f "$1")"; shift
SECTIONS=${@:-baseline phasing haphic_order main sv pairs genome summary}
cd "${OUTDIR}"
OUT=14_contact_maps; PAIRS=03_links/pairs.tsv.gz
MAIN=10_final/${SAMPLE}.hic_phased.agp; SV=13_sv/${SAMPLE}.hic_phased.sv.agp
mkdir -p ${OUT}; use_env ${ENV_HAPHIC}

# One group per chromosome: its homologs named by a pattern with {c} (chromosome) and {h} (haplotype)
per_chrom_groups() {
    for c in ${CHROMS}; do
        names=()
        for h in $(seq 1 ${PLOIDY}); do n=${1//\{c\}/$c}; names+=("${n//\{h\}/$h}"); done
        printf "%s\t%s\n" "$c" "$(IFS=,; echo "${names[*]}")"
    done
}
# Per-chromosome maps, grid and whole-genome map: $1 = output dir, $2 = name pattern, rest = contact_map.py args
chrom_maps() {
    local dir=$1 pattern=$2; shift 2
    mkdir -p ${dir}/genome
    per_chrom_groups "${pattern}" > ${dir}/groups.tsv
    python3 ${PIPE}/contact_map.py --pairs ${PAIRS} --groups ${dir}/groups.tsv --out ${dir} "$@"
    python3 ${PIPE}/combine_images.py grid ${dir}/all_chroms.png $(for c in ${CHROMS}; do echo ${dir}/$c.png; done)
    printf "genome\t%s\n" "$(cut -f2 ${dir}/groups.tsv | paste -sd,)" > ${dir}/genome/groups.tsv
    python3 ${PIPE}/contact_map.py --pairs ${PAIRS} --groups ${dir}/genome/groups.tsv --out ${dir}/genome --bin 1000000 "$@"
    mv ${dir}/genome/genome.png ${dir}/genome.png
    echo "# ${dir}: $(awk 'NR>1{h+=$4; i+=$3} END{printf "homolog crosstalk %.3f", h/(h+i)}' ${dir}/crosstalk.tsv)"
}
partition_agps() { for h in $(seq 1 ${PLOIDY}); do echo "$1/partition$h/ragtag/ragtag.scaffold.agp"; done; }
partition_prefixes() { for h in $(seq 1 ${PLOIDY}); do echo "p${h}_"; done; }

for section in ${SECTIONS}; do
case ${section} in

baseline) # An earlier annotarium phase + RagTag run, if given
    if [[ -z ${BASELINE_DIR} ]]; then echo "# baseline: no BASELINE_DIR, skipped"; continue; fi
    chrom_maps ${OUT}/baseline 'p{h}_{c}_RagTag' --agp $(partition_agps ${BASELINE_DIR}) --prefixes $(partition_prefixes)
    ;;

phasing) # The Hi-C phasing before relocation and Hi-C ordering, in RagTag order
    chrom_maps ${OUT}/phasing 'p{h}_{c}_RagTag' --agp $(partition_agps 06_ragtag) --prefixes $(partition_prefixes)
    ;;

haphic_order) # The final groups ordered by Hi-C alone
    chrom_maps ${OUT}/haphic_order '{c}h{h}' --agp 09_sort/scaffolds.agp
    ;;

main) chrom_maps ${OUT}/main '{c}_h{h}' --agp ${MAIN} ;;

sv)
    if [[ ! -f ${SV} ]]; then echo "# sv: no ${SV} (run apply_sv.sh), skipped"; continue; fi
    chrom_maps ${OUT}/sv '{c}_h{h}' --agp ${SV}
    ;;

pairs) # Chromosome pairs involved in SVs, before (main) and after (sv), at 250 kb
    if [[ ! -f ${SV} || -z ${SV_PAIRS} ]]; then echo "# pairs: needs apply_sv.sh output and SV_PAIRS, skipped"; continue; fi
    mkdir -p ${OUT}/pairs/before ${OUT}/pairs/after
    for pair in ${SV_PAIRS}; do
        a=${pair%:*}; b=${pair#*:}; name=${a}_${b}
        printf "%s\t%s\n" ${name} "$(for x in $a $b; do for h in $(seq 1 ${PLOIDY}); do printf "${x}_h${h},"; done; done | sed 's/,$//')" \
            > ${OUT}/pairs/${name}.groups.tsv
        for state in before after; do
            agp=$([[ ${state} == before ]] && echo ${MAIN} || echo ${SV})
            python3 ${PIPE}/contact_map.py --pairs ${PAIRS} --agp ${agp} --groups ${OUT}/pairs/${name}.groups.tsv \
                --out ${OUT}/pairs/${state} --bin 250000
            mv ${OUT}/pairs/${state}/crosstalk.tsv ${OUT}/pairs/${state}/${name}.crosstalk.tsv
        done
        python3 ${PIPE}/combine_images.py row ${OUT}/pairs/${name}.before_after.png \
            ${OUT}/pairs/before/${name}.png ${OUT}/pairs/after/${name}.png
    done
    ;;

genome) # Whole-genome maps of main and sv side by side
    [[ -f ${OUT}/main/genome.png ]] || { echo "# genome: run the main section first, skipped"; continue; }
    if [[ -f ${OUT}/sv/genome.png ]]; then
        python3 ${PIPE}/combine_images.py row ${OUT}/main_vs_sv.genome.png ${OUT}/main/genome.png ${OUT}/sv/genome.png
    fi
    ;;

summary) # Homolog crosstalk of every per-chromosome run, overall and per chromosome
    {
        printf "#run\tall"; for c in ${CHROMS}; do printf "\t%s" $c; done; echo
        for f in $(find ${OUT} -name crosstalk.tsv | sort); do
            first=$(echo ${CHROMS} | awk '{print $1}')
            awk -F'\t' -v c=${first} 'NR == 2 { found = ($1 == c) } END { exit !found }' $f || continue # per-chromosome runs only
            awk -F'\t' -v run=$(dirname ${f#${OUT}/}) -v chroms="${CHROMS}" 'NR>1{h[$1]+=$4; i[$1]+=$3; H+=$4; I+=$3}
                END{printf "%s\t%.3f", run, H/(H+I); n=split(chroms, cs, " "); for(k=1;k<=n;k++){c=cs[k]; printf "\t%.3f", h[c]/(h[c]+i[c])} print ""}' $f
        done
    } > ${OUT}/crosstalk_summary.tsv
    column -t ${OUT}/crosstalk_summary.tsv
    ;;

*) echo "Unknown section '${section}'" >&2; exit 1 ;;
esac
done
