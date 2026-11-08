#!/bin/bash
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Reference-guided Hi-C phasing and scaffolding of a polyploid hifiasm p_utg assembly.
# See REPLICATION.md for the method, inputs, decision points and parameters.
#
# Usage: run_pipeline.sh <config.sh> [stage ...]
#   stages (default all, in order):
#     map refplace links haphic phase validate ragtag relocate groups sort build qc juicebox dotplots
# A completed stage writes <stage>.done in OUTDIR and is skipped on rerun; delete it to redo.
# Contact maps: plot_contact_maps.sh <config.sh>. Structural variants: apply_sv.sh <config.sh>.
set -eo pipefail

[[ $# -ge 1 && -f $1 ]] || { echo "Usage: $0 <config.sh> [stage ...]" >&2; exit 1; }
PIPE=$(cd "$(dirname "$0")" && pwd)
CONFIG=$(readlink -f "$1"); shift
source "${PIPE}/lib.sh"
load_config "${CONFIG}"
STAGES=${@:-map refplace links haphic phase validate ragtag relocate groups sort build qc juicebox dotplots}
mkdir -p "${OUTDIR}"; cd "${OUTDIR}"
echo "# Output directory: ${OUTDIR}"

for stage in ${STAGES}; do
    if [[ -f ${stage}.done ]]; then echo "# ${stage}: already done"; continue; fi
    echo "# ${stage}: started $(date '+%F %T')"
case ${stage} in

map) # Hi-C to unitigs; keep allele-specific pairs (MAPQ >= 1, NM < 3) as HapHiC recommends
    mkdir -p 00_ref 01_map
    ln -sf "${UNITIGS}" 00_ref/utg.fasta
    use_env ${ENV_ASSEMBLY}; samtools faidx 00_ref/utg.fasta
    if [[ -n ${HIC_BAM} ]]; then
        ln -sf "${HIC_BAM}" 01_map/HiC.filtered.bam
    else
        [[ -f 00_ref/utg.fasta.bwt ]] || ${BWA} index 00_ref/utg.fasta > 00_ref/bwa_index.log 2>&1
        ${BWA} mem -5SP -t ${THREADS} 00_ref/utg.fasta "${HIC_R1}" "${HIC_R2}" 2> 01_map/bwa.log \
            | samblaster 2> 01_map/samblaster.log \
            | samtools view - -@ 4 -S -h -b -F 3340 -o 01_map/HiC.bam
        use_env ${ENV_HAPHIC}
        ${HAPHIC_DIR}/utils/filter_bam 01_map/HiC.bam 1 --nm 3 --threads ${THREADS} \
            | samtools view - -b -@ 4 -o 01_map/HiC.filtered.bam
    fi
    ;;

refplace) # Place unitigs on the reference: annotarium's coverage-based partition is the starting point
    mkdir -p 02_refplace
    if [[ -n ${PHASE_DIR} ]]; then
        ln -sf "${PHASE_DIR}/phasing.tsv" 02_refplace/phasing.tsv
        ln -sf "${PHASE_DIR}/minimap2.paf" 02_refplace/minimap2.paf
    else
        use_env ${ENV_BASE}
        python "${ANNOTARIUM_DIR}/annotarium.py" assembly phase -i 00_ref/utg.fasta -r "${REFERENCE}" \
            -o 02_refplace/phase --ploidy ${PLOIDY} --preset asm20 --threads ${THREADS} --minimap2 ${MINIMAP2}
        ln -sf phase/phasing.tsv 02_refplace/phasing.tsv; ln -sf phase/minimap2.paf 02_refplace/minimap2.paf
    fi
    ;;

links) # Hi-C pairs (read 1; '=' mate on the same unitig) and unitig-unitig link counts + depths
    mkdir -p 03_links; use_env ${ENV_ASSEMBLY}
    samtools view -@ 4 -f 64 01_map/HiC.filtered.bam \
        | awk -F'\t' -v OFS='\t' '{print $3, $4, ($7 == "=" ? $3 : $7), $8}' | gzip -1 > 03_links/pairs.tsv.gz
    samtools view -@ 4 01_map/HiC.filtered.bam | awk -F'\t' '
        { depth[$3]++ }
        and($2, 64) && $7 != "=" && $7 != "*" { if ($3 < $7) key = $3 "\t" $7; else key = $7 "\t" $3; links[key]++ }
        END { for (k in links) print k "\t" links[k]; for (c in depth) print c "\t.\t" depth[c] }' > 03_links/links.tsv
    awk -F'\t' '$1==$3{s++} END{printf "# %d allele-specific pairs; %.1f%% within one unitig\n", NR, 100*s/NR}' <(zcat 03_links/pairs.tsv.gz)
    ;;

haphic) # Reference-free HapHiC run: an independent check, and the inputs for Hi-C ordering later
    mkdir -p 04_haphic; use_env ${ENV_HAPHIC}
    (cd 04_haphic && ${HAPHIC_DIR}/haphic pipeline ../00_ref/utg.fasta ../01_map/HiC.filtered.bam ${NGROUPS} \
        --RE "${RE}" --remove_allelic_links ${PLOIDY} --correct_nrounds 2 \
        --threads ${THREADS} --processes ${THREADS} > haphic.log 2>&1)
    ;;

phase) # Reference-guided Hi-C phasing
    mkdir -p 05_phase; use_env ${ENV_PYTHON}
    python3 ${PIPE}/hic_phase.py --paf 02_refplace/minimap2.paf --phasing 02_refplace/phasing.tsv \
        --links 03_links/links.tsv --out 05_phase/phase --ploidy ${PLOIDY} --beam ${BEAM} \
        --minOverlap ${MIN_OVERLAP} --minOverlapFrac ${MIN_OVERLAP_FRAC}
    column -t 05_phase/phase.report.tsv
    ;;

validate) # Is the phasing signal real? Shuffled-link null and a planted-phasing simulation
    mkdir -p 05_phase/validate; use_env ${ENV_PYTHON}
    python3 ${PIPE}/null_links.py 02_refplace/phasing.tsv 03_links/links.tsv 05_phase/validate/null_links.tsv
    python3 ${PIPE}/hic_phase.py --paf 02_refplace/minimap2.paf --phasing 02_refplace/phasing.tsv \
        --links 05_phase/validate/null_links.tsv --out 05_phase/validate/null --ploidy ${PLOIDY} --beam ${BEAM} \
        --minOverlap ${MIN_OVERLAP} --minOverlapFrac ${MIN_OVERLAP_FRAC}
    paste <(cut -f1,9 05_phase/phase.report.tsv) <(cut -f9 05_phase/validate/null.report.tsv) \
        | sed '1s/.*/#chrom\tmean_support_real\tmean_support_null/' > 05_phase/validate/support_vs_null.tsv
    column -t 05_phase/validate/support_vs_null.tsv
    python3 ${PIPE}/simulate_phasing.py --paf 02_refplace/minimap2.paf --phasing 02_refplace/phasing.tsv \
        --ploidy ${PLOIDY} --workdir 05_phase/validate/simulation | tee 05_phase/validate/simulation.txt
    ;;

ragtag) # Reference order and orientation of each haplotype partition (annotarium + RagTag)
    use_env ${ENV_BASE}
    python "${ANNOTARIUM_DIR}/annotarium.py" assembly scaffold -i 05_phase/phase.phasing.tsv \
        -f 00_ref/utg.fasta -r "${REFERENCE}" -o 06_ragtag --threads ${THREADS} --minimap2 ${MINIMAP2} \
        --paf2dotplot ${PAF2DOTPLOT} --preset asm20 > 06_ragtag.log 2>&1
    ;;

relocate) # Unitigs whose Hi-C points to another chromosome
    mkdir -p 07_relocate; use_env ${ENV_PYTHON}
    python3 ${PIPE}/relocate_scan.py 05_phase/phase.phasing.tsv 03_links/links.tsv --paf 02_refplace/minimap2.paf \
        --ratio ${RELOCATE_RATIO} --minLinks ${RELOCATE_MIN_LINKS} --reliableQ ${RELOCATE_RELIABLE_Q} > 07_relocate/relocate_scan.tsv
    echo "# flagged: $(($(wc -l < 07_relocate/relocate_scan.tsv) - 1)) unitigs"
    ;;

groups) # Final (chromosome, haplotype) groups: phasing + Hi-C relocations; SV candidates listed
    mkdir -p 08_groups; use_env ${ENV_PYTHON}
    python3 ${PIPE}/make_groups.py --phasing 05_phase/phase.phasing.tsv --scan 07_relocate/relocate_scan.tsv \
        --fasta 04_haphic/01.cluster/corrected_asm.fa --clm 04_haphic/01.cluster/paired_links.clm \
        --minReliable ${MIN_RELIABLE} --RE "${RE}" --out 08_groups > 08_groups/groups_summary.txt
    cut -f3 08_groups/membership.tsv | sort | uniq -c
    ;;

sort) # Hi-C order of each group (HapHiC sort/build), used to place relocated unitigs
    mkdir -p 09_sort; use_env ${ENV_HAPHIC}
    (cd 09_sort && ${HAPHIC_DIR}/haphic sort ../04_haphic/01.cluster/corrected_asm.fa ../04_haphic/01.cluster/HT_links.pkl \
            ../08_groups/split_clms ../08_groups/groups/*.txt --processes ${THREADS} > sort.log 2>&1 \
        && ${HAPHIC_DIR}/haphic build ../04_haphic/01.cluster/corrected_asm.fa ../00_ref/utg.fasta ../01_map/HiC.filtered.bam \
            final_tours/*.tour --corrected_ctgs ../04_haphic/01.cluster/corrected_ctgs.txt > build.log 2>&1)
    ;;

build) # Pseudochromosomes: RagTag skeleton + Hi-C insertions + telomere capping
    mkdir -p 10_final; use_env ${ENV_PYTHON}
    python3 ${PIPE}/telomeres.py 00_ref/utg.fasta --motif ${TELOMERE_MOTIF} | awk -F'\t' 'NR>1 && $5!="--"' > 10_final/telomeric_unitigs.tsv
    python3 ${PIPE}/rescue_telomeric.py --membership 08_groups/membership.tsv --phasing 05_phase/phase.phasing.tsv \
        --telomeric 10_final/telomeric_unitigs.tsv --links 03_links/links.tsv \
        --minLinks ${RESCUE_MIN_LINKS} --minShare ${RESCUE_MIN_SHARE} --out 10_final/membership.tsv
    python3 ${PIPE}/build_final.py --membership 10_final/membership.tsv \
        --ragtag $(for h in $(seq 1 ${PLOIDY}); do echo 06_ragtag/partition$h/ragtag/ragtag.scaffold.agp; done) \
        --tours 09_sort/final_tours --fasta 00_ref/utg.fasta --phasing 05_phase/phase.phasing.tsv \
        --telomeric 10_final/telomeric_unitigs.tsv --out 10_final/${SAMPLE}.hic_phased > 10_final/build.log
    cp 08_groups/sv_candidates.tsv 10_final/
    ;;

qc) # Telomeres, phase-support confidence, and agreement of Hi-C order with the reference
    use_env ${ENV_PYTHON}
    python3 ${PIPE}/telomeres.py 10_final/${SAMPLE}.hic_phased.fasta --motif ${TELOMERE_MOTIF} > 10_final/${SAMPLE}.hic_phased.telomeres.tsv
    python3 ${PIPE}/confidence.py 10_final/${SAMPLE}.hic_phased.agp 10_final/membership.tsv \
        05_phase/phase.support.tsv 03_links/links.tsv > 10_final/${SAMPLE}.hic_phased.confidence.tsv
    python3 ${PIPE}/order_agreement.py 09_sort/scaffolds.agp 05_phase/phase.support.tsv > 10_final/hic_vs_reference_order.tsv
    python3 ${PIPE}/haphic_agreement.py 05_phase/phase.phasing.tsv 04_haphic/04.build/scaffolds.agp > 10_final/haphic_agreement.tsv
    echo "# phasing vs HapHiC: $(tail -1 10_final/haphic_agreement.tsv | cut -f4) of phased bp in each group's dominant HapHiC group"
    awk 'NR>1{n+=gsub(/T/,"",$5); e+=2} END{print "# telomeric ends:", n, "of", e}' 10_final/${SAMPLE}.hic_phased.telomeres.tsv
    tail -1 10_final/${SAMPLE}.hic_phased.confidence.tsv
    ;;

juicebox) make_juicebox 10_final/${SAMPLE}.hic_phased.agp 11_juicebox ;;

dotplots) make_dotplots 10_final/${SAMPLE}.hic_phased.fasta 12_dotplots ;;

*) echo "Unknown stage '${stage}'" >&2; exit 1 ;;
esac
    touch ${stage}.done
    echo "# ${stage}: finished $(date '+%F %T')"
done
