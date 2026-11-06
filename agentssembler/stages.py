"""The stage list. A gate lists decisions that must be recorded before later stages run."""
from dataclasses import dataclass, field

@dataclass
class Stage:
    id: str
    title: str
    gates: list = field(default_factory=list)          # decisions, each (id, condition or None)
    outputs: list = field(default_factory=list)        # paths moved aside by `reset`, kept when a mapping
    keep: list = field(default_factory=list)           # costly mappings, kept by `reset` unless --remap
    when: str = None                                   # setting that must be true for the stage to do anything
    hours: str = ""                                    # rough cost on 16 cores, from the two genomes tried

STAGES = [
    Stage("01_inputs", "Link and check inputs", outputs=[]),
    Stage("02_map_unitigs", "Reads and unitigs aligned", keep=["02_ont_map", "03_unitig_align"], hours="2-6 h; can run on an HPC"),
    Stage("03_depth", "Depth per unitig", outputs=["04_copynumber"], gates=[("D01", None)]),
    Stage("04_chains", "Non-redundant unitigs and chains", outputs=["05_nonredundant", "06_chains/chains*", "07_chains2/kept2.txt", "08_final_chains"], keep=["06_chains/*.paf.gz", "07_chains2/*.paf.gz"], gates=[("D03", None)], hours="1-3 h"),
    Stage("05_hic_chains", "Hi-C to chains, chimera table", outputs=["10_groups/chain_breaks.tsv"], keep=["09_hic_map"], gates=[("D04", None)], hours="about 1.6 h"),
    Stage("06_evidence", "Per-chain evidence", outputs=["12_broken", "10_groups/alleles*", "10_groups/chain_*", "10_groups/shared.candidates.txt", "16_validate/telomeres.chains.tsv", "10_groups/end_links.tsv", "11_ont_chains/chain_class.tsv", "11_ont_chains/depth*"], keep=["11_ont_chains/*.paf.gz"], hours="about 40 min"),
    Stage("07_group", "Grouping", outputs=["10_groups/groups*", "10_groups/hic_group.log", "10_groups/shared.txt", "10_groups/eval_groups.txt"], gates=[("D05", None), ("D06", None), ("D14", "collapsed")]),
    Stage("08_bridges", "Shared copies, bridges", outputs=["10_groups/final*", "12_shared", "13_scaffold/bridges"], gates=[("D07", None)], hours="minutes"),
    Stage("09_order", "Ordering per chromosome", outputs=["12_shared/order*", "12_shared/by_chrom", "12_shared/eval_order.txt", "12_shared/unhonoured.txt"], gates=[("D09", None), ("D10", None)], hours="minutes to 70 min"),
    Stage("10_arms_tails", "Arms and tail candidates", outputs=["12_shared/choices.tsv", "12_shared/arm_*", "12_shared/pre", "12_shared/tails.tsv", "12_shared/tail_*", "12_shared/build_args.txt", "13_scaffold/names.tsv", "13_scaffold/repeat_chains.fasta"], gates=[("D17", "collapsed"), ("D15", "collapsed"), ("D16", None)]),
    Stage("11_build", "Scaffolds", outputs=["13_scaffold/*.graph.*", "13_scaffold/target.fasta*", "13_scaffold/build.log"]),
    Stage("12_fill", "Read filling and final files", outputs=["15_final", "14_readfill/*/fill*", "14_readfill/merged.agp", "14_readfill/fills.fasta", "14_readfill/*.graph.*"], keep=["14_readfill/*.paf.gz"], gates=[("D11", None)], hours="about 1 h"),
    Stage("13_rawreads", "Raw reads at gaps, join audit", outputs=["17_raw_reads/fill", "17_raw_reads/cmp*", "17_raw_reads/audit.tsv", "17_raw_reads/neargap*"], keep=["17_raw_reads/raw_vs_target*", "17_raw_reads/reads.tsv.gz"], when="raw", gates=[("D18", "raw")], hours="1-2 h"),
    Stage("14_apply_joins", "Apply join decisions", outputs=["15_final_raw"], when="raw"),
    Stage("15_validate", "Validation on final files", outputs=["16_validate/maps", "16_validate/sv", "16_validate/juicebox", "16_validate/compleasm", "16_validate/reads", "16_validate/pairs*", "16_validate/*.txt", "16_validate/depth*", "16_validate/telomeres.final.tsv"], gates=[("D12", None), ("D13", None)], hours="about 1.5 h"),
]
BY_ID = {s.id: s for s in STAGES}

GUIDES = {
    "D01": "D01_haplotype_depth.md", "D03": "D03_weak_joins.md", "D04": "D04_chimeric_chains.md", "D05": "D05_group_count.md", "D06": "D06_wrong_group.md",
    "D07": "D07_telomere_tips.md", "D09": "D09_order.md", "D10": "D10_curation.md", "D11": "D11_gap_closures.md", "D12": "D12_sv_calls.md", "D13": "D13_user_decisions.md",
    "D14": "D14_shared_chains.md", "D15": "D15_arm_phase.md", "D16": "D16_tails.md", "D17": "D17_phase_switches.md", "D18": "D18_raw_read_joins.md",
}
