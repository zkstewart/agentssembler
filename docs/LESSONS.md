# What the genomes taught

**Not needed to run the pipeline.** These are the principles behind the process, one or two sentences each. The numbers are in [RECORD](RECORD.md), under the section named; the rule that acts on each lesson is linked. They come from the citrus genomes in RECORD's roster, so treat them as evidence, not laws.

## The graph

- **Dropped does not mean useless.** A unitig set aside as redundant can be the only graph link between two chains, and restoring these connectors removed most gaps before any read filling. RECORD: Busk #3, "What rebuilding with connectors changed". Acts through `restore_connectors`.
- **A chain can run straight through a fork, and a haplotype's route is then lost.** Connectors are found only between chain ends, so the fix is a break at the hub. RECORD: "Lost-route forks". Acts through the tangle check ([TANGLES](TANGLES.md), [D16](decisions/D16_tails.md)).
- **A collapsed hub hides inside a chain at the chain's own median depth.** Depth of unique reads separates a true collapse (a 2x core beside 1x arms) from allele bubbles, and finding them before D04 saves remap rounds. RECORD: "Hi-C allele rescue and collapsed hubs (292)". Acts through [D04](decisions/D04_chimeric_chains.md) and [D06](decisions/D06_wrong_group.md).

## Hi-C

- **Sequence set aside must stay in the Hi-C target, as decoys.** Otherwise organelle reads land on nuclear insertions as if unique and merge chromosomes. RECORD: "Organelle Hi-C on nuclear insertions (13Q010)". Built into stage 05; the signal is in [D05](decisions/D05_group_count.md).
- **Near-identical haplotypes are near-identical to Hi-C too.** Hi-C can name the chromosome of such a chain but rarely its haplotype, and it cannot phase arms or place a near-exact copy. RECORD: "Phasing arms", "Allele-aware Hi-C for hard unplaced chains", "Hi-C allele rescue and collapsed hubs (292)". Acts through [D04](decisions/D04_chimeric_chains.md) and [D15](decisions/D15_arm_phase.md).
- **An identical overlap is a graph link, not an allele.** Near-identical haplotypes make de 0 overlaps common, and reading them as allele pairs puts a chain on the haplotype opposite the one it continues. RECORD: "Busk rebuilt with 0.3.6 and 0.3.7 (2026-11)". Acts through [D04](decisions/D04_chimeric_chains.md).
- **Where neither the graph nor Hi-C can tell alleles apart, the choice is the user's.** RECORD: "Hi-C allele rescue and collapsed hubs (292)". Acts through rescue swaps at [D16](decisions/D16_tails.md) and [D13](decisions/D13_user_decisions.md).

## Reads

- **Read correction loses and shortens reads near gaps, but it is not why gaps exist, and neither read set is enough alone.** RECORD: "Read correction", "Why reads stop at gaps". Acts through [D18](decisions/D18_raw_read_joins.md).
- **Read evidence per site, never pooled.** One array can supply nearly all of a pooled count. RECORD: Busk #16. Acts through [TRAPS](TRAPS.md) T13 and [D18](decisions/D18_raw_read_joins.md).

## Builds and results

- **A rebuild changes more than the gaps:** orientation between Hi-C-only gaps, which satellite stacks are placed, and haplotype numbers. RECORD: "What rebuilding with connectors changed". Acts through [D12](decisions/D12_sv_calls.md).
- **A low per-haplotype gene score can be numbering.** Haplotype numbers are arbitrary across chromosomes, so with translocations "all h1 scaffolds" is not one set. RECORD: "Gene sets across translocations (13Q010)". Acts through `validate/complement_genes.py`.
- **An edit is done when the built junction shows it.** An order link sets order, not a join; a group edit resolves against the groups before edits. RECORD: "Busk rebuilt with 0.3.6 and 0.3.7 (2026-11)". Acts through [D06](decisions/D06_wrong_group.md) and [D10](decisions/D10_curation.md).
- **An automated label is a claim, not a check.** Several were wrong in ways only the evidence behind them showed. RECORD: "Automated labels that were wrong".
- **A result that looks exact or too good is a question about the test.** Every significant fix began with a result that looked wrong while every command had succeeded. RECORD: "Errors made while building this". Acts through [AGENTS](../AGENTS.md).
