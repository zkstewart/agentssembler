# D16. Should a tail be placed?

A tail is a graph path leaving one flank and stopping before the other. It goes beside the flank it links to, with the gap beyond. **The default is not to place.** Nothing is placed without a row in `edits/accepted_tails.tsv`.

**First, the tangle check** ([TANGLES](../TANGLES.md)). Do it before choosing any tail, because what it finds replaces tails:

1. **Lost-route forks** (collapsed samples): read the `LOST` rows of `12_shared/chain_forks.tsv`. Break where the lost chain is the allele of a reused arm (`over` in `12_shared/arm_usage.tsv`), or leads on to more unplaced 1x or telomeric chains. Leave a bubble arm that rejoins the same chain. `BUBBLE` rows (the chain reads under 2x over the hub) need nothing. A break goes in `edits/chain_breaks.txt` at `break_at`, then `reset 06_evidence --remap`.
2. **Small tangles:** a group end with no telomere, whose tails stop inside a small tangle short of unplaced telomeric chains. Use TTT ([TANGLES](../TANGLES.md)), then reset from the earliest stage its edits touch.

Record the check in the reason even when nothing is done, for example: `Tangle check: 3 LOST rows, 1 broken (chain0073 -> chain0081, reused arm); TTT on chr5 (17 unitigs), both ends now telomeric; 2 left as bubbles.` After a reset you return here on a new pre-build, and choose tails then.

**Read:** `12_shared/tails.tsv`, `12_shared/tail_report.tsv`, `12_shared/tail_checks.tsv`.

**`make_tails.py` labels are not safety checks.** It tests a path's copy limit, repeat class and chromosome, per chain. On 13Q010, chain0256 was labelled `only`, and all three chr5 haplotypes already hold it at 92-99.9%. Placing it would make a fourth copy in a triploid.

**Run these tests on each arm** (`tail_checks.tsv` holds the first three):

1. **Is it already placed?** A near-exact copy on the receiving scaffold is a false duplication (13Q010 chain0269: 81% an exact copy of the scaffold it would join). Copies on every haplotype of the chromosome mean no free slot.
2. **Is the budget shared?** Near-identical candidates draw on one budget. Sum their depths and compare with haplotype depth. 13Q010 chain0207 and chain0214 are 100% identical over 85% and sum to about one copy.
3. **Is it an allele?** In a set, alleles align with each other at the sample's usual divergence over most of their length. Arms identical over most of their length are the same sequence, unless summed depth shows more than one copy.
4. **Do reads cross the join?** Measured after placement (`16_validate/tail_joins.tsv`). No MAPQ >= 20 reads means the join rests on the graph link alone.
5. **Does it add anything?** Count genes missing from the receiving haplotype, as well as genes on the arm.

**A full set with unknown assignment.** If as many arms as the ploidy match one chromosome and pass tests 1-3, place them, labelled arbitrary: the chromosome-level content is right whichever arm goes where. Where the flank is collapsed, there was no phase information to lose. Record it, and end haplotype blocks there.

**Risks to weigh, besides mosaicism:**

- False duplication inflates length and duplicated gene counts, and can look like copy-number variation. This is the largest.
- Thin position evidence: one overlap, often little Hi-C (5 of 8 13Q010 tails have one pair or none).
- An arbitrary placement looks as trustworthy as a phased one downstream, so haplotype comparisons can report assignment choices as differences.
- No gain: on 13Q010 none of the 7 genes on tails was missing from the haplotype that would receive it.

**This is the user's call when in doubt.** On 13Q010, the user chose not to place the six that passed, to minimise assembly risk. On WatsonsChoice, 47 were placed (7.06 Mb distinct, 27 arbitrary), none had a copy of 30% or more on its own scaffold within 2 Mb, and all joins had reads, 7 of them with none at q20.

**Edit:** copy accepted rows from `12_shared/tails.tsv` to `edits/accepted_tails.tsv` unchanged, and give the reasoning in the decision. Then `run`.

**A tail through a tangle the check resolved is superseded.** On GammaCitrus chr5, tails reached into the tangle, and TTT later carried both haplotypes through it to their telomeres. Accept the TTT rows, not those tails.

**Beyond tails:** placing them changes every coordinate downstream, so decide before mapping in `12_fill`.

**Rescue swaps** (`12_shared/rescue_swaps.tsv`, with `hic_allele_rescue`): Hi-C gaps in the pre-build that would close if two chains swapped homologs, because the graph links the pair the other way from the rescue's Hi-C call. Neither side is proof in near-identical haplotypes, so the choice is the user's. **Ask the user now, at D16, not at D13**: before the build a swap costs a group edit and `reset 08_bridges` (stages 08-10, minutes); after it, the same swap also reruns stages 11-15 (fill, raw reads and validation, about 4 h on 292). For each row give: which links are direct graph overlaps and which pass through a connector, the junction reads on each side, and the rescue's call and ratio. On 292 the chr3_h1 swap was deferred to D13, and the user then made it, at the cost of a second build. Record the answer, quoted, in D16; if the user defers, carry the row to D13.

## Genes on tails

`12_shared/tail_genes.tsv` (with genes on the unitigs) gives each tail's genes, and how many have no full set of placed copies elsewhere (`short`). A tail with short genes is worth testing; the genes do not pass it. On 292 the user chose the five gene-carrying tails at D13 for completeness; two had strong raw-read support (134 and 104 reads at q20) and three little or none. A tail left unplaced keeps its genes in the assembly as unlocalized sequence when its chromosome is known (stage 15), so the gain from placing it is position, not genes.

`12_shared/unplaced_genes.tsv` does the same for every chain left unplaced and on no tail, sorted by `short`, with `lost_from` naming the placed chain whose LOST fork leads to it. Use it in the tangle check: a LOST row or tangle whose chains carry short genes is worth the work, one carrying none rarely is. On the 292 rebuild, 54 of 99 such chains carried genes (178 short); the one LOST target, chain0266, carried 1.

**Direction, at this stop:** the report flags a pre-build scaffold written in reverse against its homolog (`12_shared/homolog_direction.tsv`). List it in `edits/flip.txt` now: this is the last stop before stage 12 maps reads ([D12](D12_sv_calls.md)).
