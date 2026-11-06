# D16. Should a tail be placed?

A tail is a graph path leaving one flank and stopping before the other. It goes beside the flank it links to, with the gap beyond. **The default is not to place.** Nothing is placed without a row in `edits/accepted_tails.tsv`.

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

**Beyond tails:** placing them changes every coordinate downstream, so decide before mapping in `12_fill`.
