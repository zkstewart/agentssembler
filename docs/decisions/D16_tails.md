# D16. Should a tail be placed?

A tail is a graph path leaving one flank and stopping before the other. It goes beside the flank it links to, with the gap beyond. **The default is not to place.** Nothing is placed without a row in `edits/accepted_tails.tsv`.

**Read:** `12_shared/tails.tsv`, `12_shared/tail_report.tsv`, `12_shared/tail_checks.tsv`.

**`make_tails.py` labels are not safety checks.** It tests a path's copy limit, repeat class and chromosome, per chain. One chain labelled `only` was already held by every haplotype of its chromosome, at 92-99.9% identity. Placing it would have made one copy more than the ploidy.

**Run these tests on each arm** (`tail_checks.tsv` holds the first three):

1. **Is it already placed?** A near-exact copy on the receiving scaffold is a false duplication (one candidate was 81% an exact copy of the scaffold it would join). Copies on every haplotype of the chromosome mean no free slot.
2. **Is the budget shared?** Near-identical candidates draw on one budget. Sum their depths and compare with haplotype depth. Two candidates that were 100% identical over 85% of their length summed to about one copy.
3. **Is it an allele?** In a set, alleles align with each other at the sample's usual divergence over most of their length. Arms identical over most of their length are the same sequence, unless summed depth shows more than one copy.
4. **Do reads cross the join?** Measured after placement (`16_validate/tail_joins.tsv`). No MAPQ >= 20 reads means the join rests on the graph link alone.
5. **Does it add anything?** Count genes missing from the receiving haplotype, as well as genes on the arm.

**A full set with unknown assignment.** If as many arms as the ploidy match one chromosome and pass tests 1-3, place them, labelled arbitrary: the chromosome-level content is right whichever arm goes where. Where the flank is collapsed, there was no phase information to lose. Record it, and end haplotype blocks there.

**Risks to weigh, besides mosaicism:**

- False duplication inflates length and duplicated gene counts, and can look like copy-number variation. This is the largest.
- Thin position evidence: one overlap, often little Hi-C (on one genome most tails had one pair or none).
- An arbitrary placement looks as trustworthy as a phased one downstream, so haplotype comparisons can report assignment choices as differences.
- No gain: on one genome none of the genes on tails was missing from the haplotype that would receive it.

**This is the user's call when in doubt.** Give the user the tails that passed, with their test results, the distinct sequence they would add, and how many placements would be arbitrary.

**Edit:** copy accepted rows from `12_shared/tails.tsv` to `edits/accepted_tails.tsv` unchanged, and give the reasoning in the decision. Then `run`.

**Beyond tails:** placing them changes every coordinate downstream, so decide before mapping in `12_fill`.
