# D15. Which arm goes where?

**Only for `collapsed: true`.** `assign_arms2.py` chooses arms jointly with depth targets and labels what settled each one (`12_shared/arm_report.tsv`):

| Label | Meaning |
|---|---|
| `only` | One candidate path |
| `coverage` | Every other path scores lower, alone and with a change at one other junction |
| `supported` | Some alternatives tie, and reads and Hi-C exceed each by the ratio |
| `arbitrary` | Some alternative scores the same and nothing separates them |

**Most arm phase is unknowable.** On WatsonsChoice 102 of 139 bridge junctions are arbitrary. Exchanging arms between haplotypes scores the same. Both alleles are present, but which haplotype holds which is not known. **Report it; do not paper over it.** An earlier test changing one junction at a time made most of them look decided; only paired swaps show this.

**An arbitrary choice is repeatable, not right.** Among solutions with the same score, observations choose only when they decide (`--ratio`, `--minObs`); otherwise the first path by unitig names wins. Where two solutions have the same reads, as when homolog junctions swap arms, the reads cancel and Hi-C alone is compared, for the choice and its label alike. On 292 the old rule followed small count differences, so a `reset 06_evidence --remap` (new read mapping, renumbered connectors) moved two arms to the other haplotype with no new evidence. Compare arbitrary arms between builds by unitig, not by chain name.

**After raw reads, some of these are no longer unknown.** `17_raw_reads/phase_reads.tsv` counts reads linking an arm to specific sequence, on its own scaffold (`cis`) and on the homolog (`trans`). Report `confirmed` separately from `arbitrary`. On Busk reads spoke at 16 of 32 arbitrary junctions and on WatsonsChoice at 28 of 102, so most stay unknown: specific sequence is typically hundreds of kb from the arm. A `conflict` is a prompt to look. **No arm has been swapped on this evidence,** and no route for doing so has been built.

**Hi-C cannot help here.** On Busk 99.6% of Hi-C records on arm chains are exact ties between copies, with only 293 of 3.17 M showing any margin, so no filter or competitive mapping can recover them. Connectors have no Hi-C pairs at all, because they are restored after Hi-C mapping.

**Two checks:**

- A near-2x chain with a real alternative is the larger of two alleles, not one chain used twice. WatsonsChoice chain0398 and chain0564 are the example (`--minAllele`, `--maxPair`).
- A 1x chain on a path that is the only route through 2x flanks may legitimately be used twice. Say so. I first called such reuse legitimate without testing it, and ten were not clearly 2x: use `scripts/review/reused_depth.py` and `scripts/review/siblings.py`. Where a reused arm's likely allele is unplaced and reachable from the same hub, `12_shared/chain_forks.tsv` lists it as `LOST`: a chain ran through the fork. Settle it in D16's tangle check ([TANGLES](../TANGLES.md)).

**Usually no edit is needed.** To record the outcome: `--reason "x only, y coverage, z supported, w arbitrary; reused arms checked."`
