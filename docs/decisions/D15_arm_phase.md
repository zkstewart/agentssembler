# D15. Which arm goes where?

**Only for `collapsed: true`.** `assign_arms2.py` chooses arms jointly with depth targets and labels what settled each one (`12_shared/arm_report.tsv`):

| Label | Meaning |
|---|---|
| `only` | One candidate path |
| `coverage` | Every other path scores lower, alone and with a change at one other junction |
| `supported` | Some alternatives tie, and reads and Hi-C exceed each by the ratio |
| `arbitrary` | Some alternative scores the same and nothing separates them |

**Most arm phase is unknowable.** On WatsonsChoice 102 of 139 bridge junctions are arbitrary. Exchanging arms between haplotypes scores the same. Both alleles are present, but which haplotype holds which is not known. **Report it; do not paper over it.** An earlier test changing one junction at a time made most of them look decided; only paired swaps show this.

**Two checks:**

- A near-2x chain with a real alternative is the larger of two alleles, not one chain used twice. WatsonsChoice chain0398 and chain0564 are the example (`--minAllele`, `--maxPair`).
- A 1x chain on a path that is the only route through 2x flanks may legitimately be used twice. Say so. I first called such reuse legitimate without testing it, and ten were not clearly 2x: use `scripts/review/reused_depth.py` and `scripts/review/siblings.py`.

**Usually no edit is needed.** To record the outcome: `--reason "x only, y coverage, z supported, w arbitrary; reused arms checked."`
