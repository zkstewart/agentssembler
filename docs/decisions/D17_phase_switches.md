# D17. Is a gap end a phase switch?

**Only for `collapsed: true`.**

**Read:** `scripts/review/cross_links.py`, then `scripts/validate/phase_scan.py` (`15_validate` writes `phase_scan.txt`).

**Links between homologs are normally overlap at divergent ends,** not phase switches. On WatsonsChoice `cross_links.py` found 21 such links, and Hi-C to current neighbours usually beat the alternative by 5x or more.

**How to decide:** compare Hi-C to the current neighbour with Hi-C to the alternative, then scan the ratio of homolog to far own contacts along the scaffold.

**Act only on** a clear and contiguous excess. Otherwise add the region to `16_validate/curation_targets.tsv` for the human review. The weakest WatsonsChoice cases, chain0003 (1.2x) and chain0302 (ambiguous), were left as ordered and listed there.

**At this stage the question is mostly a check.** Record `--reason "cross_links: N links; weakest ...; no switch shown."`. If a chain moves between haplotypes, it is a group edit: `reset 08_bridges --keepDecision`.
