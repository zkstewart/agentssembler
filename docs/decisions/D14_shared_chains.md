# D14. Is a chain shared between haplotypes?

**Only for `collapsed: true`.** A shared chain is collapsed sequence that belongs in every haplotype of its chromosome. `make_shared_copies.py` writes a copy `<chain>.<group>` per haplotype.

**Read:** `11_ont_chains/chain_class.tsv`, `10_groups/shared.candidates.txt`, `10_groups/groups.shared.tsv`, and `scripts/review/hub_depth.py`, `scripts/review/reused_depth.py`.

**Yes when** interior depth (outside graph overlaps) is near 2x and Hi-C is tied between homologs.

**Not shared when** the chain is high copy: over about 3x, with Hi-C across all groups. Unassign it (WatsonsChoice chain0184, 262x, and chain0267, 100x).

**The score ratio alone is too strict.** Hub chains at 2x with a ratio of 1.5 were left unassigned or forced into one haplotype on WatsonsChoice. Read the call with depth and graph position. 11 chains were promoted by hand, for example:

`chain0175  shared  c8  Interior depth 64.6x; Hi-C to chr3 5.23x next chromosome (32,356 pairs), no separation between g17 and g18; 13 complete genes; graph-linked hub of the chr3 gap`

On Busk five short hubs (190-340 kb) had medians of 71-84x, depth class `specific`, with 50-130 kb overlaps at each end and 125-145x in between. The auto-candidates missed them, because class comes from the median. The report lists such chains from the middle third of their windows, **which only nominates them**.

**Then run `scripts/review/hub_depth.py`, and require real interior.** I first promoted seven on window depth, then reversed two: chain0171 and chain0206 had 4.7 and 4.4 kb outside graph overlaps, and about 360 Hi-C pairs. A depth over that little sequence says nothing. The five kept had 27-119 kb. Unassigned chains can still serve as bridge arms, so leaving one out loses little. Read degree too: chain0108 is tied in Hi-C but has degree 1,1 and depth about 60-88x, and was left alone.

**Why this matters:** promotion and 8-hop ordering, with tails, took genes complete only in unplaced sequence from 52 to 0. They were tested together, so their separate contributions are not known.

**Candidate list:** by default it comes from depth class. To replace it, write `edits/shared.txt` (one chain per line, no comments) before `07_group` runs, or reset that stage.

**Edit:** `edits/group_edits.tsv`, action `shared` with the chromosome label from `10_groups/groups.chroms.tsv`.
