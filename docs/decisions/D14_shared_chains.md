# D14. Is a chain shared between haplotypes?

**Only for `collapsed: true`.** A shared chain is collapsed sequence that belongs in every haplotype of its chromosome. `make_shared_copies.py` writes a copy `<chain>.<group>` per haplotype.

**Read:** `11_ont_chains/chain_class.tsv`, `10_groups/shared.candidates.txt`, `10_groups/groups.shared.tsv`, and `scripts/review/hub_depth.py`, `scripts/review/reused_depth.py`.

**Yes when** interior depth (outside graph overlaps) is near 2x and Hi-C is tied between homologs.

**Not shared when** the chain is high copy: over about 3x, with Hi-C across all groups. Unassign it (chains at 4 and 9 times haplotype depth, for example).

**The score ratio alone is too strict.** Hub chains at 2x with a ratio of 1.5 were left unassigned or forced into one haplotype on one genome. Read the call with depth and graph position, and promote by hand where they agree, with a reason of this kind:

`<chain>  shared  with:<chain>  Interior depth 2.3x haplotype; Hi-C to this chromosome 5x the next chromosome, no separation between its two groups; 13 complete genes; graph-linked hub of a gap`

Short hubs of a few hundred kb can have a median of 1.1-1.4x haplotype depth and depth class `specific`, with 50-130 kb overlaps at each end and 2x or more in between. The auto-candidates missed them, because class comes from the median. The report lists such chains from the middle third of their windows, **which only nominates them**.

**Then run `scripts/review/hub_depth.py`, and require real interior.** I once promoted chains on window depth, then reversed two of them: each had under 5 kb outside graph overlaps, and a few hundred Hi-C pairs. A depth over that little sequence says nothing. Those kept had tens of kb or more. Unassigned chains can still serve as bridge arms, so leaving one out loses little. Read degree too: a chain tied in Hi-C, but with degree 1,1 and depth of 1-1.4x haplotype, was left alone.

**Why this matters:** promotion and 8-hop ordering, with tails, took genes complete only in unplaced sequence from 52 to 0. They were tested together, so their separate contributions are not known.

**Candidate list:** by default it comes from depth class. To replace it, add rows to `edits/shared.txt` (one chain per line, under its header) before `07_group` runs, or reset that stage. With no rows, the automatic list is used.

**Edit:** `edits/group_edits.tsv`, action `shared`, value `with:<a chain on that chromosome>`.
