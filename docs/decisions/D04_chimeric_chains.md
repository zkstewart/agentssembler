# D04. Is a chain chimeric?

**Read:** the stage report, which crosses the two signals, then `10_groups/chain_breaks.tsv`. For candidates use `scripts/chain_profile.py` and `07_chains2/joins.verify.tsv`.

**A break needs at least two independent signals:**

1. Different Hi-C partners on each side (`partner_cos` well below 1), **and**
2. Weak reads across the join: few at MAPQ >= 20, with low local depth.

**Do not act on `cross_ratio` alone.** It is 0 at many satellite positions, whose flanks have no mappable Hi-C. It is also 0 where a tail of unique-looking sequence is near-identical to a sibling, so that Hi-C reads there get MAPQ 0. Read the `trans` counts on both sides first.

**Examples:**

- **Broken.** 13Q010 chain0010 at 15.78 Mb: `partner_cos` 0.48, no MAPQ >= 20 reads (8 at any MAPQ), low depth, and the two sides contact different chromosomes. The 1.65 Mb chain0010b was a chr9 piece on a chr6 chain.
- **Not broken.** 13Q010 chain0005 and chain0016 top the `cross_ratio` list, but they are satellite. On WatsonsChoice six chains showed `cross_ratio` 0 with a 100-170 kb tail, and no break was justified.

- **Not broken.** Busk chain0037 at 5.06 and 5.16 Mb: `partner_cos` 0.43-0.48 with plenty of `trans`, but 46-54 reads at MAPQ >= 20. Reads overrule Hi-C. Busk chain0021 at 5.91 Mb had 44 reads, 2 at MAPQ >= 20, `partner_cos` 0.71, and no clear change of partner, so no break either.

**A chain can be right and still need breaking.** The two-signal test is for chimeras. A different fault is a chain running through a fork, because a dropped unitig left the link looking one-to-one. It cannot be seen until groups exist, so `07_group` runs `scripts/review/chain_forks.py` for D06. Busk chain0061 had 705 kb at 1x then 787 kb at 2x, with Hi-C tied between homologs and the other exit leading to a chain on the homolog. Its reads are fine, since one haplotype does run that way. Break between arm and hub, then make the hub shared under D14. Do not act on a fork that rejoins its own group: that is an ordinary bubble.

**Edit:** `edits/chain_breaks.txt`, one per line: `chain position`. The position is the 1-based last base of the left part. The left part keeps its name; the right becomes `<chain>b`. An empty file is a valid outcome, and it was the outcome on WatsonsChoice.

**After:** `run` applies the breaks without remapping Hi-C.
