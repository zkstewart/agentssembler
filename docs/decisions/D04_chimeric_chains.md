# D04. Is a chain chimeric?

**Read:** the stage report, which crosses the two signals, then `10_groups/chain_breaks.tsv`. For candidates use `scripts/chain_profile.py` and `07_chains2/joins.verify.tsv`.

**A break needs at least two independent signals:**

1. Different Hi-C partners on each side (`partner_cos` well below 1), **and**
2. Weak reads across the join: few at MAPQ >= 20, with low local depth.

**Do not act on `cross_ratio` alone.** It is 0 at many satellite positions, whose flanks have no mappable Hi-C. It is also 0 where a tail of unique-looking sequence is near-identical to a sibling, so that Hi-C reads there get MAPQ 0. Read the `trans` counts on both sides first.

**Examples:**

- **Broken.** A join with `partner_cos` near 0.5, no MAPQ >= 20 reads (8 at any MAPQ), low depth, and the two sides contacting different chromosomes. The part cut off, over 1 Mb, was a piece of another chromosome.
- **Not broken.** Chains at the top of the `cross_ratio` list that are satellite. Chains with `cross_ratio` 0 at a 100-170 kb tail, where no break was justified.
- **Not broken.** Joins with `partner_cos` under 0.5 and plenty of `trans`, but about 50 reads at MAPQ >= 20. Reads overrule Hi-C.
- **Not broken.** A join with over 40 reads, only 2 of them at MAPQ >= 20, `partner_cos` about 0.7, and no clear change of partner.

**A chain can be right and still need breaking.** The two-signal test is for chimeras. A different fault is a chain running through a fork, because a dropped unitig left the link looking one-to-one. It cannot be seen until groups exist, so `07_group` runs `scripts/review/chain_forks.py` for D06. One such chain had about 0.7 Mb at 1x then about 0.8 Mb at 2x, with Hi-C tied between homologs and the other exit leading to a chain on the homolog. Its reads are fine, since one haplotype does run that way. Break between arm and hub, then make the hub shared under D14. Do not act on a fork that rejoins its own group: that is an ordinary bubble.

**Edit:** `edits/chain_breaks.txt`, one per line: `chain position`. The position is the 1-based last base of the left part. The left part keeps its name; the right becomes `<chain>b`. An empty file is a valid outcome.

**After:** `run` applies the breaks without remapping Hi-C.
