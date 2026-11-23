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

**A chain can be right and still need breaking.** The two-signal test is for chimeras. A different fault is a chain running through a fork, because a dropped unitig left the link looking one-to-one. It cannot be seen until groups exist, so `07_group` runs `scripts/review/chain_forks.py` for D06. Busk chain0061 had 705 kb at 1x then 787 kb at 2x, with Hi-C tied between homologs and the other exit leading to a chain on the homolog. Its reads are fine, since one haplotype does run that way. Break between arm and hub, then make the hub shared under D14. Do not act on a fork that rejoins its own group: that is an ordinary bubble. The scan's `break_at` column gives the position. The same fault with chains in no group shows only after placement: stage 10 repeats the scan as `12_shared/chain_forks.tsv` for D16's tangle check. There, GammaCitrus's chain0073 ran through hub utg000537l, and the other exit was the route to its unplaced 779 kb allele ([TANGLES](../TANGLES.md)).

**Edit:** `edits/chain_breaks.txt`, one per line: `chain position`. The position is the 1-based last base of the left part. The left part keeps its name; the right becomes `<chain>b`. A chain can take several rows, each position on the unbroken chain: the parts are `<chain>`, `<chain>b`, `<chain>c`, ... from left to right, so a hub can be cut out from both sides. An empty file is a valid outcome, and it was the outcome on WatsonsChoice.

**After:** `run` applies the breaks without remapping Hi-C.

## Also at D04: the Hi-C allele rescue

Read `09_hic_map/rescue_signal.txt`. Its Hi-C-blind chains have an allele partner, almost no unique Hi-C, and many pairs with an allele-ambiguous end: grouping will leave them unassigned however good the rest is. When they are a large share (the report's alarm: 1% of chain sequence or more), or carry many genes, record D04 with `--set hic_allele_rescue=true` and say why. This is the cheapest point: the setting acts in stages 06-08.

- Repeat chains are Hi-C-blind too (chain0091 on 292: 74% repeat). The rescue leaves chains over half repeat alone, so they do not count for it.
- The rescue knows which chromosome, rarely which haplotype: most of its haplotype calls are `arbitrary` or `allele_of_arbitrary`. It also brings orientation risk in the rescued region: on 292 inverted homolog blocks rose from 9 to 14, mostly in rescued chr2 sequence. Say so at D13 and hand those regions over for Juicebox.
- It leaves three kinds of chain unassigned, each named in `10_groups/groups.rescued.tsv`:
  - `allele_in_every_group`: an allele of 50 kb or more in every group of the chromosome.
  - `redundant`: allele blocks of any length cover half the chain or more in every group. 292's chain0288 was a 105 kb near-copy whose copies were each under 50 kb. Placed, it fell between a chain and its telomeric tip.
  - `collapsed`: a breakable 2x unique-read core (`10_groups/hub_cores.tsv`), so the chain holds both haplotypes there. That is for D06 to break; no haplotype call fits. A core inside one unitig cannot be broken off, so its chain is placed as usual.
  Collapsed and redundant chains still count towards their chromosome, so chains rescued in later rounds can anchor on them.
- Turned on later, it costs `reset 07_group` (from D05) or more.
