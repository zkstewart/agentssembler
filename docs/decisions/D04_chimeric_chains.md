# D04. Is a chain chimeric?

**Read:** the stage report, which crosses the two signals, then `10_groups/chain_breaks.tsv`. For candidates use `scripts/chain_profile.py` and `07_chains2/joins.verify.tsv`.

**A break needs at least two independent signals:**

1. Different Hi-C partners on each side (`partner_cos` well below 1), **and**
2. Weak reads across the join: few at MAPQ >= 20, with low local depth.

**Do not act on `cross_ratio` alone.** It is 0 at many satellite positions, whose flanks have no mappable Hi-C. It is also 0 where a tail of unique-looking sequence is near-identical to a sibling, so that Hi-C reads there get MAPQ 0. Read the `trans` counts on both sides first.

**Examples:**

- **Broken.** 13Q010 chain0011 at 15.78 Mb (chain0010 in the first builds): `partner_cos` 0.48 with 8,376 and 9,524 `trans` pairs, no MAPQ >= 20 reads (8 at any MAPQ), flank depth 10-12x against 27x, and the two sides contact different chromosomes. The 1.65 Mb chain0011b was a chr9 piece on a chr6 chain.
- **Not broken.** 13Q010 chain0055 has `cross_ratio` 0 at 14 joins and chain0117 at 7, but they are satellite (repeat fraction 0.85 and 0.72). chain0117 even passed both signals at 1.08 Mb (`partner_cos` 0.796, 3 of 11 reads at MAPQ >= 20): with 26-324 `trans` per side and `partner_cos` 0.72-0.85 at every join, its Hi-C is thin along the whole chain, with no change of partners at one point. On WatsonsChoice six chains showed `cross_ratio` 0 with a 100-170 kb tail, and no break was justified.

- **Not broken.** In the first Busk build (the 2026 rebuild had no candidates): chain0037 at 5.06 and 5.16 Mb: `partner_cos` 0.43-0.48 with plenty of `trans`, but 46-54 reads at MAPQ >= 20. Reads overrule Hi-C. Busk chain0021 at 5.91 Mb had 44 reads, 2 at MAPQ >= 20, `partner_cos` 0.71, and no clear change of partner, so no break either.

**A chain can be right and still need breaking.** The two-signal test is for chimeras. A different fault is a chain running through a fork, because a dropped unitig left the link looking one-to-one. It cannot be seen until groups exist, so `07_group` runs `scripts/review/chain_forks.py` for D06. Busk chain0061 had 705 kb at 1x then 839 kb at 2x, with Hi-C tied between homologs and the other exit leading to a chain on the homolog. Its reads are fine, since one haplotype does run that way. Break between arm and hub, then make the hub shared under D14. On the 2026 Busk rebuild `10_groups/hub_cores.d04.tsv` showed this hub at D04 already, so the break went in before grouping; the first build found it only at D06. Do not act on a fork that rejoins its own group: that is an ordinary bubble. The scan's `break_at` column gives the position. The same fault with chains in no group shows only after placement: stage 10 repeats the scan as `12_shared/chain_forks.tsv` for D16's tangle check. There, GammaCitrus's chain0073 ran through hub utg000537l, and the other exit was the route to its unplaced 779 kb allele ([TANGLES](../TANGLES.md)).

**Edit:** `edits/chain_breaks.txt`, one per line: `chain position reason`, split by whitespace, for example `chain0056 2852659 D06 HUB: utg000239l 2x, other exit to the homolog`. The position is the 1-based last base of the left part. The code reads the first two fields only, so the reason is for the record, but give it on every row. The left part keeps its name; the right becomes `<chain>b`. A chain can take several rows, each position on the unbroken chain: the parts are `<chain>`, `<chain>b`, `<chain>c`, ... from left to right, so a hub can be cut out from both sides. An empty file is a valid outcome, and it was the outcome on WatsonsChoice.

**After:** `run` applies the breaks without remapping Hi-C.

## Also at D04: collapsed hubs

`10_groups/hub_cores.d04.tsv` (collapsed samples) lists chains with a run of 2x unique-read depth inside 1x arms: hifiasm collapsed both haplotypes there. Break the `breakable` ones now, at `break_left` and `break_right`, together with any chimera breaks, after checking each edge against `07_chains2/depth.windows.tsv`. Give "D04 hub" as the reason. At D06 the same breaks cost a `reset 06_evidence --remap`. On the 292 rebuild, the hubs surfaced over two such rounds, while this table found them all before D04. Mark each freed core shared at D14. A core is not breakable when it lies inside one unitig, or when the piece the breaks would free carries long 1x arms (`core_2x_frac` under 0.75); stage 06 lists those again for D06 and D14.

## Also at D04: the Hi-C allele rescue

**Diploids only.** Above ploidy 2 leave it off, whatever the signal: since 0.3.8 `decide --set` and `run` refuse it there, and the reports no longer suggest it. On 13Q010 (triploid) its chromosome calls held, but 10 of 19 haplotype calls were wrong and the build gained no gene ([RECORD](../RECORD.md#13q010-rebuilt-with-037-2026-11)). Two of three haplotypes are often identical over a stretch, so exclusion leaves two groups and arbitrary balance picks one. A chain that continues into such a stretch reads as its allele. A bubble arm links two chains of one group and goes there by `graph`.

Read `09_hic_map/rescue_signal.txt`. Its Hi-C-blind chains have an allele partner, almost no unique Hi-C, and many pairs with an allele-ambiguous end: grouping will leave them unassigned however good the rest is. When they are a large share (the report's alarm: 1% of chain sequence or more), or carry many genes, record D04 with `--set hic_allele_rescue=true` and say why. This is the cheapest point: the setting acts in stages 06-08.

- Repeat chains are Hi-C-blind too (chain0091 on 292: 74% repeat). The rescue leaves chains over half repeat alone, so they do not count for it.
- The rescue knows which chromosome, rarely which haplotype: most of its haplotype calls are `arbitrary` or `allele_of_arbitrary`. Where Hi-C does not decide and the chain has a graph link to one haplotype's chains only, that link decides (`graph`): a tie-break, never set against a Hi-C call, since near-identical alleles can both link. It also brings orientation risk in the rescued region: on 292 inverted homolog blocks rose from 9 to 13, mostly in rescued chr2 sequence ([RECORD](../RECORD.md)). Say so at D13 and hand those regions over for Juicebox.
- **Identical end overlaps are graph links, not alleles.** `allele_table.py` sets apart blocks at de 0.002 or less that sit at an end of both chains, end to start, covering under 80% of either (`10_groups/alleles.dovetails.tsv`). Read as alleles, they had pushed a chain into the haplotype opposite the one it continues: on Busk, six rescue calls were moved by hand at D06 for this before the fix. Dovetails vote in the `graph` call and still count towards `redundant`. A dovetail into an end that another chain of that group already continues votes for nothing. Two calls come from them:
  - `bubble_arm`: such a dovetail, where the chain is also identical to part of the chain already there. It is that chain's other-haplotype arm, rejoining the same path, so it goes to the other group (Busk chain0166).
  - `fork` (left unassigned): a chain that links to every group of its chromosome, so it sits where the haplotypes meet. Tested again after the last round, since chains go longest first and a shorter one placed later can be the second link (Busk chain0210, whose `graph` call had cost a chr4_h2 bridge).
- It leaves these kinds of chain unassigned, each named in `10_groups/groups.rescued.tsv`:
  - `allele_in_every_group`: an allele of 50 kb or more in every group of the chromosome.
  - `redundant`: allele blocks of any length cover half the chain or more in every group. 292's chain0288 was a 105 kb near-copy whose copies were each under 50 kb. Placed, it fell between a chain and its telomeric tip.
  - `fork`: see above.
  - `copy`: identical (de 0.002 or less) over 80% or more of its length to one placed chain, and the two read under 1.25x haplotype depth together over that stretch (`11_ont_chains/depth.windows.tsv`): one haplotype's sequence twice. Summed near 2x it is a homozygous stretch and the allele call stands. On Busk three de-0 copies summed 1.55-1.85x (alleles, rightly) and one 0.79x (chain0181, a copy).
  - `collapsed`: a breakable 2x unique-read core (`10_groups/hub_cores.tsv`), so the chain holds both haplotypes there. That is for D06 to break; no haplotype call fits. A core inside one unitig cannot be broken off, so its chain is placed as usual.
  Collapsed and redundant chains still count towards their chromosome, so chains rescued in later rounds can anchor on them.
- Turned on later, it costs `reset 07_group` (from D05) or more.
