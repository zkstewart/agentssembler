# D03. Which weak joins are real?

**Read:** `07_chains2/joins.verify.tsv`, `07_chains2/ont_vs_chains.paf.gz.log` and the stage report. Check first that the self-alignment finished (`06_chains/chains_self.paf.gz.log`): an incomplete `chains_self.paf.gz` undermines round 2.

**The rule, applied automatically:** fewer than 3 spanning reads breaks the join. It removed 1 of 1,035 joins on 13Q010.

**Decide:** look at joins under 10 reads, for their MAPQ >= 20 count and local depth.

**Limit:** this rule missed the worst join on 13Q010, which had 8 reads, all at MAPQ 0 (see [D04](D04_chimeric_chains.md)). Read count alone is not enough. You are not expected to find every chimera here: D04 tests them another way.

**Long overlaps are barely tested by reads.** A read must cover the join to count, but where two chains overlap by more than read length (median 66 kb on Busk, 28 of 289 over 100 kb), reads confirm the sequence inside the overlap, not that the chains connect. Look at dovetail joins over about 30 kb: do any reads run from unique sequence on one side into unique sequence on the other, and do flanks have sensible depth? `07_chains2/chains.contained.tsv` lists unitigs wholly inside such an overlap.

**Edit:** add lines to `edits/join_breaks.txt`, in the form `utg000364l+ utg001421l-<TAB>reason`, then `reset 04_chains`.
