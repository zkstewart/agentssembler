# D09. Is an odd order real or a search artefact?

**Read:** `12_shared/unhonoured.txt` (no reference needed), `12_shared/order.tsv` (`position_llr`, `orientation_llr`), `12_shared/eval_order.txt` (alarm), and `scripts/review/order_review.py`.

**Signals:** a same-group graph link not honoured, an unused bridge, a telomere facing inward, low collinearity in one group, or a small `position_llr`.

**How to decide:**

1. **Is the search stuck?** After the first search, random restarts run until 3 in a row fail to beat the best (at most 40). Each group's log line (`12_shared/by_chrom/cN.log`) gives restarts run, those that improved, and those that reached the best order. If none improved, the objective itself prefers this order. If some improved, the report raises an alarm: the search has local optima there, so read that group's LLRs before trusting its junctions. On 292 the patience rule cut chr2's ordering (52 and 53 chains per group) from about 2 h to 6 min, and every order in four builds (36 chromosomes) stayed the same. A group no restart reaches is not wrong for that: on 292 these were chr8_h1 (7 chains, the t(4;8)/t(5;8) derivative) and chr3_h2 (11), where the greedy start is hard to find from a random one.
2. **Test explicit alternatives** with `scripts/score_orders.py`.
3. **Are weak chains driving it?** Hundreds of satellite pairs can outweigh megabases of good chains. Two-stage ordering (`--stats`) fixed this once. Check `min_within` (`reports/06_evidence.md`), and pass it as the fifth argument to `score_orders.py`: without it that script uses 700.
4. **Was a path missed?** `hops: 12` is the default. 3 hops missed a 4 Mb hub path on WatsonsChoice and reversed a block; 8 missed a Busk allele path once connectors were in the graph (see `SETTINGS.md`).
5. **Check the orientation of hub chains with gene evidence.** Run compleasm on the unplaced sequence later (`15_validate`): genes missing from both haplotypes but in unplaced sequence mean something is wrong with order or placement. That is how the WatsonsChoice chr3 error was found.

**Edit:**

- `edits/order_edits.tsv` to unplace or add a chain at a group's end, with a reason. Then `reset 09_order --keepDecision`.
- A changed link needs [D10](D10_curation.md).

**A link can be real and still outside the search.** `make_order_links.py` rewards only overlaps of 5 kb or more, or 3 or more reads, and `unhonoured.py` lists all same-group graph links. Busk chain0193 (146 kb, 434 pairs) had overlaps of 2.6 and 3.5 kb to chain0080 and chain0011, consistent in orientation, but was placed by a noisy Hi-C score 8 ranks away (`position_llr` 172). The same junction came back in each build. Its contacts were densest with those neighbours (`review/weak_placed.py`, `top_density_ranks`) and it was allelic to chain0034. Check each line of `unhonoured.txt` this way; curate under D10 when the signals agree. Check alleles on the homolog too: chain0080's orientation was weak in Hi-C (`orientation_llr` 29.5), and fixed by the same graph links. Both overlaps are under the builder's 5 kb cut, so give the links as `overlap:N` ([D10](D10_curation.md)) or they stay gaps with the overlap on both sides.

**Do not reorder toward the reference.** A disagreement with collinear support from reads or Hi-C is an SV candidate.
