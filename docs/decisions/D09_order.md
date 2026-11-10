# D09. Is an odd order real or a search artefact?

**Read:** `12_shared/unhonoured.txt` (no reference needed), `12_shared/order.tsv` (`position_llr`, `orientation_llr`), `12_shared/eval_order.txt` (alarm), and `scripts/review/order_review.py`.

**Signals:** a same-group graph link not honoured, an unused bridge, a telomere facing inward, low collinearity in one group, or a small `position_llr`.

**How to decide:**

1. **Is the search stuck?** Multi-start search is the default (40 starts). If it does not improve, the objective itself prefers this order.
2. **Test explicit alternatives** with `scripts/score_orders.py`.
3. **Are weak chains driving it?** Hundreds of satellite pairs can outweigh megabases of good chains. Two-stage ordering (`--stats`) fixed this once. Check `min_within`.
4. **Was a path missed?** `hops: 12` is the default. 3 hops missed a 4 Mb hub path on one genome and reversed a block; 8 missed an allele path on another once connectors were in the graph (see `SETTINGS.md`).
5. **Check the orientation of hub chains with gene evidence.** Run compleasm on the unplaced sequence later (`15_validate`): genes missing from both haplotypes but in unplaced sequence mean something is wrong with order or placement. One ordering error was found that way.

**Edit:**

- `edits/order_edits.tsv` to unplace or add a chain at a group's end, with a reason. Then `reset 09_order --keepDecision`.
- A changed link needs [D10](D10_curation.md).

**A link can be real and still outside the search.** `make_order_links.py` rewards only overlaps of 5 kb or more, or 3 or more reads, and `unhonoured.py` lists all same-group graph links. One chain of about 150 kb had overlaps of 2.6 and 3.5 kb to two chains, consistent in orientation, but was placed by a noisy Hi-C score 9 ranks away, at a cost of 0.05% of the score. Its contacts were densest with those neighbours (`review/weak_placed.py`, `top_density_ranks`) and it was allelic to the chain at that position on the homolog. Check each line of `unhonoured.txt` this way; curate under D10 when the signals agree. Check alleles on the homolog too: another chain's orientation was a tie in Hi-C, and fixed by its alleles' order.

**Do not reorder toward the reference.** A disagreement with collinear support from reads or Hi-C is an SV candidate.
