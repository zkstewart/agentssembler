# D09. Is an odd order real or a search artefact?

**Read:** `12_shared/unhonoured.txt` (no reference needed), `12_shared/order.tsv` (`position_llr`, `orientation_llr`), `12_shared/eval_order.txt` (alarm), and `scripts/review/order_review.py`.

**Signals:** a same-group graph link not honoured, an unused bridge, a telomere facing inward, low collinearity in one group, or a small `position_llr`.

**How to decide:**

1. **Is the search stuck?** Multi-start search is the default (40 starts). If it does not improve, the objective itself prefers this order.
2. **Test explicit alternatives** with `scripts/score_orders.py`.
3. **Are weak chains driving it?** Hundreds of satellite pairs can outweigh megabases of good chains. Two-stage ordering (`--stats`) fixed this once. Check `min_within`.
4. **Was a path missed?** `hops: 8` is the default because 3 hops missed a 4 Mb hub path on WatsonsChoice and reversed a block. On 13Q010 it made no difference.
5. **Check the orientation of hub chains with gene evidence.** Run compleasm on the unplaced sequence later (`15_validate`): genes missing from both haplotypes but in unplaced sequence mean something is wrong with order or placement. That is how the WatsonsChoice chr3 error was found.

**Tried and dropped:** a reference synteny term (`synteny_links.py`) changed nothing over four rounds, and using it would break the evidence rules.

**Edit:**

- `edits/order_edits.tsv` to unplace or add a chain at a group's end, with a reason. Then `reset 09_order --keepDecision`.
- A changed link needs [D10](D10_curation.md).

**Do not reorder toward the reference.** A disagreement with collinear support from reads or Hi-C is an SV candidate.
