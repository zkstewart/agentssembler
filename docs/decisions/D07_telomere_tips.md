# D07. Should a telomeric tip be attached?

**Read:** `12_shared/graph.tips.tsv`, `16_validate/telomeres.chains.tsv`, `scripts/review/tel_ends.py`, `scripts/review/gap_walks.py`.

**Attach when** a short telomeric chain has a graph link to an assigned chain's end, in a group lacking a telomere there, and no group ends up with more than 2 telomeric ends.

**Leave unplaced when** the link is weak or conflicting, or when competing tips remain unresolved.

**Examples:**

- **Accepted.** 13Q010 chain0221 to g16 via chain0139+.
- **Left unplaced.** 13Q010 chain0182 and chain0265. chain0182 is a subtelomeric repeat shared among chromosome ends: raw reads later confirmed that its anchors are unreliable.
- **Competing tips with Hi-C.** WatsonsChoice chain0530 and chain0539 both link to chain0305 with a telomere at the far end. Each went to a haplotype on 7 pairs. Report how thin that is.
- **Surplus.** WatsonsChoice chain0563 would be a third telomeric end in g14, so it was unplaced.

**Edit:**

1. `edits/group_edits.tsv`: `chain0221<TAB>group<TAB>g16<TAB>reason`
2. `edits/order_links.tsv`: `chain0139+<TAB>chain0221-<TAB>telomere_tip<TAB>reason`

Then `agentssembler reset WORK 08_bridges --keepDecision`, because group edits are applied in that stage.
