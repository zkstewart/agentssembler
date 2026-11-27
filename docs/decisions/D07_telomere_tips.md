# D07. Should a telomeric tip be attached?

**Read:** `12_shared/graph.tips.tsv`, `16_validate/telomeres.chains.tsv`, `scripts/review/tel_ends.py`, `scripts/review/gap_walks.py`.

**Attach when** a short telomeric chain has a graph link to an assigned chain's end, in a group lacking a telomere there, and no group ends up with more than 2 telomeric ends.

**Leave unplaced when** the link is weak or conflicting, or when competing tips remain unresolved.

**Ask the user when** every read at the junction is MAPQ 0, because the tip is near-identical to another placed copy. The placement then rests on the graph alone and may be mosaic ([D13](D13_user_decisions.md)). On 13Q010 the user chose, in two builds, to place one such path (chain0135 and chain0182 at chr5_h2, 1.29 Mb with a connector): it carries 19 genes otherwise short of three copies and gives the scaffold its second telomere.

**Examples:**

- **Accepted.** 13Q010 chain0221 to g16 via chain0139+.
- **Left unplaced.** 13Q010 chain0255, the competing tip on chain0139+: it lies inside chain0221 at 0.19% divergence, and their depths sum to one haplotype (9.4 and 16.7x), so it is one tip assembled twice. chain0266 links an end that already carries a telomere (surplus).
- **Competing tips with Hi-C.** WatsonsChoice chain0530 and chain0539 both link to chain0305 with a telomere at the far end. Each went to a haplotype on 7 pairs. Report how thin that is.
- **Surplus.** WatsonsChoice chain0563 would be a third telomeric end in g14, so it was unplaced.

**Edit:**

1. `edits/group_edits.tsv`: `chain0221<TAB>group<TAB>with:chain0139<TAB>reason`. `with:` names the group by a chain in it; a bare `g16` shifts when groups are renumbered (T19).
2. `edits/order_links.tsv`: `chain0139+<TAB>chain0221-<TAB>telomere_tip<TAB>reason`

Then `agentssembler reset WORK 08_bridges --keepDecision`, because group edits are applied in that stage.
