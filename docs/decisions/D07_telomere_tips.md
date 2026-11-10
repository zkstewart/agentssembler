# D07. Should a telomeric tip be attached?

**Read:** `12_shared/graph.tips.tsv`, `16_validate/telomeres.chains.tsv`, `scripts/review/tel_ends.py`, `scripts/review/gap_walks.py`.

**Attach when** a short telomeric chain has a graph link to an assigned chain's end, in a group lacking a telomere there, and no group ends up with more than 2 telomeric ends.

**Leave unplaced when** the link is weak or conflicting, or when competing tips remain unresolved.

**Ask the user when** every read at the junction is MAPQ 0, because the tip is near-identical to another placed copy. The placement then rests on the graph alone and may be mosaic ([D13](D13_user_decisions.md)). Tell the user what the placement would add: genes completed, and whether the scaffold gains a telomere.

**Examples:**

- **Accepted.** A tip with one graph link, to the end of an assigned chain in a group with no telomere at that end.
- **Left unplaced.** A tip that is a subtelomeric repeat shared among chromosome ends: raw reads later confirmed that its anchors are unreliable.
- **Competing tips with Hi-C.** Two tips linked to the same chain, which had a telomere at its far end. Each went to a haplotype on 7 pairs. Report how thin that is.
- **Surplus.** A tip that would have been a third telomeric end in its group was unplaced.

**Edit:**

1. `edits/group_edits.tsv`: `<tip><TAB>group<TAB>with:<the chain it links to><TAB>reason`
2. `edits/order_links.tsv`: `<chain>+<TAB><tip>-<TAB>telomere_tip<TAB>reason`, with the orientations from `graph.tips.tsv`

Then `agentssembler reset WORK 08_bridges --keepDecision`, because group edits are applied in that stage.
