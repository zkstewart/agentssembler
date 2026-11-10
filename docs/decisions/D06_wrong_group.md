# D06. Is a chain in the wrong group?

**Read:** the sizes and alarms in `reports/07_group.md`, `10_groups/eval_groups.txt` (an alarm, never a reason), `scripts/review/group_profile.py`, `scripts/review/weak_placed.py`.

**Signals:**

- A large chain whose group disagrees with its allele conflicts or with Hi-C.
- A group far larger than its siblings. Satellite chains once stacked on one group, leaving it a quarter larger than its largest sibling.
- A high-copy chain passing because its density to one group is still highest. One chain at about nine times haplotype depth had nearly four times as many pairs outside its group as inside it. **Read pairs outside the group as well as the best group, and check class `high` chains before ordering.**

**Subtelomeric satellite blocks.** A genome can have several chains of megabase size with a repeat fraction over 0.5, most of them with a telomere at one end. They align to several chromosomes at 3-7% divergence, so their allele conflicts are ignored and their reference hits mean little. Hi-C reaches them only thinly, and a group can look right by density while holding a false member. Work through the ends:

1. List telomeric ends per group, and unassigned chains with a telomere (the report does this).
2. For each such chain compare Hi-C per Mb of group, not raw counts, across all groups. Ratios of the best group to the next have run from 2.5 (clear) to 1.2 (thin).
3. A group has only two ends. Use a free end to choose between candidates and the graph where it exists, such as a link from the chain to the satellite end of a chain already in the group. This is the only group decision where counting ends is the deciding evidence, so say so in the reason, and call it low confidence if Hi-C alone would not decide.

**Read the `HUB` rows of `10_groups/chain_forks.txt`** (from `scripts/review/chain_forks.py`; collapsed samples only). Each is a 2x unitig inside a chain, whose other exit reaches a chain on the homolog. See D04: a break needs `reset 06_evidence --remap`, which costs a read mapping, so do this before spending time on order. Also look at embedded 2x runs in `11_ont_chains/depth.windows.tsv` for assigned chains; one with Hi-C tied between homologs is the same fault.

**Decide from** Hi-C, allele conflicts and graph links. The reference may only prompt the look.

**When evidence is thin, unassign.** Unplaced is a legitimate outcome. Reasons used:

- `Depth 9x haplotype; Hi-C spread over all groups` (high copy)
- `No graph link; about 100 pairs to its group, none to its order neighbours` (no evidence for a position)

**Edit:** `edits/group_edits.tsv`, columns `chain, action, value, reason`. Actions are `group`, `unassigned .`, `shared`. For `group` and `shared`, write the value as `with:<chain>`, meaning that chain's group or chromosome. A bare `g06` or `c4` also works, but labels are numbered by size and shift when anything upstream changes: a chain break once swapped two group labels and two chromosome labels, and the old labels would have put a hub on the wrong chromosome.

**After:** `run` continues, as `08_bridges` applies the edits. If you find this later, edit and `reset 08_bridges`.

**A low-confidence chain: walk the graph from its gap end first.** One chain of several megabases, with a repeat fraction over 0.7, had only 1.3 times as many Hi-C links to its scaffold as to another, and Hi-C was flat at its end, as at a control gap. The graph settled it: its end unitig and the next chain's arm converge on one unitig. Do not leave such a note open until validation, as I did.
