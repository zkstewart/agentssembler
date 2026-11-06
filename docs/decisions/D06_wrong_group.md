# D06. Is a chain in the wrong group?

**Read:** the sizes and alarms in `reports/07_group.md`, `10_groups/eval_groups.txt` (an alarm, never a reason), `scripts/review/group_profile.py`, `scripts/review/weak_placed.py`.

**Signals:**

- A large chain whose group disagrees with its allele conflicts or with Hi-C.
- A group far larger than its siblings. On 13Q010, satellite chains chain0056 and chain0091 stacked on one chr8 group (48.2 Mb against 39.1 and 34.1).
- A high-copy chain passing because its density to one group is still highest. WatsonsChoice chain0184 (262x) had 26,032 pairs inside g03 and 97,601 outside it. **Read pairs outside the group as well as the best group, and check class `high` chains before ordering.**

**Decide from** Hi-C, allele conflicts and graph links. The reference may only prompt the look.

**When evidence is thin, unassign.** Unplaced is a legitimate outcome. Reasons used:

- `Depth 262x; Hi-C spread over all groups` (high copy)
- `No graph link; 102 pairs to g01, none to its order neighbours` (no evidence for a position)

**Edit:** `edits/group_edits.tsv`, columns `chain, action, value, reason`. Actions are `group gNN`, `unassigned .`, `shared cN`.

**After:** `run` continues, as `08_bridges` applies the edits. If you find this later, edit and `reset 08_bridges`.
