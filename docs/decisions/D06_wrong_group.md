# D06. Is a chain in the wrong group?

**Read:** the sizes and alarms in `reports/07_group.md`, `10_groups/eval_groups.txt` (an alarm, never a reason), `scripts/review/group_profile.py`, `scripts/review/weak_placed.py`.

**Signals:**

- A large chain whose group disagrees with its allele conflicts or with Hi-C.
- A group far larger than its siblings. On 13Q010, satellite chains chain0056 and chain0091 stacked on one chr8 group (48.2 Mb against 39.1 and 34.1).
- A high-copy chain passing because its density to one group is still highest. WatsonsChoice chain0184 (262x) had 26,032 pairs inside g03 and 97,601 outside it. **Read pairs outside the group as well as the best group, and check class `high` chains before ordering.**

**Subtelomeric satellite blocks.** Busk has seven chains, 0.8-8 Mb, with repeat fraction 0.57-0.77 and most of them with a telomere at one end; they align to several chromosomes at 3-7% divergence, so their allele conflicts are ignored and their reference hits mean little. Hi-C reaches them only thinly, and a group can look right by density while holding a false member. Work through the ends:

1. List telomeric ends per group, and unassigned chains with a telomere (the report does this).
2. For each such chain compare Hi-C per Mb of group, not raw counts, across all groups. Busk ratios ran from 2.5 (clear) to 1.2 (thin).
3. A group has only two ends. Use a free end to choose between candidates and the graph where it exists: chain0082 had a link to chain0013's satellite end. This is the only group decision where counting ends is the deciding evidence, so say so in the reason, and call it low confidence if Hi-C alone would not decide.

**Read the `HUB` rows of `10_groups/chain_forks.txt`** (from `scripts/review/chain_forks.py`; collapsed samples only). Each is a 2x unitig inside a chain, whose other exit reaches a chain on the homolog. `break_at` is the break that separates hub from arm. Rows to chains in no group are not flagged here, since most are ordinary bubbles until placement is known; D16 checks them. See D04: a break needs `reset 06_evidence --remap`, which costs a read mapping, so do this before spending time on order. The `hub_total` and `hub_unique` columns are the chain's own depth over the hub. A row reading about 1x there is an allele bubble, whose unitig class reads 2x only because its two alleles are near-identical. Leave it alone (292: most of 30 rows with the rescue on).

**Read `10_groups/hub_cores.tsv`** (stage 06, `scripts/review/hub_cores.py`). It lists chains classed single-copy that hold a run of 50 kb or more at 2x unique-read depth. The arms around such a core read 1x on MAPQ 0 reads, which keeps the chain's median under the shared class. This test needs no groups. HUB rows need a homolog's chain at the other exit, so with the Hi-C allele rescue on, they surface one round of breaks at a time: 292 took two remaps and still missed seven. A core is `breakable` when each inner edge lies at a unitig boundary, since `break_chain.py` splits only between unitigs. Break each breakable chain at `break_left` and `break_right`, unless it is already shared, and check the edges against the windows file first. Then mark the core shared (D14), and run `reset 06_evidence --remap` once for all of them. The rescue leaves breakable chains unassigned with the call `collapsed`. `inside` means hifiasm collapsed the hub within one unitig, so no D04 break can free it. On 292 all 15 cores were of this kind. Such a chain stays whole, and D14 chooses between two costs: on one haplotype, its core is missing from the other; shared, its 1x arms are duplicated. Weigh the genes on each part.

**Decide from** Hi-C, allele conflicts and graph links. The reference may only prompt the look.

**When evidence is thin, unassign.** Unplaced is a legitimate outcome. Reasons used:

- `Depth 262x; Hi-C spread over all groups` (high copy)
- `No graph link; 102 pairs to g01, none to its order neighbours` (no evidence for a position)

**Edit:** `edits/group_edits.tsv`, columns `chain, action, value, reason`. Actions are `group`, `unassigned .`, `shared`, `shared_in` (above ploidy 2, or across a translocation, [D14](D14_shared_chains.md)), `pair`. For `group` and `shared`, write the value as `with:chain0013`, meaning that chain's group or chromosome. A bare `g06` or `c4` also works, but labels are numbered by size and shift when anything upstream changes: on Busk a chain break swapped g06 with g07 and c4 with c5, and old labels would have put a chr2 hub on chr3.

**When edits move a whole haplotype between groups** (D05: two chromosomes merged and another split), add a `pair` row as well: `chainX  pair  with:chainY  reason` makes the group holding chainX, after the other edits, a haplotype of chainY's chromosome. Without it the group-to-chromosome pairing is copied unchanged from clustering, and ordering, shared copies and names would follow the old pairing. Here `with:` reads the grouping after the edits. GammaCitrus needed it: chr1 and chr7 haplotypes merged in one group, and one chr9 haplotype split in two.

**After:** `run` continues, as `08_bridges` applies the edits. If you find this later, edit and `reset 08_bridges`.

**A low-confidence chain: walk the graph from its gap end first.** Busk chain0042 (4.6 Mb, repeat fraction 0.77) had 1,144 links to its scaffold against 874 to another, and Hi-C was flat at its end, as at a control gap. The graph settled it: its end unitig and the next chain's arm converge on one unitig. Do not leave such a note open until validation, as I did.
