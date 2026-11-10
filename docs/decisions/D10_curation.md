# D10. Should a junction be set by hand?

**Almost always no.** Record `--reason "No junction needs a hand-set link."` when that is so.

**Curate only when** Hi-C is near a tie **and** a higher-ranked or homolog signal decides it. Never because the reference disagrees.

**Example:** a scaffold whose audit showed a large inversion. Hi-C barely decided between layouts (0.3% of the log-likelihood, with the first 1.5 Mb of one chain unmappable). A link between two chains, taken from their alleles' positions on the homologs, gave a scaffold colinear with both homologs, with telomeres at both ends.

**Lesson from that example:** the fix was right, and my explanation was wrong. I recorded a real inversion of several megabases, and it does not exist. Test a biological story separately from the fix ([D12](D12_sv_calls.md)).

**Example:** the two chains described at the end of [D09](D09_order.md).

**Edit:** `edits/order_links.tsv`: `<chainA>+<TAB><chainB>+<TAB>curated<TAB>reason`. Then `reset 09_order --keepDecision`. Use the same chain names as in `12_shared/groups.tsv`, including `.gNN` for shared copies.
