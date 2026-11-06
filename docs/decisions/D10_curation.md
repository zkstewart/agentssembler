# D10. Should a junction be set by hand?

**Almost always no.** Record `--reason "No junction needs a hand-set link."` when that is so.

**Curate only when** Hi-C is near a tie **and** a higher-ranked or homolog signal decides it. Never because the reference disagrees.

**Example:** 13Q010 chr8_h2. The scaffold audit showed a large inversion. Hi-C barely decided between layouts (454 of about 166,000 in log-likelihood, with the first 1.5 Mb of chain0005 unmappable). A link `chain0129+`/`chain0081+` from the alleles' homolog positions gave a scaffold colinear with both homologs, with telomeres at both ends.

**Lesson from the 13Q010 example:** the fix was right, and my explanation was wrong. I recorded a real 5 Mb inversion, and it does not exist. Test a biological story separately from the fix ([D12](D12_sv_calls.md)).

**Example:** Busk chain0194 and chain0081, described under [D09](D09_order.md).

**Edit:** `edits/order_links.tsv`: `chain0129+<TAB>chain0081+<TAB>curated<TAB>reason`. Then `reset 09_order --keepDecision`. Use the same chain names as in `12_shared/groups.tsv`, including `.gNN` for shared copies.
