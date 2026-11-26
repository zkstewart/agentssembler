# D10. Should a junction be set by hand?

**Almost always no.** Record `--reason "No junction needs a hand-set link."` when that is so.

**Curate only when** Hi-C is near a tie **and** a higher-ranked or homolog signal decides it. Never because the reference disagrees.

**Example:** 13Q010 chr8_h2. The scaffold audit showed a large inversion. Hi-C barely decided between layouts (454 of about 166,000 in log-likelihood, with the first 1.5 Mb of chain0005 unmappable). A link `chain0129+`/`chain0081+` from the alleles' homolog positions gave a scaffold colinear with both homologs, with telomeres at both ends.

**Lesson from the 13Q010 example:** the fix was right, and my explanation was wrong. I recorded a real 5 Mb inversion, and it does not exist. Test a biological story separately from the fix ([D12](D12_sv_calls.md)).

**Example:** Busk chain0194 and chain0081, described under [D09](D09_order.md).

**Edit:** `edits/order_links.tsv`: `chain0129+<TAB>chain0081+<TAB>curated<TAB>reason`. Then `reset 09_order --keepDecision`. Use the same chain names as in `12_shared/groups.tsv`, including `.gNN` for shared copies.

**A link sets order and orientation only.** The builder merges sequence only where the graph has an overlap of at least 5 kb between the two ends; anywhere else the junction is a Hi-C gap. If the two ends overlap outside the graph (an identical dovetail in `10_groups/alleles.dovetails.tsv`), a `curated` link leaves that overlap on both sides of the gap, as a duplicate. Give the type as `overlap:N` instead, with N the overlap in bp, measured end to end on the chain sequences (minimap2 `-x asm5` of the two ends: one block, no gaps). The build checks the N bp at both ends and joins them with no gap at 99.8% identity or more; otherwise `13_scaffold/build.log` says so and the junction stays a gap. Busk: `chain0148+<TAB>chain0131+<TAB>overlap:40154`, a 40 kb identical overlap of two chains on a homozygous stretch whose only graph link runs through the other haplotype's copy. As `curated` it left the 40 kb on both sides of a gap. The same `reset 09_order --keepDecision` applies.
