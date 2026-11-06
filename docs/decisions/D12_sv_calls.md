# D12. Is an SV candidate real?

**A claim about biology must survive a direct test before it is written down.** I made two false inversion claims on 13Q010, and one reached a curation file.

**Candidates:** `16_validate/sv_candidates.tsv` (reference), `16_validate/sv/inv/merged.tsv` (homologs). They are a list to test, not findings.

**Apply in this order, and record each result in `16_validate/sv/sv_calls.tsv`:**

1. **Whole-scaffold direction first.** A scaffold mostly inverted against its homolog is a writing-direction flag, not an inversion. The report raises this.
2. **Dosage:** `16_validate/dosage.txt`. Do translocation partners balance to the ploidy?
3. **Homologs:** align the region to all scaffolds **except its own** (`scripts/validate/summarise_regions.py`). A self-hit otherwise hides everything else. A real inversion is `+ / - / +` against both homologs. Colinear with even one means it is not specific to this haplotype.
4. **Reads at the exact transition:** `sv/bp_reads.tsv`, or `scripts/validate/transitions.py`. Take the boundary from step 3, not the reference block's edge.
5. **Hi-C flanks:** `scripts/sv_flanks.py`, at the breakpoints and every 2 Mb for background. A ratio `cross / sqrt(L_self x R_self)` far below background, with well-covered flanks, is a warning. **A flank with under 500 self pairs says nothing**: it is satellite, not a chimera.
6. **Orientation:** `sv/inversion_test.txt`. Near-end contacts should beat far-end. Hi-C is uninformative for blocks under about 1 Mb at a 500 kb window. Where a block's end has low-Hi-C arms, near/far misleads, so read `scripts/validate/block_profile.py`.

**Grades:**

- **High:** breakpoint inside a contig or unitig, reads span it, Hi-C normal, homologs agree.
- **Medium:** one element weak or untestable.
- **Low:** few signals, usually a small block with flat Hi-C.
- **Not supported:** no coherent homolog evidence, or satellite. Withdraw claims openly.

A breakpoint inside one unitig has read-level support from the assembler. At a Hi-C-only gap it has much less. Say which it is.

**An opposite direction between homologs with flat Hi-C and no reads is undetermined, not an inversion.** WatsonsChoice chr5 chains 0275 and 0350 are the example.

**Never fix toward the reference.** Report the disagreement.
