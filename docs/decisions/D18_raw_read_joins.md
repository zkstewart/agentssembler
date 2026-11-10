# D18. Which read-derived joins survive raw reads?

**Only with raw BAMs.** Corrected reads remain first choice for fills, because single raw reads are too noisy to compare. Raw reads recover overlap closures that correction lost, and they audit joins.

**Read:** `17_raw_reads/audit.tsv` (one row per join with a suggestion), `cmp.by_length.tsv`, `cmp.by_q.tsv`, `cmp.sites.tsv`, and `reports/13_rawreads.md`.

**The suggestion is a prompt, not a verdict. Each flag means:**

| Flag | Test | What it catches |
|---|---|---|
| `under_2_clean_q20_reads` | Reads through at MAPQ >= 20 with no indel over 50 bp | Unsupported joins. Use judgement with shared flanks, where MAPQ is low for good reason |
| `array_depth` | Alignment depth within 20 kb is over 4 x haplotype | A tandem array. Anchors inside it pass MAPQ and agreement but say nothing about its length. **Reopen.** The corrected fill has accepted two such joins on one genome |
| `trim_has_no_copy_in_kept_flank` | The trimmed sequence has a copy of at least 99% identity over at least 90% of its length in the kept flank | A trim that discards real sequence. One such trim matched its flank at 96.6%. **Reject.** Check by aligning the trim to the flank's end yourself (`minimap2 -x asm20 -c`): a join was once rejected on "12.9 of 29.7 kb at 89.8%", when the whole trim aligns at 99.3% and 10 clean q20 reads cross the join |
| `length_ratio_off` | Read/target length within 2% of 1 | Collapse or a kept duplication |
| `reads_break_here` | Many reads stop aligning at the join | A bad join |
| `low_identity_at_join` | Identity at the join under 95% | Wrong sequence, or low-quality reads. Compare with `id_flank` |

**Also check:**

- A raw-read fill keeps about 1% error, with no polishing here. Say so for each.
- Low-quality raw reads can look like disagreement: at one join two reads at 80% identity looked like an 83% conflict, and the fill was correct.
- Do not read pooled depth ratios: use medians per flank. Two arrays inflated the pooled figure from about 50% to 17-23% in my first report. Pooled spanner counts mislead too: on one genome, 37 of 40 q20 raw spanners were at one array.
- Many spanners with disagreeing lengths mean an array, and they still confirm the adjacency. Say that, and leave the gap (one gap had 306 spanners with 8-14 kb of disagreement).
- `15_validate` writes `join_reads.tsv` for merged joins, which this audit does not cover. `no_q20` with many reads at low MAPQ is usually near-identical alleles: the read shows the join is possible, not which haplotype it belongs to.

**Correction behaviour to expect:** it loses and shortens reads near gaps, and it does not move reads between haplotypes. See [LESSONS](../LESSONS.md).

**Edit:** `edits/join_decisions.tsv`: `scaffold, left_unitig, right_unitig, reason`, with unitig names without orientation, for every join to reverse or reject. It works on corrected joins as well as raw ones. Then `run`.
