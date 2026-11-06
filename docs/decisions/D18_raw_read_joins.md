# D18. Which read-derived joins survive raw reads?

**Only with raw BAMs.** Corrected reads remain first choice for fills, because single raw reads are too noisy to compare. Raw reads recover overlap closures that correction lost, and they audit joins.

**Read:** `17_raw_reads/audit.tsv` (one row per join with a suggestion), `cmp.by_length.tsv`, `cmp.by_q.tsv`, `cmp.sites.tsv`, and `reports/13_rawreads.md`.

**The suggestion is a prompt, not a verdict. Each flag means:**

| Flag | Test | What it catches |
|---|---|---|
| `under_2_clean_q20_reads` | Reads through at MAPQ >= 20 with no indel over 50 bp | Unsupported joins. Use judgement with shared flanks, where MAPQ is low for good reason |
| `array_depth` | Alignment depth within 20 kb is over 4 x haplotype | A tandem array. Anchors inside it pass MAPQ and agreement but say nothing about its length. **Reopen.** WatsonsChoice chr8_h1:15713522 and chr8_h2:14209070, both accepted by the corrected fill |
| `trim_has_no_copy_in_kept_flank` | The trimmed sequence has a copy of at least 99% identity over at least 90% of its length in the kept flank | A trim that discards real sequence. 13Q010 chr2_h3 (12.9 of 29.7 kb at 89.8%), WatsonsChoice chr8_h2:5743770 (96.6%). **Reject** |
| `length_ratio_off` | Read/target length within 2% of 1 | Collapse or a kept duplication |
| `reads_break_here` | Many reads stop aligning at the join | A bad join |
| `low_identity_at_join` | Identity at the join under 95% | Wrong sequence, or low-quality reads. Compare with `id_flank` |

**Also check:**

- A raw-read fill keeps about 1% error, with no polishing here. Say so for each.
- Low-quality raw reads can look like disagreement: at chr7_h1:11467902 two reads at 80% identity looked like an 83% conflict, and the fill was correct.
- Do not read pooled depth ratios: use medians per flank. Two arrays inflated the pooled figure from about 50% to 17-23% in my first report.

**Correction behaviour to expect:** it loses and shortens reads near gaps, and it does not move reads between haplotypes. See [LESSONS](../LESSONS.md).

**Edit:** `edits/join_decisions.tsv`: `scaffold, left_unitig, right_unitig, reason`, with unitig names without orientation, for every join to reverse or reject. It works on corrected joins as well as raw ones. Then `run`.
