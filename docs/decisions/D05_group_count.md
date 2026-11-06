# D05. Why not the expected number of groups?

**Read:** `10_groups/hic_group.log`, the merge scores it prints, and `reports/07_group.md`.

**Do not force the count.** Look at what blocks the merge.

**Causes seen, and their fixes:**

| Signal | Cause | Fix |
|---|---|---|
| One group too many, or a chromosome with too few haplotypes | Satellite chains' alignments look like allelic conflicts and block merges | Built in: held out of seeding (`--stats`), conflicts ignored above repeat fraction 0.5 (`--repeat`). Check `max_repeat` and `min_within` |
| One chromosome split by haplotype and two others merged | Divergent homologs: Hi-C between them is as weak as between chromosomes | `hap_first: true` clusters directly to chromosomes x ploidy. **It pairs groups into chromosomes by allele alignment, for ploidy 2 only** |

**Independent support to report:**

- With `hap_first`, a cliff in merge scores right after the target (WatsonsChoice: 114.8 to 38.9). Groups pair clearly by allelic bp (7.4-30.8 Mb against at most 0.1 Mb for the next group).
- Bridges and links falling within one group (13Q010: 33/33 and 38/42).
- No group with more than 2 telomeric ends.

**The expected number itself may be wrong.** Chromosome number from a reference is an assumption. If the data keep disagreeing, tell the user.

**Changing settings:** edit `settings` in `WORK/config.yaml`, then `reset 07_group`.
