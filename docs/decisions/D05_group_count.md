# D05. Why not the expected number of groups?

**Read:** `10_groups/hic_group.log`, the merge scores it prints, and `reports/07_group.md`.

**Do not force the count.** Look at what blocks the merge.

**Causes seen, and their fixes:**

| Signal | Cause | Fix |
|---|---|---|
| One group too many, or a chromosome with too few haplotypes | Satellite chains' alignments look like allelic conflicts and block merges | Built in: held out of seeding (`--stats`), conflicts ignored above repeat fraction 0.5 (`--repeat`). Check `max_repeat` and `min_within` |
| Two chromosomes merged and another split | A point of thousands of pairs between two chains: organelle reads on nuclear insertions, with the organelle missing from the target | Stage 05 maps with decoys to prevent this ([LESSONS](../LESSONS.md)). Check `09_hic_map/decoys.txt` holds the organelle, and look for such a point in `pairs.tsv.gz` |
| One chromosome split by haplotype and two others merged | Divergent homologs: Hi-C between them is as weak as between chromosomes | `hap_first: true` clusters directly to chromosomes x ploidy. **It pairs groups into chromosomes by allele alignment, for ploidy 2 only** |

**Independent support to report:**

- With `hap_first`, a cliff in merge scores right after the target (a threefold drop on one genome). Groups pair clearly by allelic bp (megabases to tens of megabases, against about 0.1 Mb for the next group).
- The cliff can be absent. On another genome the last merges scored no better than the next would have, while pairing still left no doubt (tens of megabases allelic against well under 1 Mb for the next group). A missing cliff is a reason to look, and the allele pairing is the stronger test.
- Bridges and links falling within one group (all bridges and about 90% of links, on one genome).
- No group with more than 2 telomeric ends.

**The expected number itself may be wrong.** Chromosome number from a reference is an assumption. If the data keep disagreeing, tell the user.

**Changing settings:** edit `settings` in `WORK/config.yaml`, then `reset 07_group`.
