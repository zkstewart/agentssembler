# D05. Why not the expected number of groups?

**Read:** `10_groups/hic_group.log`, the merge scores it prints, and `reports/07_group.md`.

**Do not force the count.** Look at what blocks the merge.

**Causes seen, and their fixes:**

| Signal | Cause | Fix |
|---|---|---|
| One group too many, or a chromosome with too few haplotypes | Satellite chains' alignments look like allelic conflicts and block merges | Built in: held out of seeding (`--stats`), conflicts ignored above repeat fraction 0.5 (`--repeat`). Check `max_repeat` and `min_within` |
| Two chromosomes merged and another split | A point of thousands of pairs between two chains: organelle reads on nuclear insertions, with the organelle missing from the target | Stage 05 maps with decoys to prevent this ([LESSONS](../LESSONS.md)). Check `09_hic_map/decoys.txt` holds the organelle, and look for such a point in `pairs.tsv.gz` |
| One chromosome split by haplotype and two others merged | Divergent homologs: Hi-C between them is as weak as between chromosomes | `hap_first: true` clusters directly to chromosomes x ploidy, then gathers the groups into chromosomes, `ploidy` at a time, by allele alignment weighted by identity. It is the default at ploidy 2. **Above ploidy 2 it is off by default and has not been run on real data:** set it in `settings` when chromosome-first grouping shows this signal |

**Independent support to report:**

- With `hap_first`, a cliff in merge scores right after the target (WatsonsChoice: 114.8 to 38.9). Groups pair clearly by allelic bp (7.4-30.8 Mb against at most 0.1 Mb for the next group).
- The cliff can be absent. On Busk (`...20:15.5 19:8.5 18:7.1 17:6.8`) the last merges were no better than the next, while pairing still left no doubt (16-30 Mb allelic against at most 0.3 Mb for the next group). A missing cliff is a reason to look, and the allele pairing is the stronger test.
- Bridges and links falling within one group (13Q010: 33/33 and 38/42).
- No group with more than 2 telomeric ends.
- With `hap_first` above ploidy 2, the log prints `set` lines in place of `pair` lines: the allelic bp of the weakest pair inside each set, and the best to any group outside it. A `WARNING` on a set means it holds fewer haplotypes than `ploidy`.

**An allopolyploid entered as twice the chromosomes at ploidy 2** ([SETTINGS](../SETTINGS.md)) has homeologs that also align. Pairing counts each alignment's bp times (1 - divergence)^10, so a homolog at 1% outranks a homeolog at 8% of the same length by about two to one. Check the `next best for either` figure on each `pair` line: if it is close to the pair's own, the subgenomes are too similar for this to separate, and the pairs need another test before ordering.

**The expected number itself may be wrong.** Chromosome number from a reference is an assumption. If the data keep disagreeing, tell the user.

**Changing settings:** edit `settings` in `WORK/config.yaml`, then `reset 07_group`.
