# What the two genomes taught

Results that shaped the process. They come from one triploid and one diploid plant genome, so treat them as evidence, not laws.

## Graph-first works

| | Graph-first 13Q010 | WatsonsChoice |
|---|---|---|
| Ploidy, chromosomes | 3, 9 | 2, 9 |
| Placed | 921.3 Mb | 681.7 Mb |
| Unplaced | 40.9 Mb | 22.0 Mb |
| Gaps | 63 | 71 |
| Contig N50 | 16.80 Mb | 14.60 Mb |
| Telomeric ends | not rerun on the final version | 30 of 36 |
| compleasm complete | not rerun on the final version | 99.43%, 99.64% |

The final Juicebox review had not been done on either at the time of writing.

## Read correction

Raw SUP reads against corrected reads, at reads of 5 kb or more:

| | WatsonsChoice | 13Q010 |
|---|---|---|
| Reads kept | 60% | 83% |
| Bases kept | 65% | 88.5% |
| q20 spanning reads kept, at gaps | 38% | 82% |
| q20 spanning reads kept, control sites | 92% | 96% |
| Bases kept, reads of 100 kb or more | 21% | 89% |
| q20 in both forms, placed differently | 215 of 1.04 M | 267 of 0.90 M |

- **Correction loses and shortens reads near gaps.** At the same quality, reads near a WatsonsChoice gap are kept about half as often (34% against 63% at mean Q20-25). Of 151 raw spanning reads at q20, 50 had no corrected version and 45 no longer spanned, 34 of them shortened by more than 10%.
- **It does not move reads between haplotypes** on any scale measured.
- **It is not the main cause of gaps.** Raw depth also falls, from 31.6x at 20-50 kb to 4.4x at the edge, and 48 of 83 WatsonsChoice gaps have no spanning raw read. Why reads stop there is unknown: GC and complexity are ordinary.
- **Neither read set is enough alone.** Raw alone closed 11 WatsonsChoice gaps against 6 corrected, with 3 in common. On 13Q010 raw alone closed 15 against 19, and raw reads were inconsistent at 5 gaps that corrected reads filled.

## Automated labels that were wrong

| Label | What happened |
|---|---|
| `reads_overlap` | Closed two tandem-array gaps on WatsonsChoice, which reads anchored in the array cannot support. Raw-read audit reopened them |
| `only` tail | chain0256 on 13Q010 was already on all three haplotypes |
| Read count at joins | The worst join on 13Q010 had 8 reads, all at MAPQ 0 |
| `make_names` | Needed a baseline; a stale flag reversed chr6_h2 |
| High group density | A 262x chain passed the group test on density alone |

## Errors of mine worth learning from

- I stated coverage from reasoning, not the histogram (9x per haplotype; it was 27x).
- I recorded a 5 Mb inversion on 13Q010 that does not exist, and another reached a curation file.
- I called reused arms legitimate before testing them.
- I said an ordering error explained an inversion. It survived correct ordering.
- A hard-coded `!= 3` in `dosage.py` listed the wrong windows on a diploid. The fractions it printed were right, but the list of runs was not. Now fixed.
- I gave a pooled corrected/raw depth ratio of 17-23% before finding two arrays had inflated it. The median is about 50%.
- Writing raw reads to an uncompressed FASTA filled the host disk and crashed the machine.

## Open questions

- **Why do reads stop at gaps?** See [IMPROVING](IMPROVING.md).
- **Is collapsed handling sound above ploidy 2?** It is generalised in code, but only ploidy 2 has been run.
- **Juicebox edits:** the route back into the AGP has not been exercised.
