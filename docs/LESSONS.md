# What the genomes taught

Results that shaped the process. They come from one triploid and two diploid citrus genomes, so treat them as evidence, not laws.

## Graph-first works

| | Graph-first 13Q010 | WatsonsChoice | Busk |
|---|---|---|---|
| Ploidy, chromosomes | 3, 9 | 2, 9 | 2, 9 |
| Haplotype depth | 27x | 28x | 62x |
| Placed | 921.3 Mb | 681.7 Mb | 615.3 Mb |
| Unplaced | 40.9 Mb | 22.0 Mb | 11.9 Mb |
| Gaps | 63 | 71 | 7 |
| Contig N50 | 16.80 Mb | 14.60 Mb | 34.47 Mb |
| Telomeric ends | 49 of 54 | 30 of 36 | 34 of 36 |
| compleasm complete | see `examples/` | 99.43%, 99.64% | 99.68%, 99.71% |
| Connectors restored | no | no | yes |

No Juicebox review had been done on any at the time of writing. The user adopted Busk as it stands. **Do not read Busk's lead as the method alone:** it has over twice the depth, as well as connectors. Counted but not built: 28 of 71 WatsonsChoice gaps and 28 of 63 13Q010 gaps have a graph path through unused unitigs, against 40 of 60 on Busk.

## Dropped does not mean useless

A unitig is redundant when another copy covers its sequence, but it may still be the only graph link between two chains. Restoring those took Busk from 60 gaps to 8 before any read filling, with no bridges across chromosomes, and grouping unchanged. Their joins hold up about as well as joins inside chains: 88% have q20 reads through them at 10 kb, against 92% within chains. **What this does not show:** which allele a connector carries relative to its neighbours. Most such arms have arbitrary phase.

## A short haplotype can mean a tangle in unplaced sequence

chr7_h2 is 4.85 Mb shorter than chr7_h1. Both gap edges open into a connected component of 27 unplaced chains (7.24 Mb, mostly repeats), which touches no other scaffold in the graph, and read tails enter it. No single chain has enough Hi-C to place. Summed Hi-C was misleading at first: one chain with a repeat fraction of 0.90 gave 77% of its links and pointed to chr4. Without it, chr7 leads, though barely separating haplotypes (347 against 281). Compare homolog lengths, and look at what lies behind the gap of the shorter one.

## A low per-haplotype gene score can be numbering

13Q010's h1 and h3 looked far less complete than any diploid, and 917 gene copies appeared lost from their chromosome's homolog set. In 784 of 843 two-copy genes the third copy lay on another chromosome, in blocks matching the translocations already called. Haplotype numbers are arbitrary across chromosomes, so "all h1 scaffolds" need not be one set. Choosing the numbering per chromosome that completes the most genes gave 99.32, 99.50 and 98.82%, from 91.84, 98.18 and 93.73%. The real shortfall is 65 gene copies of 8,415 (0.77%; Busk 0.21%). This holds only if the translocations are real, and the sets are not shown to be parental genomes. `validate/complement_genes.py` reports it.

## Read correction

Raw SUP reads against corrected reads, at reads of 5 kb or more:

| | WatsonsChoice | 13Q010 | Busk |
|---|---|---|---|
| Reads kept | 60% | 83% | 85% |
| Bases kept | 65% | 88.5% | 89% |
| q20 spanning reads kept, at gaps | 38% | 82% | 8%, but 37 of 40 reads at one array |
| q20 spanning reads kept, control sites | 92% | 96% | 96% |
| Bases kept, reads of 100 kb or more | 21% | 89% | 81% |
| q20 in both forms, placed differently | 215 of 1.04 M | 267 of 0.90 M | not measured |

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
- I recorded D14 on Busk before reading the new report's nominations, and had to reverse two.
- I reported "raw 40, corrected 3" spanners on Busk from a pooled count, when 37 of the 40 were at one array.
- A report line read `homolog fraction 0.0000` on every genome until Busk, because of a parser bug. A value that exact should have looked wrong sooner.

## Open questions

- **Why do reads stop at gaps?** See [IMPROVING](IMPROVING.md).
- **What is chain0192, 0205, 0212, 0222 on Busk?** About 1,500x, 514 kb, with evenly spread Hi-C. Organelle is the obvious guess, and it is untested.
- **Is collapsed handling sound above ploidy 2?** It is generalised in code, but only ploidy 2 has been run.
- **Juicebox edits:** the route back into the AGP has not been exercised.
