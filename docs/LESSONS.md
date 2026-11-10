# What the genomes taught

Results that shaped the process. They come from one triploid and two diploid citrus genomes, so treat them as evidence, not laws. Each genome's decisions and final figures are in [`examples/`](../examples/README.md), and the fuller record behind these findings is in [RECORD](RECORD.md).

## Dropped does not mean useless

A unitig is redundant when another copy covers its sequence, but it may still be the only graph link between two chains. Restoring those (connectors) took Busk from 60 gaps to 8 before any read filling, with no bridges across chromosomes, and grouping unchanged. On the other two genomes it halved the gaps. Their joins hold up about as well as joins inside chains: 88% have q20 reads through them at 10 kb, against 92% within chains. **What this does not show:** which allele a connector carries relative to its neighbours. Most such arms have arbitrary phase.

## A rebuild changes more than the gaps

Seen when two adopted builds were rebuilt after a change to the method:

- **Orientation of blocks between Hi-C-only gaps can flip,** where a graph join now sets what Hi-C set before. Two SV calls were withdrawn for this reason. Re-test every SV call touching such a block.
- **Satellite stacks can become placed** by graph bridges: about 12.6 Mb on one chromosome. Check the stack in Juicebox.
- **Haplotype numbers change,** though the content is the same. Map names by shared unitigs before comparing builds, and warn anyone with annotation keyed to old names.

## Sequence set aside must stay in the Hi-C target, as decoys

Stage 04 sets high-copy chains aside (`07_chains2/repeat_chains.txt`), and the chloroplast is one of them. When it was missing from the Hi-C target on 13Q010, chloroplast Hi-C reads landed on two nuclear plastid insertions as if unique: 21,394 pairs in one 500 kb bin pair. Grouping merged chr1 with chr2 twice and split chr9 to keep the count. `05_hic_chains` now maps Hi-C to the chains plus every set-aside unitig as `decoy_*`, and drops pairs touching a decoy. **Signal to look for:** thousands of pairs between two chains at one point (`awk` over `pairs.tsv.gz` in 500 kb bins).

## A low per-haplotype gene score can be numbering

13Q010's h1 and h3 looked far less complete than any diploid, and 917 gene copies appeared lost from their chromosome's homolog set. In 784 of 843 two-copy genes the third copy lay on another chromosome, in blocks matching the translocations already called. Haplotype numbers are arbitrary across chromosomes, so "all h1 scaffolds" need not be one set. Choosing the numbering per chromosome that completes the most genes gave 99.32, 99.50 and 98.82%, from 91.84, 98.18 and 93.73%. This holds only if the translocations are real, and the sets are not shown to be parental genomes. `validate/complement_genes.py` reports it.

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
- **It is not the main cause of gaps.** Raw depth also falls towards a gap edge, and 48 of 83 WatsonsChoice gaps have no spanning raw read. See [IMPROVING](IMPROVING.md), "Why reads stop at gaps".
- **Neither read set is enough alone.** Raw alone closed 11 WatsonsChoice gaps against 6 corrected, with 3 in common. On 13Q010 raw alone closed 15 against 19, and raw reads were inconsistent at 5 gaps that corrected reads filled.

## Automated labels that were wrong

| Label | What happened |
|---|---|
| `reads_overlap` | Closed two tandem-array gaps on WatsonsChoice, which reads anchored in the array cannot support. Raw-read audit reopened them |
| `only` tail | chain0256 on 13Q010 was already on all three haplotypes |
| Read count at joins | The worst join on 13Q010 had 8 reads, all at MAPQ 0 |
| High group density | A 262x chain passed the group test on density alone |

## Errors worth learning from

Each error made while building this is told in the guide where it applies: D01 (coverage from reasoning), D10 and D12 (an inversion that does not exist), D14 (shared chains recorded too early), D15 (reused arms called legitimate untested), D18 (pooled ratios). Three more, with no guide of their own:

- I said an ordering error explained an inversion. It survived correct ordering.
- I accepted a D06 grouping without checking the satellite stack the guide names.
- A report line read `homolog fraction 0.0000` on every genome until a parser bug was found. A value that exact should have looked wrong sooner.

## Open questions

- **Are a gap locus's own reads rare, or do they go to a cleaner sibling copy?** See [IMPROVING](IMPROVING.md).
- **Is collapsed handling sound above ploidy 2?** It is generalised in code, but only ploidy 2 has been run.
- **Juicebox edits:** the route back into the AGP has not been exercised.
