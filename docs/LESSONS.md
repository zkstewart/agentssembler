# What the genomes taught

Results that shaped the process. They come from three related genomes, so treat them as evidence, not laws. **This is the blinded copy:** the genomes are not named, and figures that would identify a chain, a scaffold or a decision's outcome on them are rounded or left out.

## Dropped does not mean useless

A unitig is redundant when another copy covers its sequence, but it may still be the only graph link between two chains. Restoring those (connectors) cut one genome's gaps more than sevenfold before any read filling, with no bridges across chromosomes, and grouping unchanged. On the other two genomes it halved the gaps. Their joins hold up about as well as joins inside chains: 88% have q20 reads through them at 10 kb, against 92% within chains. **What this does not show:** which allele a connector carries relative to its neighbours. Most such arms have arbitrary phase.

## A rebuild changes more than the gaps

Seen when two adopted builds were rebuilt after a change to the method:

- **Orientation of blocks between Hi-C-only gaps can flip,** where a graph join now sets what Hi-C set before. Two SV calls were withdrawn for this reason. Re-test every SV call touching such a block.
- **Satellite stacks can become placed** by graph bridges: over 10 Mb on one chromosome. Check the stack in Juicebox.
- **Haplotype numbers change,** though the content is the same. Map names by shared unitigs before comparing builds, and warn anyone with annotation keyed to old names.

## Sequence set aside must stay in the Hi-C target, as decoys

Stage 04 sets high-copy chains aside (`07_chains2/repeat_chains.txt`), and the chloroplast is one of them. When it was missing from the Hi-C target on one genome, chloroplast Hi-C reads landed on two nuclear plastid insertions as if unique: over 20,000 pairs in one 500 kb bin pair. Grouping merged two chromosomes and split another to keep the count. `05_hic_chains` now maps Hi-C to the chains plus every set-aside unitig as `decoy_*`, and drops pairs touching a decoy. **Signal to look for:** thousands of pairs between two chains at one point (`awk` over `pairs.tsv.gz` in 500 kb bins).

## A low per-haplotype gene score can be numbering

On a genome with translocations, some haplotype sets looked far less complete than others, and hundreds of gene copies appeared lost from their chromosome's homolog set. For over 90% of the genes one copy short, the missing copy lay on another chromosome, in blocks matching the translocations already called. Haplotype numbers are arbitrary across chromosomes, so "all h1 scaffolds" need not be one set. Choosing the numbering per chromosome that completes the most genes brought every set to about 99%, from as low as 92%. This holds only if the translocations are real, and the sets are not shown to be parental genomes. `validate/complement_genes.py` reports it.

## Read correction

Raw SUP reads against corrected reads, at reads of 5 kb or more:

| | Genome 1 | Genome 2 | Genome 3 |
|---|---|---|---|
| Reads kept | 60% | 83% | 85% |
| Bases kept | 65% | 88.5% | 89% |
| q20 spanning reads kept, at gaps | 38% | 82% | 8%, but 37 of 40 reads at one array |
| q20 spanning reads kept, control sites | 92% | 96% | 96% |
| Bases kept, reads of 100 kb or more | 21% | 89% | 81% |
| q20 in both forms, placed differently | 215 of 1.04 M | 267 of 0.90 M | not measured |

- **Correction loses and shortens reads near gaps.** At the same quality, reads near a gap on genome 1 are kept about half as often (34% against 63% at mean Q20-25). Of 151 raw spanning reads at q20, 50 had no corrected version and 45 no longer spanned, 34 of them shortened by more than 10%.
- **It does not move reads between haplotypes** on any scale measured.
- **It is not the main cause of gaps.** Raw depth also falls towards a gap edge, and over half of genome 1's gaps have no spanning raw read. See [IMPROVING](IMPROVING.md), "Why reads stop at gaps".
- **Neither read set is enough alone.** On one genome raw reads alone closed about twice as many gaps as corrected reads alone, with few in common. On another they closed fewer, and were inconsistent at some gaps that corrected reads filled.

## Automated labels that were wrong

| Label | What happened |
|---|---|
| `reads_overlap` | Closed two tandem-array gaps, which reads anchored in the array cannot support. Raw-read audit reopened them |
| `only` tail | A chain so labelled was already on every haplotype of its chromosome |
| Read count at joins | The worst join on one genome had 8 reads, all at MAPQ 0 |
| High group density | A chain at 9x haplotype depth passed the group test on density alone |

## Errors worth learning from

Each error made while building this is told in the guide where it applies: D01 (coverage from reasoning), D10 and D12 (an inversion that does not exist), D14 (shared chains recorded too early), D15 (reused arms called legitimate untested), D18 (pooled ratios). Three more, with no guide of their own:

- I said an ordering error explained an inversion. It survived correct ordering.
- I accepted a D06 grouping without checking the satellite stack the guide names.
- A report line read `homolog fraction 0.0000` on every genome until a parser bug was found. A value that exact should have looked wrong sooner.

## Open questions

- **Are a gap locus's own reads rare, or do they go to a cleaner sibling copy?** See [IMPROVING](IMPROVING.md).
- **Is collapsed handling sound above ploidy 2?** It is generalised in code, but only ploidy 2 has been run.
- **Juicebox edits:** the route back into the AGP has not been exercised.
