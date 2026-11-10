# Record

**Not needed to run the pipeline. An agent running an assembly can skip this file.**

This is the evidence behind the process, kept for anyone assessing it: where following the guides exposed problems, which rules survived being doubted, the experiments behind entries in [IMPROVING](IMPROVING.md), and the errors made on the way. Chain and scaffold names apply to the build named, and to no other.

## Busk: where following the guides disagreed with the process

Busk was the first genome run by following the guides and making each decision, not replaying earlier ones. Each row is a point where following the repository as it stood would have given a worse result, or where it failed. Every change listed is now in the repository. Busk's decisions are in [`examples/`](../examples/README.md).

| # | What happened | Result if followed as written | Change |
|---|---|---|---|
| 1 | `build_chains.py` crashed on an overlap longer than a unitig | Stopped at stage 04 | Overlap carries into the next unitig. WatsonsChoice had one in the safe direction by luck, and 4 more unwalked |
| 2 | That crash left a partial AGP (21 of 266 chains), and resume treated it as finished | A 70-minute mapping on wrong input, wasted | `step` moves failed output aside; `r_chains` compares chain counts between stages |
| 3 | **"Redundant" unitigs were the graph path between chain ends.** 40 of 60 gaps had one | 60 gaps. With them: 8 before filling | `restore_connectors.py`, on by default |
| 4 | `hops: 8` missed a 9-chain allele path once connectors were present: one false duplicate, two chains unplaced | Those errors | Default 12; 16 gave the same |
| 5 | Chains ran through collapsed forks, so 2x hubs went to one haplotype only | 13 of 14 h2-only genes on chr2_h2, 10 of them in one hub | `review/chain_forks.py`, `review/embedded_runs.py`, D04 fork section. h1 missing genes went from 19 to 7 |
| 6 | D14 accepted nominations with too little measured interior | Two wrong shared chains (under 5 kb interior), which I recorded first and reversed | D14: require interior; `r_group` lists nominations as candidates |
| 7 | D17 was asked before its input existed | An uninformed decision | Moved to the stage 10 gate; guide states the exact command |
| 8 | Group edits named labels, which are numbered by size and shift after a break | Old edits would have misplaced chain0082 and the chr2 hub. Caught by comparing with archived groups | `with:<chain>` edits |
| 9 | Subtelomeric satellites grouped by Hi-C with the wrong chromosome | Wrong telomere count per group | D06 telomere ledger; `r_group` alarm |
| 10 | D01 expected a visible 2x peak; Busk has none, and a 4-40x shelf instead | Misreading as a low-coverage class | D01, D05 wording; report gives histograms |
| 11 | compleasm's tools were not on PATH and the stage carried on | Validation finishing without gene content | Tools checked; failure stops the stage |
| 12 | The report printed `Specific-pair homolog fraction: 0.0000`. The parser skipped every row starting with `c`, so every `chr` row. **The same bug would have hit every genome** | A perfect-looking false number. Real values are 3.2-4.8% | Fixed |
| 13 | Merged joins were never read-tested at final coordinates | New joins taken on trust | `validate/join_reads.py`: 671 of 741 have q20 reads; connector joins 116 of 132, others 555 of 609 |
| 14 | Nothing asked where reads go past a gap edge. 14 of 14 q20 tails at chr7_h2 enter chain0099 | A 7.24 Mb tangle unnoticed behind a 4.85 Mb short haplotype | `review/gap_tails.py`, `review/unplaced_components.py`; D11 |
| 15 | A low-confidence chain, chain0042, stayed low-confidence by Hi-C alone | Unresolved 4.6 Mb | D06: walk the graph from the chain's gap end first. It settled the question |
| 16 | A pooled "raw 40, corrected 3 q20 spanners" looked like a strong correction effect. 37 of the 40 are at one array | A wrong claim about correction, which I made in passing | Per-site reading, already a rule. Followed late |

### Where the process held

- **The 5 kb overlap rule.** I doubted it for two genuine overlaps of 2.6 and 3.5 kb that were left as gaps. Raw reads closed one at exactly the overlap length, and the other has one read. But a 2.6 kb overlap also links chain0019 and chain0038, on homologs, with Hi-C 8x against. Keep the rule; curate with evidence.
- **Refusing array closures.** The chr6_h2 gap had 306 raw spanners, and length disagreed by 8-14 kb.
- **Evidence order.** Hi-C pointed the tangle at chr4 until one chain with a repeat fraction of 0.90 was excluded; the graph had pointed at chr7_h2 throughout.

## The first builds, before connectors

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

### What rebuilding with connectors changed

Both earlier genomes were rebuilt with connectors and adopted. WatsonsChoice went from 71 to 35 gaps (contig N50 14.60 to 20.43 Mb) and 13Q010 from 63 to 32 (16.80 to 24.68 Mb), with telomeres and gene completeness unchanged or better. Three effects to expect:

- **Orientation of blocks between Hi-C-only gaps can flip.** A connector sets the orientation by graph where Hi-C did before. On WatsonsChoice this withdrew two SV calls (chr1 19.5-23.3 Mb, chr5 11.1-11.8 Mb); on 13Q010 a 1.9 Mb chr8 block reversed. Re-test SV calls touching any block whose orientation changed.
- **Satellite stacks get placed.** On 13Q010, about 12.6 Mb of chr8 satellite chains, unplaced before, are now joined by graph bridges, and chr8_h1 grew to 48.2 Mb. A chain in the stack with no graph link (chain0192) was unplaced at D09 instead of being kept on Hi-C. Check the stack in Juicebox.
- **Haplotype numbers change.** Most 13Q010 scaffolds were renumbered (old chr7_h3 is new chr7_h1), though the content is the same. Map names by shared unitigs before comparing builds, and warn anyone with annotation keyed to old names.

A tip path whose junction reads are all MAPQ 0, because the tip is near-identical to another placed copy, is a user choice, not an automatic placement. On 13Q010 the user chose to place one (chain0136, chain0183 at chr5_h2), accepting possible mosaicism: it completed 17 genes to three copies and gave chr5_h2 its second telomere.

## Experiments

### Why reads stop at gaps

**Result, on three genomes: nothing physically stops reads there.** Use `review/edge_reads.py` and `review/edge_seq_depth.py` to repeat it on a new genome.

| Test | Result | So |
|---|---|---|
| Read ends and read starts at edges, against controls | No excess in the median, and starts rise with ends: WatsonsChoice 2.05 and 1.57 per kb at the edge, in the mean | Not a block on one strand, which could not produce read starts |
| Depth with distance | Fades over about 100 kb, q20 depth first: 31 at 100-200 kb to 2 at the edge, any MAPQ 34 to 5 | Not a hard edge. The flank becomes repetitive |
| Reads carrying 5 kb of flank, alone as target, at least 90% cover and 97% identity | Median 1.43-1.85 x haplotype depth against 1.25-1.56 at controls. Under 0.5 x at 2-14% of edges against 4-7% | Reads carrying this sequence exist at a normal rate. Windows over 10 x were left out as repeats |
| Where those reads are placed, WatsonsChoice | Here at q20 4.8% (controls 62.6%), here at low MAPQ 12.4% (5.8%), homolog 53.6% (30.1%), same scaffold elsewhere 19.2% (1.3%), unplaced 8.3% (0%) | Gap flanks have copies at 97% or more: on the homolog, nearby on the scaffold, and in unplaced sequence |

**What this means for more data:** the sequence is readable, so a different chemistry is not obviously needed. A read must anchor beyond a duplicated flank, and depth fades over about 100 kb, so read length is the likely limit. That is an inference: no longer library has been tried.

**Not settled:** whether reads from the gap locus itself are rare, or exist and go to a cleaner sibling copy. Depth at the homolog copy (median 1.88 x) barely differs from controls (1.75 x), and that control figure is inflated by shared sequence, so this test could not tell.

**My errors on the way:** a loose criterion counted transposon fragments and gave 0 of 259 low-depth edges as depleted, which could not have come out otherwise. Busk has only 7 edges left in the strict table.

### Phasing arms

**Reads can say something at a minority of junctions, and it is built, as a report.** `review/phase_reads.py` runs in stage 13 with no new mapping. Busk: of 32 arbitrary junctions, 14 confirmed, 2 conflict, 16 no reads or thin. WatsonsChoice: of 102, 27 confirmed, 1 mixed, 74 no reads or thin. Specific sequence is typically hundreds of kb from the arm, and reads of that length are rare (WatsonsChoice: 1,435 of at least 250 kb). The conflicts have 2 reads each, and no arm has been swapped.

**Calling heterozygous sites is not needed for this,** though I first described it that way: reads across the arm into specific sequence already carry the signal.

**Hi-C cannot phase arms,** at least on Busk. Hi-C does separate homologs when pairs exist: chains it grouped favour their own scaffold by a median of 6.8-8.4 x. But arms get under 20 pairs at 31 of 32 Busk and 89 of 102 WatsonsChoice junctions, and the cause is not the filter. Of 3.17 M Busk records on arm chains, 99.6% tie exactly between copies, 293 show any margin, and 388 are lost to `NM`.

### Allele-aware Hi-C for hard unplaced chains: no

**Result: it cannot place 13Q010's 958 kb unplaced chain (chain0135 in the first build, chain0136 in the connector rebuild).** Scaffold names below are the first build's; its chr7_h3 is chr7_h1 in the rebuild. The user later chose to place the chain at chr5_h2's start by graph alone, accepting possible mosaicism ([LESSONS](LESSONS.md)). I expected about 3% divergence from its copy. That was wrong: it is 99.97% identical over 928 kb to the first 957 kb of chr7_h3 (`minimap2 -x asm20`, de 0.0003), and 99.2% to the end of chr4_h3. Both it and the chr7_h3 copy have about haplotype depth (primary 11.6-29.1 and 26.5-27.3, haplotype 27), so they look like two real copies, not a duplicated haplotig.

| Pairs with one end on the chain and the other elsewhere | Count | Top scaffold, against its share of length |
|---|---:|---|
| All | 184,632 | chr4_h3 15.8% (3.5%), chr7_h3 14.6% (3.6%): the copies' neighbours, from tied reads |
| Chain end MAPQ 0 | 182,366 (98.8%) | |
| Chain end favoured by an alignment score margin of at least 5, `NM` at most 2, mate MAPQ at least 20 | 185 | chr7_h3 7.6% (3.6%); no chromosome above 17% |
| Both ends MAPQ at least 20 | 49 | chr7_h3 8 pairs |

At 0.03% divergence a 150 bp read carries a difference about once in 20 reads, and the pairs that do carry one spread over all chromosomes. With the arm result on Busk (99.6% exact ties), the margin adds too few pairs to change any placement. Not built. The scratch tool was a stream of the Hi-C BAM through `awk` (9 minutes for 20 GB) and a tally by scaffold.

### A homolog synteny term in ordering: no effect

`synteny_links.py` scored a junction by whether the two chains' alleles lie close together on a homolog, and `hic_order.py --synteny` added a small bonus for each homolog in support. It used no reference. Homolog evidence ranks below Hi-C and above the reference, so a tie-breaking bonus was consistent with the evidence order. It changed no order over four rounds, and both were removed from the code.

## Errors made while building this

- I stated coverage from reasoning, not the histogram (9x per haplotype; it was 27x).
- I recorded a 5 Mb inversion on 13Q010 that does not exist, and another reached a curation file.
- I called reused arms legitimate before testing them.
- I said an ordering error explained an inversion. It survived correct ordering.
- A hard-coded `!= 3` in `dosage.py` listed the wrong windows on a diploid. The fractions it printed were right, but the list of runs was not. Now fixed.
- I gave a pooled corrected/raw depth ratio of 17-23% before finding two arrays had inflated it. The median is about 50%.
- Writing raw reads to an uncompressed FASTA filled the host disk and crashed the machine.
- I recorded D14 on Busk before reading the new report's nominations, and had to reverse two.
- I called chain0135/0136 about 3% divergent from its placed match. It is 0.03%.
- I accepted a 13Q010 D06 grouping without checking the chr8 satellite stack the guide names.
- I reported "raw 40, corrected 3" spanners on Busk from a pooled count, when 37 of the 40 were at one array.
- A report line read `homolog fraction 0.0000` on every genome until Busk, because of a parser bug. A value that exact should have looked wrong sooner.
