# Record

**Not needed to run the pipeline. An agent running an assembly can skip this file.**

This is the evidence behind the process, kept for anyone assessing it: where following the guides exposed problems, which rules survived being doubted, the experiments behind entries in [IMPROVING](IMPROVING.md), and the errors made on the way. Chain and scaffold names apply to the build named, and to no other.

**Where things go.** This file is the one home for numbers and for what happened. Each other document keeps its own role and links here:

- an idea not yet tested is in [IMPROVING](IMPROVING.md);
- once tested, its evidence comes here, under the genome it was run on;
- a rule that follows goes into a decision guide or [TRAPS](TRAPS.md);
- a lesson that holds beyond one case gets one line in [LESSONS](LESSONS.md);
- IMPROVING keeps a single line with the outcome.

## Genomes run so far

Other documents name these genomes rather than count them.

| Genome | Ploidy, chromosomes | Haplotype depth | `collapsed` | What it tested first | Decisions |
|---|---|---|---|---|---|
| 13Q010 | 3, 9 | 27x | false | The graph-first build; ploidy 3; organelle decoys; gene sets across translocations | [`examples/13Q010`](../examples/README.md) |
| WatsonsChoice | 2, 9 | 28x | true | `hap_first` grouping; read correction against raw reads; the raw-read join audit | [`examples/WatsonsChoice`](../examples/README.md) |
| Busk | 2, 9 | 62x | true | The first run made by following the guides; connectors | [`examples/Busk`](../examples/README.md) |
| GammaCitrus | 2, 9 | – | true | Tangles, TTT and lost-route forks | Kept outside the repository (confidential) |
| irma | – | – | – | Hi-C density only ([SETTINGS](SETTINGS.md)) | – |
| 292 | 2, 9 | 84x | true | Hi-C allele rescue, collapsed hubs, unlocalized scaffolds, genes on the unitigs | Not yet adopted |

All are citrus, so the findings are evidence, not laws.

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

Both earlier genomes were rebuilt with connectors and adopted. WatsonsChoice went from 71 to 35 gaps (contig N50 14.60 to 20.43 Mb) and 13Q010 from 63 to 32 (16.80 to 24.68 Mb), with telomeres and gene completeness unchanged or better. Connector joins hold up about as well as joins inside chains: 88% have q20 reads through them at 10 kb, against 92% within chains. Which allele a connector carries relative to its neighbours is not shown, and most such arms have arbitrary phase. Three effects to expect:

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

### TTT on graph tangles (GammaCitrus)

**Result: one clear gain, from the smallest tangle.** TTT ran after the pipeline had finished, on the finished build. Six tangles were found with `ttt/tangle_explore.py` (anchors: unitigs of 300 kb or more at 0.6-1.4x haplotype depth). Reads were 186,250 corrected reads touching them, aligned with GraphAligner to the whole graph. Coverage was stage 03's depth, with `--median-unique` from D01. The method is in [TANGLES](TANGLES.md).

| Tangle | Unitigs | Mb | What TTT did | Carried |
|---|---|---|---|---|
| chr5 telomeric | 17 | 2.77 | Clean 2-2. Both pairings scored the same. Both haplotypes reached unplaced 1x chains leading uniquely to telomeric chains. Paths identical with 30 seeds | Yes: both haplotypes telomere to telomere |
| chr4 hub block | 31 | 4.69 | Same topology as the baseline. chain0177 reuse forced. 2 extra chains entered from inside hub chains | No |
| chr2c distal | 49 | 7.27 | Same tip assignment as D07, swap tied. 4 bubble alleles at directly joined chain ends | No: Hi-C misplaced them, and they would have replaced joins with gaps |
| chr2a | 56 | 6.89 | Not separable from chr2b by its boundaries. Alone, on a graph with chr2b's links cut, it preferred the unplaced allele chain0081 over the reused arm, 114 to 24 | Yes, as the weakest change: 2 chain breaks and an order link the bridges could not see |
| chr2a+chr2b | 254 | about 31.5 | No valid path. About 10 banned Z connections needed. MILP stalled at a 49% gap | No |
| chr3 telomeric | 30 | 4.16 | One-sided. Not run | – |

**Acceptance:**

- Every traversal adjacency was a GFA link.
- Multiplicity against depth disagreed often. TTT gave multiplicity 2 to 0.5-1.0x nodes both haplotypes had to cross, and its coverage scale hit its lower bound (60 of 72).
- Reads at q20 spanned 0-8 chain junctions each, so most joins rest on the graph alone, as the baseline's do.
- hifiasm `rd:i` coverage was worse: it set three real alleles to 0.

**Build:** a second work directory with the baseline's edits plus TTT rows.

**Carried, and the result:**

- chr5 on both haplotypes. The pairing was broken by Hi-C at 41 against 24 pairs.
- The chr2a swap, by two chain breaks.

**Effect:** telomere to telomere rose from 4 to 6 of 18. Complete genes per haplotype rose from 95.76 and 95.26% to 97.18 and 96.83%, with 40-44 of each haplotype's 46 new genes on the chr5 arms. Unplaced sequence fell by 4.5 Mb, and there were 2 more gaps. SVs were unchanged.

**Dropped:**

- The chr2c and chr4 additions. Grouped, Hi-C put them at the wrong place and dragged a chr2 block with them.

**The chr2a swap is weak.** TTT's support was GraphAligner reads on utg000989l, which turned out to be an exact copy of the arm's first ~130 kb (as is utg000659l). In raw reads, 157 q20 tails from the hub run ~28 kb past the duplicated segment into the arm, and none go into chain0081.

**What it exposed:** stage 04's redundancy removal dropped unitigs that real haplotype paths need (utg000989l, utg000659l, utg000601l). Connectors restored two of them only after a chain break. A path through the others stays a gap: chr5_h1 has one. See "Lost-route forks" below.

### Lost-route forks

**Result: the fork scan had seen both of GammaCitrus's lost routes, and labelled them harmless.** chain0073 runs through hub utg000537l, whose other exit leads through dropped utg000989l to chain0081. chain0107 runs through hub utg000467l, whose other exits lead through dropped utg000913l and utg000601l. Both chains were in no group. `chain_forks.py` compared "unassigned" with "unassigned" as the same group, and so read both as ordinary bubbles.

**Why not flag every fork to an ungrouped chain at stage 07:** on GammaCitrus that is 66 rows, most of them bubbles inside collapsed chains in the chr2 tangle. At stage 10, placement is known, and the rule becomes specific. A row is flagged `LOST` when all of these hold:

- the hub is 2x;
- the chain is placed or on a tail;
- the other exit reaches a chain that is neither placed nor on a tail, at about 1x and at least 50 kb.

| Genome | Stage 07 HUB rows, before and after the fix | LOST rows at stage 10 (distinct chain to target) |
|---|---|---|
| GammaCitrus | 2, 2 | 12, including chain0073 -> chain0081 (779 kb), chain0107 -> chain0173 (h2's chr5 route), and the alleles of the reused arms chain0093, chain0138 and chain0185 |
| Busk | 1, 1 | 3 |
| 13Q010 | 1, 1 | 0 |

`break_at` for chain0073 is 778960, the break made by hand in the TTT build. Of GammaCitrus's other rows, chain0083 -> chain0193 and chain0089 -> chain0240 are bubble arms at directly joined chain ends. Added as TTT proposed, Hi-C misplaced them, so the guide says to leave such rows. Not yet run through a full build: the rows are a check for D16, and each break is a judgement.

### Hi-C allele rescue and collapsed hubs (292)

292 is a diploid at 84x per haplotype, whose haplotypes are near-identical over much of chr1-chr4. Two builds: the first without the rescue, and a rebuild with it (agentssembler 0.3.4). chr5-chr9 came out byte-identical in both. 292 has not been adopted yet, so nothing below is a final figure.

**The rescue's gain and its cost.** Stage 05's signal flagged 150 Hi-C-blind chains, 27.1 Mb carrying 391 gene copies: chains with an allele partner and almost no unique Hi-C, which grouping leaves unassigned however good the rest is. The rescue placed about 26 Mb of them by allele-ambiguous Hi-C.

| | First build | Rebuild with the rescue |
|---|---|---|
| Assembly | 615.5 Mb | 638.2 Mb |
| Gaps | 24 | 48: the new ones are Hi-C joins between rescued chains |
| compleasm single-copy, h1 / h2 | 92.16% / 92.16% | 95.54% / 96.43% |
| Gene copies short of 2 on the scaffolds | 405 | 187 |
| Inverted homolog blocks | 9 | 13, mostly in rescued chr2 |

Most rescued haplotype calls were `arbitrary` or `allele_of_arbitrary`: Hi-C named the chromosome, rarely the haplotype. The two chr2 homologs ordered their 12 shared rescued chains differently, at Hi-C gaps with no reads.

**Collapsed hubs.** With the rescue on, `chain_forks.txt` had 30 HUB rows, against 2 in the first build. Depth of unique reads on each chain split them:

- three were true collapses (chain0102, chain0119, chain0090): a continuous block at about 2x in which every read is unique, beside 1x arms whose reads are MAPQ 0;
- the rest read about 1x on the chain itself: allele bubbles, with both alleles present as separate chains.

Breaking the three, then two more that the next round exposed (chain0101, chain0115), cost two `reset 06_evidence --remap` rounds. Run on stage 04's mapping, before D04, `review/hub_cores.py` finds every one of these hubs in one pass. It reproduces the hand-made breaks for chain0101, chain0102, chain0104 (chr5), chain0119 and chain0012 at the same positions.

Breaks must fall at unitig boundaries (`break_chain.py`). A break is useful only if the piece it frees is mostly 2x: without that test, chain0159 and chain0119 would have carried a 1x arm into the shared piece. At stage 06 of the rebuild, seven remaining cores lie inside single unitigs, where hifiasm collapsed both haplotypes. They cannot be broken with the current tools. Each such chain goes on one haplotype or is shared, a D14 choice.

**Rescue swaps.** Of the rebuild's 46 Hi-C gaps, 2 (chr1_h1 17.67 Mb, chr3_h1 10.25 Mb) would close if two chains swapped homologs: the graph linked them the other way from the rescue's Hi-C call. The first build had no gap at either site. At chr1_h1, no q20 read crosses into either allele, so reads cannot settle it. At chr3_h1 the first build's bridges had join reads, but in near-identical haplotypes a read cannot show which homolog it came from. Neither the graph nor Hi-C is proof there, so the choice went to the user.

**Rescue errors caught by hand, then built in:**

- chain0288, a 17x near-copy of the telomeric tip chain0078, was placed between a chain and that tip. The `redundant` call now leaves such chains out.
- 2x chains with a collapsed core were given one haplotype. The `collapsed` call now leaves breakable ones out.

**Unlocalized scaffolds.** Handing over unplaced chains whose chromosome is known (`unloc/`) raised compleasm from 92% to 97-98% per haplotype on the first build, without changing a scaffold. On the rebuild, 97.18% and 97.50%.

**Direction.** chr2_h1 came out written in reverse against chr2_h2. Found at stage 15, the flip cost `reset 12_fill`. The same check now runs at stage 11 (`review/homolog_direction.py`), from the unitig self-alignments, and flags chr2_h1 in 0.1 s.

**Bugs found afterwards** by an external review, in code these builds used:

- `allele_table.py` wrote the unitig-level strand for 1,051 of 3,396 allele blocks, mirroring projected rescue pairs in those blocks;
- allele conflicts were judged against the longer chain, missing 104 of them.

Both are fixed. Whether they explain the chr2 order disagreement is not yet tested.

### Read correction (WatsonsChoice, 13Q010, Busk)

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
- **It is not the main cause of gaps.** Raw depth also falls towards a gap edge, and 48 of 83 WatsonsChoice gaps have no spanning raw read. See "Why reads stop at gaps" above.
- **Neither read set is enough alone.** Raw alone closed 11 WatsonsChoice gaps against 6 corrected, with 3 in common. On 13Q010 raw alone closed 15 against 19, and raw reads were inconsistent at 5 gaps that corrected reads filled.

### Organelle Hi-C on nuclear insertions (13Q010)

Stage 04 sets high-copy chains aside (`07_chains2/repeat_chains.txt`), and the chloroplast is one of them. When it was missing from the Hi-C target on 13Q010, chloroplast Hi-C reads landed on two nuclear plastid insertions as if unique: 21,394 pairs in one 500 kb bin pair. Grouping merged chr1 with chr2 twice and split chr9 to keep the count. `05_hic_chains` now maps Hi-C to the chains plus every set-aside unitig as `decoy_*`, and drops pairs touching a decoy. **Signal to look for:** thousands of pairs between two chains at one point (`awk` over `pairs.tsv.gz` in 500 kb bins).

### Gene sets across translocations (13Q010)

13Q010's h1 and h3 looked far less complete than any diploid, and 917 gene copies appeared lost from their chromosome's homolog set. In 784 of 843 two-copy genes the third copy lay on another chromosome, in blocks matching the translocations already called. Haplotype numbers are arbitrary across chromosomes, so "all h1 scaffolds" need not be one set. Choosing the numbering per chromosome that completes the most genes gave 99.32, 99.50 and 98.82%, from 91.84, 98.18 and 93.73%. This holds only if the translocations are real, and the sets are not shown to be parental genomes. `validate/complement_genes.py` reports it.

## Automated labels that were wrong

| Label | What happened |
|---|---|
| `reads_overlap` | Closed two tandem-array gaps on WatsonsChoice, which reads anchored in the array cannot support. Raw-read audit reopened them |
| `only` tail | chain0256 on 13Q010 was already on all three haplotypes |
| Read count at joins | The worst join on 13Q010 had 8 reads, all at MAPQ 0 |
| High group density | A 262x chain passed the group test on density alone |

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
- I said `break_chain.py` cuts anywhere, so `hub_cores.py` suggested breaks inside unitigs that the stage would have refused. On 292, 12 of 15 suggested breaks were of that kind.
- A report alarm I added for those hubs used its list before creating it, and would have stopped the stage 05 report whenever a breakable hub existed.
