# Ideas

## Tried on Busk, and adopted

- **Connectors:** dropped unitigs that are the graph path between chain ends. 60 gaps to 8. See [LESSONS](LESSONS.md). Since applied to WatsonsChoice (71 to 35 gaps) and 13Q010 (63 to 32) by rebuilding, and both rebuilds were adopted.
- **Tails past a gap edge and unplaced components:** `review/gap_tails.py` (stages 12, 13), `review/unplaced_components.py` (stage 15), with alarms in the reports. They find what is behind a gap; they do not place it.

## Out of scope

The user's decisions, so that nobody re-proposes them as next steps:

- **Base accuracy and polishing.** Beyond what this process should cover. Assemblies here are structural, so say that QV was not measured.
- **Telomere-bearing raw reads at bare ends.** Refinement best left to manual effort after the assembly. Note the bare ends in the hand-over.

## Not yet tried

What has not been squeezed from the unitigs, reads and Hi-C. **None of these has been tested**, so each is a hypothesis with a test. **The process is already large: an idea joins it only if a test shows a gain worth the added size.** A test that says no is a good result, and it goes under "Not worth repeating".

| # | Idea | Would improve | Cost |
|---|---|---|---|
| 6 | Reads under 5 kb, and local assembly at gaps | A few more closures | Medium |
| 7 | hifiasm's `r_utg` graph at gap ends | A few more closures | Low |

## 2. Why reads stop at gaps: tested

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

## 4. Phase arms: tested

**Reads can say something at a minority of junctions, and it is built, as a report.** `review/phase_reads.py` runs in stage 13 with no new mapping. Busk: of 32 arbitrary junctions, 14 confirmed, 2 conflict, 16 no reads or thin. WatsonsChoice: of 102, 27 confirmed, 1 mixed, 74 no reads or thin. Specific sequence is typically hundreds of kb from the arm, and reads of that length are rare (WatsonsChoice: 1,435 of at least 250 kb). The conflicts have 2 reads each, and no arm has been swapped.

**Calling heterozygous sites is not needed for this,** though I first described it that way: reads across the arm into specific sequence already carry the signal.

**Hi-C cannot phase arms,** at least on Busk. Hi-C does separate homologs when pairs exist: chains it grouped favour their own scaffold by a median of 6.8-8.4 x. But arms get under 20 pairs at 31 of 32 Busk and 89 of 102 WatsonsChoice junctions, and the cause is not the filter. Of 3.17 M Busk records on arm chains, 99.6% tie exactly between copies, 293 show any margin, and 388 are lost to `NM`.

## 5. Allele-aware Hi-C for hard unplaced chains: tested, no

**Result: it cannot place 13Q010's 958 kb unplaced chain (chain0135 in the first build, chain0136 in the connector rebuild).** Scaffold names below are the first build's; its chr7_h3 is chr7_h1 in the rebuild. The user later chose to place the chain at chr5_h2's start by graph alone, accepting possible mosaicism ([LESSONS](LESSONS.md)). I expected about 3% divergence from its copy. That was wrong: it is 99.97% identical over 928 kb to the first 957 kb of chr7_h3 (`minimap2 -x asm20`, de 0.0003), and 99.2% to the end of chr4_h3. Both it and the chr7_h3 copy have about haplotype depth (primary 11.6-29.1 and 26.5-27.3, haplotype 27), so they look like two real copies, not a duplicated haplotig.

| Pairs with one end on the chain and the other elsewhere | Count | Top scaffold, against its share of length |
|---|---:|---|
| All | 184,632 | chr4_h3 15.8% (3.5%), chr7_h3 14.6% (3.6%): the copies' neighbours, from tied reads |
| Chain end MAPQ 0 | 182,366 (98.8%) | |
| Chain end favoured by an alignment score margin of at least 5, `NM` at most 2, mate MAPQ at least 20 | 185 | chr7_h3 7.6% (3.6%); no chromosome above 17% |
| Both ends MAPQ at least 20 | 49 | chr7_h3 8 pairs |

At 0.03% divergence a 150 bp read carries a difference about once in 20 reads, and the pairs that do carry one spread over all chromosomes. With the arm result on Busk (99.6% exact ties), the margin adds too few pairs to change any placement. Not built. The scratch tool was a stream of the Hi-C BAM through `awk` (9 minutes for 20 GB) and a tally by scaffold.

## 6. Reads under 5 kb, and local assembly at gaps

Reads under 5 kb were never mapped (WatsonsChoice: 13.3 of 42.7 Gb, and on 13Q010 8.8 of 39.3 Gb). They cannot span by the anchor rule, but they may extend a flank, or tile across a gap with longer reads. A local assembly of all reads near each gap, with hifiasm or Flye, would use that, and it needs the tool to be installed. At `no_span` gaps few reads overhang the edge, so expect little, and run it after idea 2.

**Test:** assemble raw reads around the 10 `no_span` gaps with the most edge reads. Accept nothing that fails the same audit as any other join ([D18](decisions/D18_raw_read_joins.md)).

## 7. hifiasm's raw-unitig graph at gap ends

Both builds used the processed `p_utg` graph. hifiasm's `r_utg` graph keeps bubbles that `p_utg` pops, so gap ends may connect there. It costs one minimap2 run of chain ends against `r_utg` and a graph walk. Any path found needs read support like any other. 13Q010 has an `r_utg.gfa` from its hifiasm run, unused so far.

## Not worth repeating

- **Loosening `readfill.py`:** `--minIdentity 0.95 --maxTrim 15000` gained 1 gap and lost 2.
- **More hops on 13Q010:** 8 hops give the same 33 bridges as 3. That is before connectors: on Busk, with connectors, 12 hops are needed.
- **Raw reads to place 13Q010 unplaced sequence:** one q20 link in total. The evidence there is Hi-C or nothing.
- **Allele-aware Hi-C (idea 5):** near-identical copies are near-identical in Hi-C reads too. See above.
- **A reference synteny term in ordering:** no effect over four rounds, and it breaks the evidence order anyway.
