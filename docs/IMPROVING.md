# Ideas

## Tried on Busk, and adopted

- **Connectors:** dropped unitigs that are the graph path between chain ends. 60 gaps to 8. See [LESSONS](LESSONS.md). **Not yet applied to WatsonsChoice or 13Q010:** 28 of 71 and 28 of 63 of their gaps have such a path. Rebuilding is the user's decision, since it changes adopted assemblies.
- **Tails past a gap edge and unplaced components:** `review/gap_tails.py` (stages 12, 13), `review/unplaced_components.py` (stage 15), with alarms in the reports. They find what is behind a gap; they do not place it.

## Out of scope

The user's decisions, so that nobody re-proposes them as next steps:

- **Base accuracy and polishing.** Beyond what this process should cover. Assemblies here are structural, so say that QV was not measured.
- **Telomere-bearing raw reads at bare ends.** Refinement best left to manual effort after the assembly. Note the bare ends in the hand-over.

## Not yet tried

What has not been squeezed from the unitigs, reads and Hi-C. **None of these has been tested**, so each is a hypothesis with a test. **The process is already large: an idea joins it only if a test shows a gain worth the added size.** A test that says no is a good result, and it goes under "Not worth repeating".

| # | Idea | Would improve | Cost |
|---|---|---|---|
| 2 | Find why reads stop at gaps | Whether any gap can close with this library | Low |
| 4 | Phase arms through collapsed flanks with long raw reads and Hi-C | Fewer arbitrary arms, less mosaicism | Count first: low. Building it: medium, plus permanent size |
| 5 | Allele-aware Hi-C for hard unplaced chains | Placement of chains such as 13Q010 chain0135 | Count first: low. Building it: medium, plus permanent size |
| 6 | Reads under 5 kb, and local assembly at gaps | A few more closures | Medium |
| 7 | hifiasm's `r_utg` graph at gap ends | A few more closures | Low |

## 2. Find why reads stop at gaps

Raw depth falls from about 32x to 4x at gap edges, with ordinary GC and complexity, and 48 of 83 WatsonsChoice gaps have no spanning raw read. Whether more sequencing could help depends on the cause. Existing PAFs can separate three of them:

- **Read ends pile up at one coordinate, on one strand:** a sequence the pore or polymerase cannot pass, such as a hairpin or a damage site. Compare strand-specific read-end positions at each edge with control sites.
- **Reads continue but stop aligning:** look at clipped tails, including alignments under 5 kb. So far only overhangs of at least 2 kb with an anchor of at least 5 kb were counted.
- **Sequence is missing from the graph:** map the tails of reads ending there to each other.

Extract the flanks too, and look for inverted repeats within about 100 kb of the edge.

**This decides whether the answer is more analysis or a different library,** so run it before spending on anything else at gaps.

**First data from Busk, not a test of this idea:** of 7 remaining gaps, 5 have at most 2 q20 tails per gap, so reads simply stop. Two q20 tails at chr8_h2 map back across their own gap on both strands, which fits a foldback or inverted repeat. One read at chr6_h2 (`674e9354`) runs 97 kb forward and then back over the same stretch.

## 4. Phase arms through collapsed flanks

102 of 139 WatsonsChoice bridge junctions and 27 of 47 tails have arbitrary phase. Reads of 100 kb or more are the ones correction lost most (21% of their bases kept), and they have never been used to phase. Next to a collapsed flank, a long raw read reaching a heterozygous site beyond the flank links an arm to the haplotype, where Hi-C has too little to say. Hi-C pairs carrying heterozygous sites can extend the same evidence, as in SNP-level phasers such as WhatsHap and HapCUT2.

**Test:** count, per arbitrary junction, the raw reads carrying a heterozygous site on both sides. If most junctions have at least 2 reads, `assign_arms2.py` gains an evidence term for those. If few do, this is settled, and the labels stand.

**This would also bear on the tail decision,** since an arbitrary tail with phase evidence is no longer arbitrary.

## 5. Allele-aware Hi-C for hard unplaced chains

Chains with a near-identical copy get MAPQ 0 Hi-C, and the pipeline's filters drop those pairs. That is why 13Q010 chain0135 (958 kb, 21 genes, 18 of them missing from both h1 and h2) has so little evidence. It differs from its copy by about 3%, so many pairs do carry a difference. Competitive mapping to the two copies with a score margin, in place of MAPQ, would recover pairs that are currently discarded.

**Test:** on chain0135 alone, count pairs whose mate is unique elsewhere and whose read favours chain0135 by at least 2 mismatches. Compare the chromosomes hit with a shuffled control. Keep it unplaced unless one chromosome clearly stands out.

## 6. Reads under 5 kb, and local assembly at gaps

Reads under 5 kb were never mapped (WatsonsChoice: 13.3 of 42.7 Gb, and on 13Q010 8.8 of 39.3 Gb). They cannot span by the anchor rule, but they may extend a flank, or tile across a gap with longer reads. A local assembly of all reads near each gap, with hifiasm or Flye, would use that, and it needs the tool to be installed. At `no_span` gaps few reads overhang the edge, so expect little, and run it after idea 2.

**Test:** assemble raw reads around the 10 `no_span` gaps with the most edge reads. Accept nothing that fails the same audit as any other join ([D18](decisions/D18_raw_read_joins.md)).

## 7. hifiasm's raw-unitig graph at gap ends

Both builds used the processed `p_utg` graph. hifiasm's `r_utg` graph keeps bubbles that `p_utg` pops, so gap ends may connect there. It costs one minimap2 run of chain ends against `r_utg` and a graph walk. Any path found needs read support like any other. 13Q010 has an `r_utg.gfa` from its run, which the reference-guided gap filling may have read but the graph-first build did not.

## Not worth repeating

- **Loosening `readfill.py`:** `--minIdentity 0.95 --maxTrim 15000` gained 1 gap and lost 2.
- **More hops on 13Q010:** 8 hops give the same 33 bridges as 3. That is before connectors: on Busk, with connectors, 12 hops are needed.
- **Raw reads to place 13Q010 unplaced sequence:** one q20 link in total. The evidence there is Hi-C or nothing.
- **A reference synteny term in ordering:** no effect over four rounds, and it breaks the evidence order anyway.
