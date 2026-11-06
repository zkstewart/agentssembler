# Ideas

## Tried on Busk, and adopted

- **Connectors:** dropped unitigs that are the graph path between chain ends. 60 gaps to 8. See [LESSONS](LESSONS.md). **Not yet applied to WatsonsChoice or 13Q010:** 28 of 71 and 28 of 63 of their gaps have such a path. Rebuilding is the user's decision, since it changes adopted assemblies.
- **Tails past a gap edge and unplaced components:** `review/gap_tails.py`, `review/unplaced_components.py`. They find what is behind a gap; they do not place it.

## Not yet tried

What has not been squeezed from the unitigs, reads and Hi-C. **None of these has been tested**, so each is a hypothesis with a test, ranked by the value of what the test would tell you for its cost.

| # | Idea | Would improve | Cost |
|---|---|---|---|
| 1 | Measure base accuracy, then polish | Base-level quality, never measured | Low to medium |
| 2 | Find why reads stop at gaps | Whether any gap can close with this library | Low |
| 3 | Telomere-bearing raw reads at bare ends | Telomeric ends | Low |
| 4 | Phase arms through collapsed flanks with long raw reads and Hi-C | Fewer arbitrary arms, less mosaicism | Medium |
| 5 | Allele-aware Hi-C for hard unplaced chains | Placement of chains such as 13Q010 chain0135 | Medium |
| 6 | Reads under 5 kb, and local assembly at gaps | A few more closures | Medium |
| 7 | hifiasm's `r_utg` graph at gap ends | A few more closures | Low |

## 1. Measure base accuracy, then polish

Neither assembly has a k-mer QV. Everything so far measures structure and gene content. The Hi-C reads are an independent short-read set, so Merqury can give a QV and a per-scaffold error map from data in hand. Hi-C's k-mer spectrum is uneven, so read it as a floor. The raw reads (dorado 1.3.0 SUP, 5.2.0 and 5.0.0 models) can then polish, with the masked targets already built for reads over shared sequence.

**Test:** QV before and after, on one chromosome first. Keep a polish only if QV and compleasm frameshifts both improve. The one fill taken from a raw read (WatsonsChoice `chr5_h1:10518828`, about 1% error) gets fixed at the same time.

## 2. Find why reads stop at gaps

Raw depth falls from about 32x to 4x at gap edges, with ordinary GC and complexity, and 48 of 83 WatsonsChoice gaps have no spanning raw read. Whether more sequencing could help depends on the cause. Existing PAFs can separate three of them:

- **Read ends pile up at one coordinate, on one strand:** a sequence the pore or polymerase cannot pass, such as a hairpin or a damage site. Compare strand-specific read-end positions at each edge with control sites.
- **Reads continue but stop aligning:** look at clipped tails, including alignments under 5 kb. So far only overhangs of at least 2 kb with an anchor of at least 5 kb were counted.
- **Sequence is missing from the graph:** map the tails of reads ending there to each other.

Extract the flanks too, and look for inverted repeats within about 100 kb of the edge.

**This decides whether the answer is more analysis or a different library,** so run it before spending on anything else at gaps.

**First data from Busk, not a test of this idea:** of 7 remaining gaps, 5 have at most 2 q20 tails per gap, so reads simply stop. Two q20 tails at chr8_h2 map back across their own gap on both strands, which fits a foldback or inverted repeat. One read at chr6_h2 (`674e9354`) runs 97 kb forward and then back over the same stretch.

## 3. Telomere-bearing raw reads at bare ends

WatsonsChoice lacks 6 telomeric ends of 36. 13Q010 lacks them on `chr5_h1`, `chr5_h3`, `chr6_h1`, `chr6_h3` and `chr7_h3`. Telomeres are at a read's very end, where correction trims hardest. Search raw reads directly for the motif at their ends, align the rest, and look for reads anchoring uniquely to a bare scaffold end. No search of raw reads for the motif has been done.

**Test:** at least 3 reads per end, unique anchors of at least 10 kb, and agreement in how far the sequence extends. Subtelomeric repeats are shared among ends (13Q010 chain0182), so MAPQ must carry the weight. Only corrected reads and unitigs have been searched.

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
