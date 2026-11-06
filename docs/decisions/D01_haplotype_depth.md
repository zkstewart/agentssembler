# D01. What is the haplotype depth, and is there collapsed sequence?

**Sets:** `hap_depth`, `repeat_depth`, `collapsed`. Every later depth threshold is derived from `hap_depth`, so this is the most consequential number in the run.

**Read:** the histogram in `reports/03_depth.md`, `04_copynumber/depth.windows.tsv`, `04_copynumber/class_depth.tsv`, and `00_inputs/unitigs.gfa` for graph degree.

**Decide:**

1. **Haplotype depth is the main peak of the histogram.** Do not derive it from genome size, and do not assume one peak.
2. **Check by adding depths of near-identical copies.** Copies summing to one haplotype's depth are one haplotype represented twice, not shared sequence. Class depth, not a unitig's own depth, gives copy number.
3. **`repeat_depth`** sits above any plausible copy number: about 4.4 x haplotype depth worked on both genomes.
4. **`collapsed: true`** when there is real sequence near twice haplotype depth **in long unitigs**, usually with graph degree 2,2. That sequence is homozygous and merged in the graph, so it belongs in more than one haplotype. It turns on restored arms, shared copies, masked mapping and joint arm choice. There need not be a visible peak: on Busk it is about 1% of the sequence and invisible in the histogram, but 300-850 kb unitigs at 2x with bubbles on both sides are unambiguous. The report lists them.

**A low shelf is not automatically low coverage.** On Busk a quarter of bases sat at 4-40x in the primary histogram, and the shelf vanished in class depth: two copies at 30x are one haplotype at 60x. Stage 04 removes them. Check the histogram for unitigs of 1 Mb or more too, since they show the peak with least of this noise.

**Traps:**

- GFA lists each link in both directions, so naive degree counts double.
- Arms overlap flanks by much of their length, so median arm depth is depressed. Do not read arm depth without excluding overlaps.
- Satellite arrays can sit at about 0.6 x haplotype depth. A collapsed satellite would read as 1x, so depth alone cannot find collapse there.
- A reference hit does not give copy number. On WatsonsChoice 20 unitigs of GC-rich satellite all hit one repeat window at 82% identity, and it looked like 4 copies of chr2. Check span against unitig length.

**Examples:**

| Sample | Reads | Peaks | `hap_depth` | `repeat_depth` | `collapsed` |
|---|---|---|---|---|---|
| 13Q010 (triploid) | 26.1 Gb | 27x | 27 | 120 | false |
| WatsonsChoice (diploid) | 19.15 Gb | 27-28x and 56x | 28 | 123 | true |
| Busk (diploid) | 41.1 Gb | 62x, a shelf at 4-40x that is 1x as a class, no visible 2x peak | 62 | 273 | true |

On 13Q010 I first said 28x total and 9x per haplotype. That was wrong, and it came from reasoning about coverage instead of reading the histogram.

**Record:** `agentssembler decide WORK D01 --by NAME --reason "..." --set hap_depth=28 --set repeat_depth=123 --set collapsed=true`

**Changing it later:** `reset 04_chains`.

**Untested:** collapsed handling in a sample with more than two haplotypes. `mask_shared.py` and the fill stage accept any number, but only ploidy 2 has been run.
