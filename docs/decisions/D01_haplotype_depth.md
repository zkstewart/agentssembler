# D01. What is the haplotype depth, and is there collapsed sequence?

**Sets:** `hap_depth`, `repeat_depth`, `collapsed`. Every later depth threshold is derived from `hap_depth`, so this is the most consequential number in the run.

**Read:** the histogram in `reports/03_depth.md`, `04_copynumber/depth.windows.tsv`, `04_copynumber/class_depth.tsv`, and `00_inputs/unitigs.gfa` for graph degree.

**Decide:**

1. **Haplotype depth is the main peak of the histogram.** Do not derive it from genome size, and do not assume one peak.
2. **Check by adding depths of near-identical copies.** Copies summing to one haplotype's depth are one haplotype represented twice, not shared sequence. Class depth, not a unitig's own depth, gives copy number.
3. **`repeat_depth`** sits above any plausible copy number: about 4.4 x haplotype depth worked on both genomes.
4. **`collapsed: true`** when there is real sequence near twice haplotype depth **in long unitigs**, usually with graph degree 2,2. That sequence is homozygous and merged in the graph, so it belongs in more than one haplotype. It turns on restored arms, shared copies, masked mapping and joint arm choice. There need not be a visible peak: it can be about 1% of the sequence and invisible in the histogram, while unitigs of hundreds of kb at 2x with bubbles on both sides are unambiguous. The report lists them.

**A low shelf is not automatically low coverage.** On one genome a quarter of bases sat below the main peak in the primary histogram, and the shelf vanished in class depth: two copies at half the peak are one haplotype. Stage 04 removes them. Check the histogram for unitigs of 1 Mb or more too, since they show the peak with least of this noise.

**Traps:**

- GFA lists each link in both directions, so naive degree counts double.
- Arms overlap flanks by much of their length, so median arm depth is depressed. Do not read arm depth without excluding overlaps.
- Satellite arrays can sit at about 0.6 x haplotype depth. A collapsed satellite would read as 1x, so depth alone cannot find collapse there.
- A reference hit does not give copy number. Once, 20 satellite unitigs all hit one repeat window of the reference at low identity, and it looked like 4 copies of a chromosome. Check span against unitig length.

**Examples:**

| Histogram | Long unitigs near 2x | `hap_depth` | `collapsed` |
|---|---|---|---|
| One peak | None | The peak | false |
| Two peaks, the second at twice the first | Yes | The first peak | true |
| One peak, a low shelf that is 1x as a class, no visible 2x peak | Yes, about 1% of the sequence | The peak | true |

I once took the peak as total depth and divided it by the ploidy. That was wrong, and it came from reasoning about coverage instead of reading the histogram.

**Record:** `agentssembler decide WORK D01 --by NAME --reason "..." --set hap_depth=40 --set repeat_depth=176 --set collapsed=true`

**Changing it later:** `reset 04_chains`.

**Untested:** collapsed handling in a sample with more than two haplotypes. `mask_shared.py` and the fill stage accept any number, but only ploidy 2 has been run.
