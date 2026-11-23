# Ideas

**Not needed to run the pipeline.** What comes next and what is ruled out: ideas not yet tested, open questions, and one line for each idea already tested. Evidence goes to [RECORD](RECORD.md) once a test is run.

## Out of scope

The user's decisions, so that nobody re-proposes them as next steps:

- **Base accuracy and polishing.** Beyond what this process should cover. Assemblies here are structural, so say that QV was not measured.
- **Telomere-bearing raw reads at bare ends.** Refinement best left to manual effort after the assembly. Note the bare ends in the hand-over.

## Not yet tried

**Neither of these has been tested**, so each is a hypothesis with a test. **The process is already large: an idea joins it only if a test shows a gain worth the added size.** A test that says no is a good result, and it goes under "Not worth repeating".

| Idea | Would improve | Cost |
|---|---|---|
| Reads under 5 kb, and local assembly at gaps | A few more closures | Medium |
| hifiasm's `r_utg` graph at gap ends | A few more closures | Low |

### Reads under 5 kb, and local assembly at gaps

Reads under 5 kb are never mapped (WatsonsChoice: 13.3 of 42.7 Gb, and on 13Q010 8.8 of 39.3 Gb). They cannot span by the anchor rule, but they may extend a flank, or tile across a gap with longer reads. A local assembly of all reads near each gap, with hifiasm or Flye, would use that, and it needs the tool to be installed. At `no_span` gaps few reads overhang the edge, so expect little.

**Test:** assemble raw reads around the 10 `no_span` gaps with the most edge reads. Accept nothing that fails the same audit as any other join ([D18](decisions/D18_raw_read_joins.md)).

### hifiasm's raw-unitig graph at gap ends

The process uses the processed `p_utg` graph. hifiasm's `r_utg` graph keeps bubbles that `p_utg` pops, so gap ends may connect there. It costs one minimap2 run of chain ends against `r_utg` and a graph walk. Any path found needs read support like any other.

## Open questions

- **Are a gap locus's own reads rare, or do they go to a cleaner sibling copy?** See "Why reads stop at gaps" in [RECORD](RECORD.md).
- **Is collapsed handling sound above ploidy 2?** Only ploidy 2 has been run. Above it, a shared chain gets a dosage and goes to that many haplotypes ([D14](decisions/D14_shared_chains.md)), and arm and tail targets rise to the ploidy. Unit tests on made-up data cover this. Whether Hi-C, alleles and graph links choose the right haplotypes on a real autotetraploid, with its short chains, is not known.
- **Does anything work outside heterozygous citrus?** `ploidy: 1`, allopolyploids and mostly collapsed genomes have config recipes ([SETTINGS](SETTINGS.md)) and no run behind them.
- **Juicebox edits:** the route back into the AGP has not been exercised.

## Tested, and what it showed

One line each: what was decided. The numbers are in [RECORD](RECORD.md), under the section named.

- **Why reads stop at gaps:** nothing physically stops them. Flanks have near-identical copies elsewhere, so read length is the likely limit; no longer library has been tried. RECORD: "Why reads stop at gaps".
- **Phasing arms with reads:** possible at a minority of junctions. Built as a report, `review/phase_reads.py`, in stage 13 ([D15](decisions/D15_arm_phase.md)). RECORD: "Phasing arms".
- **TTT on small tangles with boundaries:** a gain on GammaCitrus, so it is part of D16's tangle check ([TANGLES](TANGLES.md)). RECORD: "TTT on graph tangles".
- **Hi-C allele rescue for Hi-C-blind chains:** a gain in genes placed on 292, at the cost of Hi-C gaps and arbitrary haplotype calls. Built as the `hic_allele_rescue` setting, decided at [D04](decisions/D04_chimeric_chains.md). RECORD: "Hi-C allele rescue and collapsed hubs (292)".

## Not worth repeating

The numbers are in [RECORD](RECORD.md) where a section is named; the others are small enough to stand here.

- **Loosening `readfill.py`:** `--minIdentity 0.95 --maxTrim 15000` gained 1 gap and lost 2.
- **Fewer or more hops than 12:** 3 missed a real path on one genome, and 8 missed one on another once connectors were in the graph. 16 found the same as 12.
- **Raw reads to place unplaced sequence:** one q20 link in total on 13Q010. The evidence there is Hi-C or nothing.
- **Allele-aware Hi-C for unplaced chains and arms:** near-identical copies are near-identical in Hi-C reads too; the few reads with a margin spread over all chromosomes. RECORD: "Allele-aware Hi-C for hard unplaced chains", "Phasing arms".
- **TTT as an early stage, before chains:** planned for after D01 on GammaCitrus, and not built. Its one gain came from one small tangle. The largest tangle (254 unitigs) was beyond it, and most other output confirmed the baseline or was bubble alleles that chains cannot express.
- **TTT on a large tangle or one with dead ends:** no valid path, as every route needs banned connections, and the MILP grew by about 340 MB a minute.
- **hifiasm's `rd:i` as TTT coverage:** it dropped real alleles to multiplicity 0. Use stage 03's corrected-read depth.
- **A homolog synteny term in ordering:** a small bonus for junctions whose alleles are adjacent on a homolog had no effect over four rounds. Removed from the code.
