# Ideas

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

## Tested, and what it showed

The numbers behind each entry here and under "Not worth repeating" are in [RECORD](RECORD.md).

**Why reads stop at gaps: nothing physically stops them.** On three genomes, read ends show no excess at gap edges, and reads carrying the flank's sequence exist at a normal rate. Depth fades over about 100 kb, q20 depth first, because the flank has copies at 97% identity or more on the homolog, nearby on the scaffold and in unplaced sequence, and reads are placed there. So the sequence is readable, and read length is the likely limit; no longer library has been tried. Not settled: whether reads from the gap locus itself are rare, or exist and go to a cleaner sibling copy. `review/edge_reads.py` and `review/edge_seq_depth.py` repeat the test on a new genome.

**Phasing arms with reads: at a minority of junctions.** Built as a report, `review/phase_reads.py`, in stage 13. The numbers and limits are in [D15](decisions/D15_arm_phase.md).

**Solving tangles with TTT: worth it for small tangles with boundaries, so it is now part of D16's tangle check.** On GammaCitrus it carried both chr5 haplotypes through a 17-unitig 2-2 tangle to telomeric arms the pipeline left unplaced, making 46 more genes complete per haplotype. Elsewhere it confirmed the pipeline, could not decide phase, or found no path. When and how: [TANGLES](TANGLES.md).

## Not worth repeating

- **Loosening `readfill.py`:** `--minIdentity 0.95 --maxTrim 15000` gained 1 gap and lost 2.
- **Fewer or more hops than 12:** 3 missed a real path on one genome, and 8 missed one on another once connectors were in the graph. 16 found the same as 12.
- **Raw reads to place unplaced sequence:** one q20 link in total on 13Q010. The evidence there is Hi-C or nothing.
- **Allele-aware Hi-C for unplaced chains and arms:** near-identical copies are near-identical in Hi-C reads too. For a 958 kb chain 99.97% identical to a placed copy, 98.8% of its pairs were MAPQ 0, and the 185 with any alignment score margin spread over all chromosomes. On arms, 99.6% of records tie exactly.
- **TTT as an early stage, before chains:** planned for after D01 on GammaCitrus, and not built. Its one gain came from one small tangle. The largest tangle (254 unitigs) was beyond it, and most other output confirmed the baseline or was bubble alleles that chains cannot express.
- **TTT on a large tangle or one with dead ends:** no valid path, as every route needs banned connections, and the MILP grew by about 340 MB a minute.
- **hifiasm's `rd:i` as TTT coverage:** it dropped real alleles to multiplicity 0. Use stage 03's corrected-read depth.
- **A homolog synteny term in ordering:** a small bonus for junctions whose alleles are adjacent on a homolog had no effect over four rounds. Removed from the code.
