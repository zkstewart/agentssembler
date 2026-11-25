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

## Waiting for a regression build

**Agreed or built, but each changes results, so it needs a full build compared with the current code on the same sample.** The next new assembly will be that test, together with the `shared-pairs` branch. Run both builds from the same edits, compare with `agentssembler compare`, `utilities/agp_compare.py` and `utilities/hic_crosstalk.py compare`, and explain each difference at its site before adopting.

| Change | Where | What to compare |
|---|---|---|
| **`shared-pairs` branch**: rescued pairs shared by two alleles go to one group by a stable hash, not always to the first | `make_shared_copies.py`, `project_rescued.py` (branch only) | Orders of shared and rescued chains; Hi-C gaps; `hic_score` |
| **Allele blocks clamped to a trimmed overlap**: today a block running into one is dropped | `allele_table.py` | Allele pairs, conflicts in grouping, groups |
| **Ordering search**: the cost of `climb()`, how the link bonus scales, and when low-information chains join the backbone | `hic_order.py` | Orders (D09), orientation LLRs, time per chromosome |
| **Read filling**: sturdier cut points, trimming with several components, faster all-vs-all, excision in collapsed builds | `readfill.py` | Gaps filled, fill audits (D11, D18), validation |
| **Raw-read `--minIdentity` for stage 13 fills** | `readfill.py` via stage 13 | Fills gained against wrong fills (D18 audit) |
| **Allele conflict as a fraction of the shorter chain** (built, on main since 97fbaa1; both chains had to pass before) | `hic_group.py` | Groups. On the 292 rebuild's inputs, 3 of 490 chains move before the rescue: chain0107 joins g10 (the side the D13 swap took it from), chain0194 and chain0228 become unassigned. The review-v2 trial ran with it; `hic_order` orders on the rebuild's inputs are unchanged |
| **`hic_group.py` agglomeration below O(N^2)** | `hic_group.py` | Groups identical; time on a genome with many more chains |
| **One unmasked raw-read pass to pick reads near gaps** | stage 13 | `neargap` reads, fills and audits identical; time |
| **`raw_phase_early` on by default** (built, off; on two 292 builds it matched stage 13 on 16 of 16 arms, including both conflicts; see RECORD) | stage 10 | Whether an early conflict changes a D15 choice on a new genome, and the time (about 2.7 h on 292, repeated on every reset to stage 10) |

Decided not to do (owner, 2026-11-24): projecting the stage 02 read mapping onto chains instead of mapping again ("not worth the risk of bugs or behaviour non-identical to a real mapping run"). Splitting chains inside a unitig, POA consensus for fills, MILP ordering and three-way swaps stay deferred until a genome needs them.

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
- **TTT as an early stage, before chains:** planned for after D01 on GammaCitrus, and not built. Its one gain came from one small tangle. The largest tangle (254 unitigs) was beyond it, and most other output confirmed the baseline or was bubble alleles that chains cannot express. RECORD: "TTT on graph tangles (GammaCitrus)".
- **TTT on a large tangle or one with dead ends:** no valid path, as every route needs banned connections, and the MILP grew by about 340 MB a minute. RECORD: "TTT on graph tangles (GammaCitrus)".
- **hifiasm's `rd:i` as TTT coverage:** it dropped real alleles to multiplicity 0. Use stage 03's corrected-read depth.
- **A homolog synteny term in ordering:** a small bonus for junctions whose alleles are adjacent on a homolog had no effect over four rounds. Removed from the code. RECORD: "A homolog synteny term in ordering: no effect".
