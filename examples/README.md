# Worked examples

Real hand decisions from the three genomes, to show the form and the depth of reasoning expected. **Chain names apply to one run only: never copy these into another sample.**

## 13Q010 (triploid, 9 chromosomes, 27x per haplotype)

32 gaps, contig N50 24.68 Mb. This is a rebuild of an earlier adopted build, and haplotype numbers differ between the two: `sv_calls.tsv` gives the earlier name for each call. **`13Q010/sv_calls.tsv` and `13Q010/curation_targets.tsv` predate the templates** in `config/handover_templates/` and use other columns: take the reasoning from them, and the columns from the template or from the Busk and WatsonsChoice tables.

| Decision | Outcome | File |
|---|---|---|
| D01 | `hap_depth` 27, `repeat_depth` 120, `collapsed` false | |
| D04 | `chain0010 15784630`: two signals agreed | `13Q010/chain_breaks.txt` |
| D05, D06 | Groups wrong on the first pass (chr1 with chr2, chr9 split): chloroplast Hi-C reads on nuclear plastid insertions. Fixed by decoys in stage 05, not by edits | |
| D07 | chain0222 telomere tip attached; a 1.16 Mb tip path (chain0136, chain0183) left for the user, who chose to place it at D13 | `13Q010/group_edits.tsv`, `13Q010/order_links.tsv` |
| D09 | A satellite chain with no graph link, placed on 15 Hi-C pairs between two gaps, unplaced | `13Q010/order_edits.tsv` |
| D10 | The earlier build's hand link `chain0129+ chain0081+` (Hi-C near a tie) is now a graph bridge through a connector; no curation needed | |
| D16 | No tails placed: two false copies found, the rest added no missing genes. The user chose minimal risk | `13Q010/arm_copies.tsv` |
| D18 | A join the earlier build rejected was re-tested and kept | `13Q010/join_changes.tsv` |
| D12 | SV calls with evidence, one new inversion | `13Q010/sv_calls.tsv` |
| | Sites for the human review | `13Q010/curation_targets.tsv` |

## WatsonsChoice (diploid, 9 chromosomes, 28x per haplotype)

35 gaps, contig N50 20.43 Mb. Also a rebuild of an earlier adopted build.

| Decision | Outcome | File |
|---|---|---|
| D01 | `hap_depth` 28, `repeat_depth` 123, `collapsed` true | |
| D04 | No breaks; `cross_ratio` alone flagged six chains | |
| D05 | `hap_first`: chromosome-first split chr5 by haplotype and merged chr2 with chr6 | |
| D06, D07, D14 | Group edits, each with a reason | `WatsonsChoice/group_edits.tsv` |
| D11 | `chr6_h2` retried at `--minMapq 1`, closing one gap with 28 of 28 reads (earlier build) | |
| D16 | 47 tails placed before the stronger tests existed; those tests later found no tandem false duplication among them (earlier build) | |
| D18 | 3 chr8 joins rejected: two tandem arrays, one trim with no copy in the kept flank | `WatsonsChoice/rejected_joins.tsv` |
| D12 | SV calls; two withdrawn because connectors now set orientations that came from Hi-C alone | `WatsonsChoice/sv_calls.tsv` |
| | Sites for the human review | `WatsonsChoice/curation_targets.tsv` |

## Busk (diploid, 9 chromosomes, 62x per haplotype)

The first genome decided by following the guides, not by replaying earlier decisions: 7 gaps, contig N50 34.47 Mb. Where that exposed problems in the process is in [RECORD](../docs/RECORD.md).

| Decision | Outcome | File |
|---|---|---|
| D01 | `hap_depth` 62, `repeat_depth` 273, `collapsed` true | |
| D03 | 2 of 523 joins broken by rule, 3 weak joins kept | |
| D04 | 2 fork breaks | `Busk/chain_breaks.txt` |
| D06, D14 | Group edits, each with a reason | `Busk/group_edits.tsv` |
| D09, D10 | Two chains unplaced; three curated links where Hi-C was near a tie | `Busk/order_edits.tsv`, `Busk/order_links.tsv` |
| D16 | 4 candidates, none placed | |
| D11, D18 | 0 fills from corrected reads; 1 raw overlap closure, passing audit | |
| D12 | No inverted block over 400 kb; 2 translocation candidates and 1 dosage run not supported | `Busk/sv_calls.tsv` |
| D13 | The user adopted the assembly as it stands, and left chain0099 unplaced | |
| | Sites for the human review | `Busk/curation_targets.tsv` |
