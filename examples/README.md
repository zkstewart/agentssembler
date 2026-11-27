# Worked examples

Real hand decisions from 13Q010, WatsonsChoice and Busk, to show the form and the depth of reasoning expected. **Chain names apply to one run only: never copy these into another sample.**

## 13Q010 (triploid, 9 chromosomes, 27x per haplotype)

The 2026 rebuild with 0.3.7, and the regression test for ploidy above 2: 27 scaffolds, 937.6 Mb placed, 32 gaps, 59 contigs, contig N50 24.68 Mb, telomeric ends 50 of 54, compleasm complete per set (best numbering, across the translocations) 99.22%, 99.50%, 99.47%. Gaps and contig N50 match the earlier adopted build; chain names and haplotype numbers do not. Where each build exposed problems in the process is in [RECORD](../docs/RECORD.md).

| Decision | Outcome | File |
|---|---|---|
| D01 | `hap_depth` 27, `repeat_depth` 120, `collapsed` false | |
| D03 | The weak-join rule removed 1 of 1,033 joins and missed the worst, which D04 broke | |
| D04 | `chain0011 15784630`: two signals agreed (the earlier build's break). Hi-C allele rescue left off: 0.9% of chain sequence | `13Q010/chain_breaks.txt` |
| D05, D06 | 27 groups, three per chromosome, chromosome-first; no edits. The reference alarm (chr5 in four groups) is the chr2/chr5 translocation | |
| D07 | chain0221 telomere tip attached by graph; a 1.29 Mb graph-only tip path (chain0135, chain0182, MAPQ 0 throughout) left for the user, who chose at D13 to place it, as in the earlier build | `13Q010/group_edits.tsv`, `13Q010/order_links.tsv` |
| D09 | A satellite chain with no graph link, placed on 15 Hi-C pairs between two gaps, unplaced | `13Q010/order_edits.tsv` |
| D10 | No hand links: the earlier build's link `chain0129+ chain0081+` is a graph bridge through a connector | |
| D16 | 3 candidates, none placed: none carries a short gene; one the user chose to leave unplaced | |
| D18 | 2 raw-read joins audited, both kept, one with trim cover 1.0 at 99.05% | |
| D12 | 7 translocations, 4 haplotype-specific inversions, a chr8 rearrangement against the reference shared by all three haplotypes; every call of the earlier build reproduced | `13Q010/sv_calls.tsv` |
| D13 | The user adopted the build, placed the chr5_h2 start path, and left one tail unplaced | |
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

The first genome decided by following the guides, not by replaying earlier decisions. This is the 2026 rebuild with 0.3.7: 12 gaps, 30 contigs, contig N50 33.12 Mb, 616.1 Mb placed, compleasm complete per haplotype 99.68%. The first build had 7 gaps and contig N50 34.47 Mb; where each build exposed problems in the process is in [RECORD](../docs/RECORD.md).

| Decision | Outcome | File |
|---|---|---|
| D01 | `hap_depth` 62, `repeat_depth` 273, `collapsed` true | |
| D03 | 2 of 522 joins broken by rule, 3 weak joins kept | |
| D04 | No chimeras; 1 hub break, seen in `hub_cores.d04.tsv` before grouping. Hi-C allele rescue on | `Busk/chain_breaks.txt` |
| D06, D14 | Group edits, each with a reason: three telomeric satellites, five short hubs shared, and the user's placements from D13 | `Busk/group_edits.tsv` |
| D09, D10 | A chain's two shared copies unplaced; three `overlap:N` links, two graph overlaps under the 5 kb cut and one identical overlap outside the graph | `Busk/order_edits.tsv`, `Busk/order_links.tsv` |
| D16 | 5 candidates, none placed; TTT found no solution for the chr7_h2 tangle | |
| D11, D18 | 0 fills from corrected reads; no raw-read join needed | |
| D12 | A rearrangement on chr8 against the reference, on both haplotypes (high); 2 candidates not supported; 2 size differences | `Busk/sv_calls.tsv` |
| D13 | The user adopted the rebuild with a code fix (identical end overlaps are graph links), placed chain0099, and closed a chr4_h1 gap through a homozygous path | |
| | Sites for the human review | `Busk/curation_targets.tsv` |
