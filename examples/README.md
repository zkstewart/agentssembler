# Worked examples

Real hand decisions from the three genomes, to show the form and the depth of reasoning expected. **Chain names apply to one run only: never copy these into another sample.**

## 13Q010 (triploid, 9 chromosomes, 27x per haplotype)

| Decision | Outcome | File |
|---|---|---|
| D01 | `hap_depth` 27, `repeat_depth` 120, `collapsed` false | |
| D03 | One join broken by rule, none by hand | |
| D04 | `chain0010 15784630`: two signals agreed | |
| D07 | chain0221 attached to g16; chain0182 and chain0265 left unplaced | |
| D10 | `chain0129+ chain0081+`, Hi-C near a tie | `13Q010/curation.tsv` |
| D16 | No tails placed: two false copies found, the rest added no missing genes. The user chose minimal risk | `13Q010/arm_copies.tsv` |
| D18 | All 19 existing joins pass; one new closure kept, one rejected | `13Q010/join_changes.tsv` |
| D12 | SV calls with evidence | `13Q010/sv_calls.tsv` |

## WatsonsChoice (diploid, 9 chromosomes, 28x per haplotype)

| Decision | Outcome | File |
|---|---|---|
| D01 | `hap_depth` 28, `repeat_depth` 123, `collapsed` true | |
| D04 | No breaks; `cross_ratio` alone flagged six chains | |
| D05 | `hap_first`: chromosome-first split chr5 by haplotype and merged chr2 with chr6 | |
| D06, D07, D14 | 23 group edits, each with a reason | `WatsonsChoice/group_edits.tsv` |
| D11 | `chr6_h2` retried at `--minMapq 1`, closing one gap with 28 of 28 reads | |
| D16 | 47 tails placed before the stronger tests existed; those tests later found no tandem false duplication among them | |
| D18 | 7 closures kept, 1 rejected, 2 earlier closures reopened | `WatsonsChoice/rejected_joins.tsv` |
| D12 | SV calls, including ones withdrawn | `WatsonsChoice/sv_calls.tsv` |
| | Sites for the human review | `WatsonsChoice/curation_targets.tsv` |

## Busk (diploid, 9 chromosomes, 62x per haplotype)

The first genome decided by following the guides. [`Busk/what_changed.md`](Busk/what_changed.md) lists each decision and every place the process had to change.
