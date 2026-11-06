# Testing

Two levels. `tests/test_units.py` runs in seconds on small made-up data. `tests/regress.py` replays a saved assembly through the driver, seeding costly mappings and replaying every hand decision, and compares outputs with those of the saved run.

```bash
python3 tests/test_units.py
python3 tests/regress.py tests/specs/WatsonsChoice.yaml --saved SAVED_DIR --work WORK_DIR
```

## Regression results

| Genome | Compared | Same | Different |
|---|---|---|---|
| WatsonsChoice (diploid, all 15 stages) | 34 | 32 | 2, both harmless |
| 13Q010 (triploid, stages 1-14) | 22 | 20 | 2, both harmless |

WatsonsChoice reproduced the adopted `15_final_v3` AGP and FASTA index, and its validation telomere table, exactly. 13Q010 reproduced its scaffold AGP and FASTA and the adopted `15_final_v2` AGP and FASTA index exactly.

Differences understood:

- WatsonsChoice `07_chains2/restored_arms.tsv`: one arm overlaps two flanks by the same length, and the table names a different one of the pair. Nothing downstream changed.
- WatsonsChoice `16_validate/pairs.specific.tsv.gz`: same lines in a different order, from parallel lifting.
- 13Q010 `12_shared/order_links.tsv`: one extra link (`chain0139+ chain0221- direct`), because group edits now precede link writing.
- 13Q010 `12_shared/order.tsv`: statistics columns differ; the resulting scaffolds are identical.

## What the regressions do not show

- **Seeded mappings.** Mappings were copied from saved runs, so `map_paf` was exercised on only a few mappings. On 13Q010, `06_chains` was also seeded, because the saved run predates a `build_chains.py` fix (567 chains against 568).
- **Decisions were replayed, not made.** The regressions show the driver does what the hand-made edits did. They say nothing about whether an agent following the guides would reach the same decisions. Busk is the first such test: see `examples/Busk`.
- **Not run in regression:** compleasm through the driver, and any Juicebox review route.

Two genomes reproduced is not a proof of generality. Thresholds derive from `hap_depth`, and were only tried at about 28x until Busk (62x).
