# Testing

Two levels. `tests/test_units.py` runs in seconds on small made-up data. `tests/regress.py` replays a saved assembly through the driver, seeding costly mappings and replaying every hand decision, and compares outputs with those of the saved run.

```bash
python3 tests/test_units.py
python3 tests/regress.py tests/specs/WatsonsChoice.yaml --saved SAVED_DIR --work WORK_DIR
```

## Regression results

Last run after the Busk changes, with `restore_connectors: false` and `hops: 8` in both specs, since those builds predate connectors. Each reproduces the build it replays, not what the current defaults would now build.

| Genome | Compared | Same | Different |
|---|---|---|---|
| WatsonsChoice (diploid, collapsed, all 15 stages) | 34 | 33 | 1, harmless |
| 13Q010 (triploid, not collapsed, all 15 stages) | 22 | 20 | 2, both harmless |

WatsonsChoice reproduced the adopted `15_final_v3` AGP and FASTA index, and its validation telomere table, exactly. 13Q010 reproduced its scaffold AGP and FASTA and the adopted `15_final_v2` AGP and FASTA index exactly.

Differences understood:

- WatsonsChoice `07_chains2/restored_arms.tsv`: one arm overlaps two flanks by the same length, and the table names a different one of the pair. Nothing downstream changed.
- 13Q010 `12_shared/order_links.tsv`: one extra link (`chain0139+ chain0221- direct`), because group edits now precede link writing.
- 13Q010 `12_shared/order.tsv`: statistics columns differ; chains, orientations and the resulting scaffolds are identical.

**If a saved build's reads are gone,** `--missingReadsOk --to 11_build` stands in an empty read file and seeds the saved bridge choices. Do not read bridge choices from such a run: stage 08 chooses bridges with reads. Without seeding, an empty read file gave no read support at the 19 bridges, chose different arms and changed the scaffolds, which looked like a regression. With them seeded, 13Q010's scaffolds were identical.

## What the regressions do not show

- **Seeded mappings.** Mappings were copied from saved runs, so `map_paf` was exercised on only a few mappings. On 13Q010, `06_chains` was also seeded, because the saved run predates a `build_chains.py` fix (567 chains against 568).
- **Decisions were replayed, not made.** The regressions show the driver does what the hand-made edits did. They say nothing about whether an agent following the guides would reach the same decisions. Busk is the first such test: see `examples/Busk`.
- **Not run in regression:** compleasm through the driver, which only the Busk run exercised, and any Juicebox review route.

- **Connectors are not regressed.** Both specs turn them off, so only the Busk run exercised them.

Two genomes reproduced is not a proof of generality. Thresholds derive from `hap_depth`, and were only tried at about 28x until Busk (62x). Busk has no earlier graph-first build to compare with, so it shows the process runs and what an agent decided, not that the result is correct.
