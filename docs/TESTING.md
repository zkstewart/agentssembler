# Testing

`tests/test_units.py` runs in seconds on small made-up data:

```bash
python3 tests/test_units.py
```

`tests/regress.py` replays a saved assembly through the driver, seeding costly mappings and replaying every hand decision, and compares outputs with those of the saved run. It needs a spec naming the seeds, decisions and comparisons.

**This is the blinded copy of the repository.** The specs, the regression results and the worked examples name the genomes the process was developed on and the decisions made for them, so they are kept on the `main` branch only.

## What the tests do not show

- **Decisions were replayed, not made.** A regression shows the driver does what the hand-made edits did. It says nothing about whether an agent following the guides would reach the same decisions.
- **Not run in regression:** compleasm through the driver, connectors, and any Juicebox review route.
- Thresholds derive from `hap_depth`, and have been tried at two depths only. Three related genomes are not a proof of generality.
