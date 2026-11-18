# Testing

Two levels. `tests/test_units.py` runs in seconds on small made-up data. `tests/regress.py` replays a saved assembly through the driver, seeding costly mappings and replaying every hand decision, and compares outputs with those of the saved run.

```bash
python3 tests/test_units.py
python3 tests/regress.py tests/specs/WatsonsChoice.yaml --saved SAVED_DIR --work WORK_DIR
```

## Regression results

**Both specs replay the first builds, which the connector rebuilds have since superseded.** They compare against the old `unitig_assembly` folders, with connectors off. That is still a valid test of the driver, but:

- Those folders are no longer the adopted assemblies. If they are deleted to free disk, the regressions cannot run.
- The adopted builds (both connector rebuilds, and Busk) have no spec. Specs replaying them would test the current defaults, and would be worth more than these two.

Last run after the Busk changes, with `restore_connectors: false` and `hops: 8` in both specs, since those builds predate connectors. Each reproduces the build it replays, not what the current defaults would now build.

| Genome | Compared | Same | Different |
|---|---|---|---|
| WatsonsChoice (diploid, collapsed, all 15 stages) | 34 | 33 | 1, harmless |
| 13Q010 (triploid, not collapsed, all 15 stages) | 22 | 20 | 2, both harmless |

WatsonsChoice reproduced the adopted `15_final_v3` AGP and FASTA index, and its validation telomere table, exactly. 13Q010 reproduced its scaffold AGP and FASTA and the adopted `15_final_v2` AGP and FASTA index exactly.

**Run after the D16 tangle check was added** (lost-route forks, `chain_forks.py` fix), through `11_build` only:

| Genome | Compared | Same | Different |
|---|---|---|---|
| GammaCitrus (diploid, collapsed; spec kept outside the repo, as the genome is confidential) | 32 | 30 | 2, expected |
| 13Q010 (`--missingReadsOk --to 11_build`) | 19 | 17 | 2, both the known ones below |

- **Both rebuilt the saved scaffolds exactly:** the stage 11 AGP, the scaffold FASTA and the unplaced FASTA.
- **GammaCitrus `10_groups/chain_forks.txt`:** the same 100 fork rows. Rows to chains in no group are now labelled `ungrouped`, not `same_group` or `other`, and there is a new `break_at` column.
- **GammaCitrus `12_shared/arm_usage.tsv`:** the same 112 rows in a different order.
- **The new stage 10 scan ran only on GammaCitrus,** since it needs a collapsed sample. It wrote `12_shared/chain_forks.tsv` with 12 `LOST` pairs and raised the stage 10 alarm. 13Q010 does not exercise it.
- **13Q010's inputs are mostly gone.** Its corrected reads exist only gzipped elsewhere, and its raw BAM and Hi-C reads are deleted. `init` refuses a missing input, so the run used a local copy of the spec with empty Hi-C stand-ins. That is safe because `09_hic_map/pairs.tsv.gz` is seeded, so stage 05 maps nothing. The WatsonsChoice saved build is deleted, so that spec cannot run.

Differences understood:

- WatsonsChoice `07_chains2/restored_arms.tsv`: one arm overlaps two flanks by the same length, and the table names a different one of the pair. Nothing downstream changed.
- 13Q010 `12_shared/order_links.tsv`: one extra link (`chain0139+ chain0221- direct`), because group edits now precede link writing.
- 13Q010 `12_shared/order.tsv`: statistics columns differ; chains, orientations and the resulting scaffolds are identical.

**If a saved build's reads are gone,** `--missingReadsOk --to 11_build` stands in an empty read file and seeds the saved bridge choices. Do not read bridge choices from such a run: stage 08 chooses bridges with reads. Without seeding, an empty read file gave no read support at the 19 bridges, chose different arms and changed the scaffolds, which looked like a regression. With them seeded, 13Q010's scaffolds were identical.

## Ploidy above 2, and what has not been re-run

Three changes followed the citrus builds: `hap_first` at any ploidy, a dosage for shared chains, and depth bounds that rise with ploidy ([SETTINGS](SETTINGS.md)). **No regression has been run since.**

- **Unit tests** cover them on made-up data: `test_hap_first_any_ploidy` and `test_shared_dosage`.
- **Ploidy 2 should be unchanged.** The derived depth values are the same numbers, and `hic_group.py` gave byte-identical groups, chromosomes, shared table and log to the version before it on 60 simulated diploids (15 seeds, both grouping modes, with and without shared chains; `PYTHONHASHSEED` fixed, since rows of equal length are otherwise written in set order). A simulation is not a genome: run a diploid spec before relying on this.
- **13Q010 will differ.** At ploidy 3 the highest tail depth is now 4 x haplotype depth, not 3, and the `high` depth class starts at 3.68 x, not 2.68. Chains between those bounds can now be tail candidates. Re-run its spec and read the differences in `12_shared/tails.tsv`.
- **`min_within`** is pinned to 700 in both specs. On copies of four builds' `chain_stats.tsv`, the driver measured 727, 611, 708 and 3,025 (13Q010, WatsonsChoice, irma, Busk), matching the table in [SETTINGS](SETTINGS.md). `test_hic_density` covers the arithmetic.
- **The length scale** ([SETTINGS](SETTINGS.md)) is pinned to 1 in both specs, so they replay with the fixed windows their builds used. A new citrus build measures 0.97 to 1.03: on the four builds measured, that moved at most one chain of about 1 Mb across the seed length, so expect a near match with an old build, not an exact one. A work directory already past stage 05 keeps scale 1. `test_length_scale` covers the arithmetic only.
- **Ploidy 4 and 6** grouped, gathered and chose dosage haplotypes correctly on simulated data with clean Hi-C. That shows the code does what it says, and nothing about a real autopolyploid.

## Spike-in test of `hic_crosstalk.py`

`tests/crosstalk_spikein.py` asks whether the Hi-C score in `utilities/` finds errors known to be there. It writes four spiked copies of an assembly (nothing wrong; phase switches; misjoins and one inversion; duplication, fragmentation and removal), scores each, and checks every category against the truth. `demo` does this on a simulated genome in about 20 minutes on 16 threads. `mutate` and `evaluate` do it on a real assembly, which is the stronger test: the script's help gives the commands.

```bash
python3 tests/crosstalk_spikein.py demo OUT_DIR --annotarium /path/to/annotarium.py -t 16
```

Run on simulated data (diploid, 5 chromosomes, 132 Mb, 2% divergence, 6 M pairs of 150 bp, about 27,600 scored pairs per Mb), with `hic_crosstalk.py` at its defaults. All 12 checks pass:

| Spiked | Result |
|---|---|
| Nothing (error-free assembly) | C 100%, no bin flagged |
| Renamed, reordered, reverse-complemented | Identical to the above |
| Phase switches of 0.2, 0.5, 1 and 3 Mb, and one 6.2 Mb arm | All 5 found on both haplotypes, none also called a misjoin. Blocks up to 1 Mb are flagged whole; the 3 Mb block and the arm only for 0.1-0.3 Mb at each breakpoint |
| Misjoins of 0.2, 0.5, 1 and 3 Mb | All 4 found over the whole inserted segment, and at the site each was taken from. None also called a switch |
| A 3 Mb inversion in place | Not found. Order and orientation are not scored, as the README says |
| 6 Mb duplicated as an extra backbone | D rose by 18 Mb: all three copies of the region count, not the extra one. Nothing flagged |
| One backbone cut into 3 Mb pieces, one removed | O rose by their length and U by their homologs' length, exactly. Nothing flagged |

**The first run of this test failed two checks, and the scan was changed** (`SCORE_VERSION` 2). A bin used to be flagged on one ratio: contacts with one other sequence over its own far contacts. That ratio also rises when the own far contacts are merely missing:

- **A phase switch was also reported as a misjoin.** With only switches spiked, X was 6.52% switch and 5.23% misjoin, though contacts with other chromosomes were at the noise level in those bins.
- **The flanks of a long region with no unique mapping were flagged as misjoins:** 2.9 Mb beside the 6 Mb duplicate, with no misjoin present.

A bin must now also receive more than `--fold` times the median share of its contacts from the other sequence. Both checks pass, and every spiked event is still found. **A simulated triploid** (198 Mb, 9 M pairs, switches between two of the three haplotypes) also passes all 12 checks, with D rising by all four copies of the duplicated region. Every bin flagged now was flagged before, so X can only fall. Scores under the old rule stay in `xtalk_q<mapq>/`; new ones go to `xtalk_q<mapq>.v2/`.

**Limits the test showed, which the change does not remove:**

- **X is not the length that is wrong, and how much of a long switch it covers depends on the rest of the assembly.** The baselines are the assembly's own medians, so many errors raise them. In the diploid run, with a switch on every chromosome, the median share doubled (0.068 against 0.032), and 0.1 Mb per haplotype of the 6.2 Mb arm was flagged. In the triploid run, where a third of the backbones were untouched, 3.7-4.6 Mb of it was. An assembly with more errors is therefore judged against a more lenient threshold. Report the four baseline rows of `summary.tsv` beside the scores of assemblies being compared, and count events from the run files, where one switch gives a run on each haplotype, and at each end of a long block.
- **With a fixed `--genome_size`, extra backbone sequence comes off O,** so a duplication hides an equal length of missing sequence.
- **Depth matters.** On a small trial with about 3% of the scored pairs per Mb, under the old rule, the switch scan missed one haplotype of two switches: too few bins reached `--min_hom`, and the baseline was set by the affected bins.

**What this does not show:** the simulated genome has no repeats, and its Hi-C library is clean. Real data may behave worse, which is why the real-assembly route exists. It has not been run.

## What the regressions do not show

- **Seeded mappings.** Mappings were copied from saved runs, so `map_paf` was exercised on only a few mappings. On 13Q010, `06_chains` was also seeded, because the saved run predates a `build_chains.py` fix (567 chains against 568).
- **Decisions were replayed, not made.** The regressions show the driver does what the hand-made edits did. They say nothing about whether an agent following the guides would reach the same decisions. Busk is the first such test: see `examples/`.
- **Not run in regression:** compleasm through the driver, which only the Busk run exercised, and any Juicebox review route.

- **Connectors are not regressed.** Both specs turn them off, so only the Busk run exercised them.

Two genomes reproduced is not a proof of generality. Thresholds derive from `hap_depth`, and were only tried at about 28x until Busk (62x). Busk has no earlier graph-first build to compare with, so it shows the process runs and what an agent decided, not that the result is correct.
