# The work directory

What a work directory holds, which stage writes each part, and what you add by hand. The directory's own name and place are yours to choose. Below, `S` stands for `sample` in the config, `hK` for a haplotype (`h1`, `h2`, ...), and `chrN` for a chromosome name.

## Read this first: directory numbers are not stage numbers

Stages are numbered 01 to 15. Directories are numbered 00 to 17. **The two sets of numbers are unrelated.** `06_chains/` is written by stage `04_chains`; there is no `06_evidence/` directory, and stage `06_evidence` writes `12_broken/`, `11_ont_chains/` and more. Use this table, not the numbers.

| Stage | Writes into | Stops for |
|---|---|---|
| `01_inputs` | `00_inputs/` | |
| `02_map_unitigs` | `02_ont_map/`, `03_unitig_align/` | |
| `03_depth` | `04_copynumber/` | D01 |
| `04_chains` | `04_copynumber/copy_number.tsv`, `05_nonredundant/`, `06_chains/`, `07_chains2/`, `08_final_chains/` | D03 |
| `05_hic_chains` | `09_hic_map/`, `10_groups/chain_breaks.tsv` | D04 |
| `06_evidence` | `12_broken/`, `11_ont_chains/`, `10_groups/` (alleles, chain statistics, end links, shared candidates), `16_validate/telomeres.chains.tsv` | |
| `07_group` | `10_groups/` (`groups*`, `shared.txt`, `eval_groups.txt`, `chain_forks.txt`) | D05, D06, D14 |
| `08_bridges` | `10_groups/final*`, `12_shared/`, `13_scaffold/bridges/` | D07 |
| `09_order` | `12_shared/` (`order*`, `by_chrom/`) | D09, D10 |
| `10_arms_tails` | `12_shared/` (`choices.tsv`, `arm_*`, `tails*`, `pre/`), `13_scaffold/names.tsv`, `13_scaffold/repeat_chains.fasta` | D17, D15, D16 |
| `11_build` | `13_scaffold/` | |
| `12_fill` | `14_readfill/`, `15_final/` | D11 |
| `13_rawreads` | `17_raw_reads/` | D18 |
| `14_apply_joins` | `15_final_raw/` | |
| `15_validate` | `16_validate/` | D12, D13 |

Three things look like mistakes and are not:

- **There is no `01_` directory.** Hi-C is mapped once, to the chains, by stage `05_hic_chains` into `09_hic_map/`.
- **There are two `12_` directories**, `12_broken/` (stage 06) and `12_shared/` (stages 08 to 10).
- **`16_validate/` is begun at stage 06** (one telomere table) and finished at stage 15, after `17_raw_reads/`.

## Top level

| Path | Made by | Contents |
|---|---|---|
| `config.yaml` | `init` | A copy of the config you gave. This copy is the one every later command reads. `init` never replaces it: to change a setting, edit it and `reset` from the earliest stage affected |
| `README.md` | `init`, then you | Status of the build, from `config/handover_templates/README.md`. Keep it current, and complete it at D13 |
| `state/` | driver | `STAGE.done` per finished stage, holding the settings it ran with; `base_dir`, the directory `init` was run in, from which relative config paths resolve |
| `decisions/` | `decide` | `Dnn.yaml` per recorded decision: who, when, reason, settings, and a checksum of every edit file at that moment |
| `edits/` | `init`, then you | Every hand edit, each row with its reason. See below |
| `reports/` | driver | `STAGE.md` per finished stage. `02_map_unitigs.md` and `01_inputs.md` are short by design |
| `logs/` | driver | `STAGE.log` per stage, appended on each run |
| `archive/` | `reset` | Only after a `reset`. See [Re-running](#re-running) |
| `00_inputs/` to `17_raw_reads/` | stages | Below |

A copy of the config beside these, such as `sample.yaml`, is harmless. Do not add other files or folders at the top level without listing them in `README.md`.

## Files you write

Everything else is written by a stage, and **must not be edited by hand**.

| File | Decision | Notes |
|---|---|---|
| `edits/join_breaks.txt` | D03 | |
| `edits/chain_breaks.txt` | D04, D06 | |
| `edits/shared.txt` | D14 | Replaces the whole candidate list. Collapsed samples only |
| `edits/group_edits.tsv` | D05, D06, D07, D14 | |
| `edits/order_links.tsv` | D07, D09, D10 | |
| `edits/order_edits.tsv` | D09 | |
| `edits/names.tsv` | stage 10 | Replaces the automatic names |
| `edits/accepted_tails.tsv` | D16 | Rows copied unchanged from `12_shared/tails.tsv` |
| `edits/fill_retry.tsv` | D11 | |
| `edits/flip.txt` | D11, D12 | |
| `edits/join_decisions.tsv` | D18 | |
| `16_validate/sv/sv_calls.tsv` | D12 | Stage 15 puts the header there |
| `16_validate/curation_targets.tsv` | D13, D17 | Stage 15 puts the header there |
| `README.md` | D13 | |

`init` copies a template with a header for each `edits/` file from `config/edits_templates/`. A file holding only comment lines changes nothing.

## Marks used below

- **[C]** only when `collapsed: true`. **[N]** only when `collapsed: false`.
- **[ref]** only with `inputs.reference`. **[raw]** only with `inputs.raw_bams`.
- **[opt]** only when the optional tool is set: `plot_python`, `juicer_tools` or `compleasm`.
- **[edit]** only when the named edit file has rows.

Every mapping `X.paf.gz` made by the driver has a `X.paf.gz.log` beside it, and is `X.paf.gz.part` while running. A mapping you linked or copied in has no log. A file ending `.failed` is the partial output of a step that crashed: read the stage log, then delete it.

## `00_inputs/` (stage 01)

Symlinks to the inputs under fixed names, and two indexes. Raw BAMs are not linked: stage 13 reads them from the config's paths.

`unitigs.fasta`, `unitigs.gfa`, `ont_corrected_reads.fasta`, `hic_R1.fq.gz`, `hic_R2.fq.gz`, `unitigs.fasta.fai`, `reference.fasta` and `reference.fasta.fai` **[ref]**.

## `02_ont_map/`, `03_unitig_align/` (stage 02)

- `02_ont_map/ont_vs_unitigs.paf.gz`: corrected reads to all unitigs.
- `03_unitig_align/unitigs_self.paf.gz`: unitigs against themselves.
- `03_unitig_align/unitigs_vs_ref.paf.gz` **[ref]**.

These are the mappings most often made elsewhere and put in place. See [HPC.md](HPC.md).

## `04_copynumber/` (stages 03 and 04)

`depth.windows.tsv`, `depth.unitigs.tsv`, `class_depth.tsv` (stage 03, read at D01). `copy_number.tsv` (stage 04, after D01).

## `05_nonredundant/`, `06_chains/`, `07_chains2/`, `08_final_chains/` (stage 04)

Three rounds of chains. Each round has `chains.agp`, `chains.fasta`, `chains.fasta.fai`, `chains.joins.tsv`, `chains.contained.tsv`, `chains.circular.txt`.

- `05_nonredundant/kept.txt`: unitigs kept.
- `06_chains/`: first chains, plus `chains_self.paf.gz`, `ont_vs_chains.paf.gz`, `depth.windows.tsv`, `depth.unitigs.tsv`.
- `07_chains2/`: round 2, plus `kept2.txt`, `dovetails.gfa`, `links.gfa`, `redundant_chains.txt`, `repeat_chains.txt`, `ont_vs_chains.paf.gz`, `joins.verify.tsv` (read at D03), `restored_arms.tsv` **[C]**.
- `08_final_chains/`: chains after join breaks, plus `breaks.txt`. **These are the chains Hi-C maps to.**

## `09_hic_map/` (stage 05)

`chains.fasta` and `chains.fasta.fai` (symlinks to `08_final_chains/`), `decoys.txt`, `targets.fasta` (chains plus `decoy_*` unitigs) with its bwa index (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`), `bwa_index.log`, `bwa.log`, `samblaster.log`, `HiC.bam`, `HiC.filtered.bam`, `pairs.tsv.gz`.

Only `pairs.tsv.gz` is read later, by any stage or review script. When it was linked in from another build, nothing else is here. The two BAMs are large (about 20 GB together on the genomes tried). `HiC.bam` only saves the remap if stage 05 is ever re-run with other filters, and `HiC.filtered.bam` can be remade from it. Deleting either is the user's call; say what was deleted in `README.md`.

## `10_groups/` (stages 05 to 08)

| File | Stage |
|---|---|
| `chain_breaks.tsv` | 05, read at D04 |
| `alleles.pairs.tsv`, `alleles.blocks.tsv`, `chain_stats.tsv`, `chain_repeat.tsv`, `end_links.tsv`, `shared.candidates.txt` | 06 |
| `groups.tsv`, `groups.chroms.tsv`, `groups.shared.tsv` **[C]**, `hic_group.log`, `shared.txt` **[C]**, `eval_groups.txt` **[ref]**, `chain_forks.txt` **[C]** | 07 |
| `final.tsv`, `final.chroms.tsv`, `final.shared.tsv` **[C]** | 08, the groups after `edits/group_edits.tsv` |

Groups are `g01`, `g02`, ... and chromosomes `c1`, `c2`, ..., both numbered by size (T19).

## `11_ont_chains/`, `12_broken/` (stage 06)

- `12_broken/`: chains after D04 breaks, with restored arms and connectors. `breaks.txt`, `chains.agp`, `chains.fasta`, `chains.fasta.fai`, `pairs.tsv.gz` and `links.tsv` (Hi-C lifted to these chains, no remapping), `connectors.tsv` and `connector_chains.tsv` (unless `restore_connectors: false`), `restored_chains.tsv` **[C]**.
- `11_ont_chains/`: `ont_vs_chains.paf.gz` (reads to `12_broken/chains.fasta`), `depth.windows.tsv`, `depth.unitigs.tsv`, `chain_class.tsv`.

**From here on, "chains" means `12_broken/chains.*`.**

## `12_shared/` (stages 08 to 10)

The working set for ordering. When **[C]**, shared chains appear here as one copy per haplotype, named `<chain>.<group>`, and every file is rewritten for them. When **[N]**, the first seven files are symlinks or copies of those in `12_broken/`, `10_groups/` and `16_validate/`.

| Files | Stage |
|---|---|
| `groups.tsv`, `chains.agp`, `chains.fasta.fai`, `pairs.tsv.gz`, `end_links.tsv`, `chain_stats.tsv`, `telomeres.tsv`, `copies.log` **[C]** | 08 |
| `graph.bridges.tsv`, `graph.tips.tsv` (read at D07), `base.bridges.tsv` **[C]** | 08 |
| `order_links.base.tsv`, `order_links.tsv`, `by_chrom/cN.groups.tsv`, `by_chrom/cN.order.tsv`, `by_chrom/cN.log`, `order.raw.tsv`, `order.tsv`, `groups.edited.tsv` **[edit: order_edits]**, `eval_order.txt` **[ref]**, `unhonoured.txt` | 09 |
| `choices.tsv`, `arm_report.tsv` `arm_depth.tsv` `arm_usage.tsv` **[C]**, `build_args.txt` | 10 |
| `pre/S.agp`, `pre/S.fasta`, `pre/S.junctions.tsv`, `pre/S.placement.tsv`, `pre/S.unplaced.fasta`, `pre/target.fasta` | 10, a trial build to find tails |
| `tails.tsv`, `tail_report.tsv`, `tail_checks.tsv` (when there are candidates) | 10, read at D16 |

## `13_scaffold/` (stages 08, 10, 11)

- `bridges/`: `choices.tsv`, `ids.txt`, `targets.fasta`, `reads.paf.gz`, `mm2.log`, `seqkit.log` (stage 08).
- `names.tsv`, `repeat_chains.fasta` (stage 10).
- `S.graph.agp`, `S.graph.fasta`, `S.graph.fasta.fai`, `S.graph.junctions.tsv`, `S.graph.placement.tsv`, `S.graph.unplaced.fasta`, `target.fasta` (scaffolds plus unplaced), `target.fasta.fai`, `build.log` (stage 11).

Scaffolds are named `chrN_hK`. Unplaced sequence is named `unplaced_<chain>` or `repeat_*`.

## `14_readfill/`, `15_final/` (stage 12)

A fill writes six files from one prefix: `.agp`, `.fills.fasta`, `.report.tsv`, `.links.tsv`, `.membership.tsv`, `.excised.tsv`.

**[N]** one fill over all scaffolds:

- `ont_vs_target.paf.gz`, `S.graph.*` (the six files), `fill.log`, `fills.fasta`, `gap_tails.tsv`.

**[C]** one fill per haplotype, on masked targets:

- `target.masked.tsv`, `target.hK.fasta`, `target.hK.fasta.fai`, `ont_vs_target.hK.paf.gz`.
- `hK/in.agp`, `hK/in.fasta`, `hK/in.fasta.fai`, `hK/names.txt`, `hK/ont.paf.gz`, `hK/fill.hK.*` (the six files), `hK/fill.log`.
- `retry_<scaffold>/` **[edit: fill_retry]**, with `in.agp`, `fill.*`, `fill.log`.
- `merged.agp`, `fills.fasta`, `gap_tails.tsv`.

`15_final/`: `S.final.agp`, `S.final.fasta`, `S.final.fasta.fai`, `S.final.unplaced.fasta`, `S.final.unplaced.fasta.fai`, `stats.txt`. **This is the assembly when there are no raw reads.**

## `17_raw_reads/`, `15_final_raw/` (stages 13 and 14) **[raw]**

`17_raw_reads/`:

- `reads.tsv.gz` (name, length, quality of every raw read), `corrected.lens.tsv.gz`.
- `raw_vs_target.paf.gz` **[N]**, or `raw_vs_target.hK.paf.gz` and `raw.hK.paf.gz` **[C]**.
- `cmp.by_length.tsv`, `cmp.by_q.tsv`, `cmp.sites.tsv`, `cmp.spanners.tsv`, `cmp.taper.tsv`, `compare.log`.
- `neargap.ids`, `neargap.fasta`: the only raw reads written to disk.
- `fill/in.fasta`, `fill/in.fasta.fai`, `fill/candidate.agp`, `fill/candidate.fasta`, `fill/candidate.fasta.fai`, `fill/fills.fasta`.
- **[N]** `fill/raw.paf.gz`, `fill/fill.*` (the six files), `fill/fill.log`, `fill/audit.paf.gz`.
- **[C]** `fill/target.masked.tsv`, `fill/raw.hK.paf.gz`, `fill/hK/` (as in `14_readfill/hK/`, with `raw.paf.gz` and `fill.*`), `fill/audit.masked.tsv`, `fill/audit.hK.paf.gz`.
- `audit.tsv` (read at D18), `gap_tails.tsv`, `phase_reads.tsv` **[C]**.

`15_final_raw/`: the same six files as `15_final/`. **This is the assembly when there are raw reads.** `15_final/` stays as it was.

## `16_validate/` (stage 15, and one file from stage 06)

- `telomeres.chains.tsv` (stage 06), `telomeres.final.tsv`, `stats.txt`.
- `identity.agp`, `pairs.scaffold.tsv.gz`, `pairs.specific.tsv.gz`, `lift.log`, `homolog_groups.tsv`, `phase_scan.txt`, `unplaced_components.txt`, `join_reads.tsv`.
- `dosage.txt`, `sv_candidates.tsv` **[ref]**.
- `tail_joins.tsv`, `tail_joins.log` **[edit: accepted_tails]**.
- `maps/specific/` and `maps/scaffold/`, each with `crosstalk.tsv` and `chrN.png`; `maps/specific.log`, `maps/scaffold.log` **[opt: plot_python]**.
- `reads/`: **[N]** `ont.paf.gz`; **[C]** `target.masked.tsv`, `mask.log`, `all.hK.paf.gz`, `ont.hK.paf.gz`.
- `depth_flat.txt` **[N]**; `depth.windows.tsv`, `depth_classes.log` **[C]**.
- `sv/inv/<scaffold>_vs_<scaffold>.paf` (each haplotype against the first), `sv/inv/merged.tsv`, `sv/inversion_test.txt` and `sv/bp_reads.tsv` (when there are inverted blocks), `sv/final_vs_ref.paf` **[ref]**.
- `sv/sv_calls.tsv`, `curation_targets.tsv`: **yours to fill.**
- `juicebox/`: `out_JBAT.hic`, `out_JBAT.assembly`, `out_JBAT.sizes`, `out_JBAT.components.tsv`, `juicer_tools.log`, `done.ok` **[opt: juicer_tools]**.
- `compleasm/`: `hK/` and `unplaced/` (compleasm's own output, with `summary.txt`), `hK.log`, `unplaced.log`, `complement.txt` **[opt: compleasm]**.

## Re-running

**Within a work directory, use `agentssembler reset WORK STAGE`.** It moves that stage's outputs, and those of every later stage, to `archive/<YYYYmmdd_HHMMSS>/` under their own paths, with the decisions they informed in `archive/<...>/decisions/`. Costly mappings stay in place unless `--remap`. Then `run` rebuilds under the usual names.

- **Do not make `12_shared_v2/`, `15_final_v3/` or the like by hand.** Stages read and write fixed paths, so a suffixed folder is invisible to them, and a reader cannot tell which one is current.
- **`archive/` is the record of what a re-run replaced.** Compare groups by chain against it after any reset above stage 07 (T19). It can be large. Delete it only when the user agrees, and write in `README.md` what was reset, why, and that the archive was deleted.
- **A reset archives both hand tables**, `16_validate/sv/sv_calls.tsv` and `16_validate/curation_targets.tsv`, and stage 15 starts new ones with the header only. Its report names the archived copies. **Treat them as a list to re-check, not rows to restore:** a re-run can change coordinates and haplotype numbers ([LESSONS](LESSONS.md)). Only when nothing before stage 15 was reset are the final files the same, and a row may be carried over once you have confirmed that.
- **To try something that changes the method**, such as a setting that alters the chains, make a second work directory beside the first with `init`, and link the reusable mappings in before `run` ([PROCESS.md](PROCESS.md), "Reusing mappings"). Name it for what differs, and say in both `README.md` files which is current.
- **To keep the output of a `run` command**, redirect it into `logs/`, for example `logs/run.06.out`. The driver does not do this.

## Not in the work directory

Anything you compute to reach a decision, with a script from `scripts/review/`, `scripts/validate/` or your own. Print it, quote the numbers in the decision's reason, and do not leave files among the stage outputs. If a result must be kept, put it under `notes/` at the top level and list it in `README.md`.
