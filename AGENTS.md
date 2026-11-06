# Instructions for an agent running an assembly

You are building haplotype-resolved pseudochromosomes from a hifiasm unitig graph, corrected nanopore reads and Hi-C. The optional inputs are raw basecalled reads and a related reference. Most stages run themselves. Your job is the judgement at each stop, and honest reporting.

## Loop

```bash
agentssembler init sample.yaml WORK    # once
agentssembler run WORK                 # runs until a decision is needed, or a stage fails
# read reports/STAGE.md and the named guide; inspect evidence; make edits; then
agentssembler decide WORK D04 --by NAME --reason "..."
agentssembler run WORK                 # continue
```

`agentssembler status WORK` shows where things stand. See [docs/PROCESS.md](docs/PROCESS.md) for each stage and [docs/decisions/](docs/decisions/README.md) for each fork.

## Rules that decide the forks

**Evidence order.** When sources disagree, the higher one wins:

1. Graph overlaps
2. Reads spanning a junction with unique anchors
3. Hi-C
4. Alignment to homologous haplotypes
5. The reference: names chromosomes and audits the finished result only

**Working rules:**

- The reference may tell you **where to look**, never **what to change**. Every fix needs evidence from levels 1-4.
- A disagreement with the reference is reported as an SV candidate. It is never fixed toward the reference.
- Leave sequence unplaced when the data do not support a position. Unplaced is a result, not a failure.
- Every hand edit goes in an `edits/` file with its reason, and is applied by re-running a stage. Never edit an output by hand.
- A claim about biology must survive a direct test before it is written down. See [D12](docs/decisions/D12_sv_calls.md).
- Do not stop at numbers that look good. Every significant fix on the two genomes tried began with a result that looked wrong while every command had succeeded.

## What only you can do

- Read the histogram, not a formula ([D01](docs/decisions/D01_haplotype_depth.md)).
- Weigh two signals against each other ([D04](docs/decisions/D04_chimeric_chains.md), [D09](docs/decisions/D09_order.md)).
- Refuse tempting joins ([D16](docs/decisions/D16_tails.md), [D18](docs/decisions/D18_raw_read_joins.md)). Every automated label here has been wrong at least once: `only`, `reads_overlap`, read counts.
- Say what you do not know. Arm phase is mostly unknowable on a diploid, and that belongs in the report.

## What belongs to the user

Ask at [D13](docs/decisions/D13_user_decisions.md): which assembly to adopt, how much risk to accept for completeness, and Juicebox review. Also ask when chromosome number or ploidy is uncertain. Do not record a user's decision on their behalf.

## Reporting

- Report failures and skipped steps plainly, and correct your earlier claims openly when they turn out wrong.
- Give counts with their denominators and compare with a control where you can.
- Prefer medians per site to pooled ratios: a pair of repeat arrays once inflated a pooled figure threefold.
- Haplotype numbers mean nothing across chromosomes. Never interpret `h1` across chromosomes.
- Declare that the reference shaped where you looked, if it did.

## Disk and long jobs

A full disk can crash the machine, especially when the work disk is a virtual disk on another drive, as under WSL: `df` inside then overstates free space. Set `compute.guard_path` to the host drive, for example `/mnt/d`. Mappings stop when free space falls below `compute.min_free_gb`. Never write raw reads to an uncompressed FASTA: stage 13 streams them.

Do not use `pkill -f PATTERN` from a command containing that pattern: it kills your own shell. Find the PID with `pgrep`, then `kill` it.

Read [docs/TRAPS.md](docs/TRAPS.md) before starting. It is short, and every entry cost real time.
