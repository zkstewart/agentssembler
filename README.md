# agentssembler

A staged, decision-gated process for building haplotype-resolved pseudochromosomes from a hifiasm unitig graph, nanopore reads and Hi-C, driven by an agent with judgement.

Most of an assembly is mechanical. Some of it is not: nearly every significant fix found so far began with a result that looked wrong while every command had succeeded. So the hands-off stages are scripted, and the driver stops at each point where evidence must be weighed. It writes a report, points to a guide for that decision, and records the answer with a reason.

**Inputs:** hifiasm unitigs (`p_utg`, FASTA and GFA), corrected nanopore reads, Hi-C reads. **Optional:** raw basecalled BAMs, a related reference.

**Output:** scaffolds named `<chrom>_h<k>`, unitig-level AGP, unplaced sequence, a table of every decision and hand edit, validation including Juicebox files, and a list of SV candidates with their evidence.

**The evidence rule:** graph overlaps outrank reads, reads outrank Hi-C, Hi-C outranks homolog alignment, and the reference ranks last. The reference names chromosomes and raises alarms. It never places sequence. A disagreement with it is reported as a structural variant candidate, not fixed.

## Status

This process comes from three citrus genomes: one triploid and two diploid. The driver reproduces the first two builds ([docs/TESTING.md](docs/TESTING.md)). Busk is the first genome decided by following the guides, from stage 04 onward, with its first three stages run beforehand. Following them changed the process in 16 places ([examples/Busk/what_changed.md](examples/Busk/what_changed.md)), so expect a new genome to find more. WatsonsChoice and 13Q010 were then rebuilt through the driver with connectors, and those rebuilds are now the adopted assemblies (35 and 32 gaps, against 71 and 63). **No genome has yet gone from raw inputs to a final assembly through the driver in one pass, and the Juicebox route back into the AGP has not been exercised.**

## Start

```bash
conda env create -f environment.yml && conda activate agentssembler
export PATH=$PWD/bin:$PATH
cp config/sample.example.yaml sample.yaml   # then fill it in
agentssembler init sample.yaml work/
agentssembler run work/
```

`run` goes until a decision is needed, then stops:

```text
STOP: decision needed: D01
Report: work/reports/03_depth.md
Guide:  docs/decisions/D01_haplotype_depth.md
Record: agentssembler decide work/ D01 --by NAME --reason '...'
```

Read the report and guide, then record your answer and carry on:

```bash
agentssembler decide work/ D01 --by Claude --reason "Main peak 28x, second peak 56x in long unitigs of graph degree 2,2." --set hap_depth=28 --set repeat_depth=123 --set collapsed=true
agentssembler run work/
```

## Commands

| Command | Does |
|---|---|
| `init CONFIG WORK` | Create a work directory |
| `run WORK [--to STAGE] [--only STAGE]` | Run until the next decision |
| `status WORK` | Stages done, decisions open |
| `decide WORK Dnn --by NAME --reason TEXT [--set k=v]` | Record a decision |
| `report WORK [STAGE]` | Print reports |
| `reset WORK STAGE [--remap] [--keepDecision]` | Undo a stage and everything after it, moving outputs to `archive/` and setting aside the decisions they informed. Mappings stay unless `--remap` |
| `env WORK` | Print a stage's environment, to run it elsewhere |

## Where things are

| Path | Contents |
|---|---|
| [AGENTS.md](AGENTS.md) | **Start here if you are an agent.** Rules and the loop |
| [docs/PROCESS.md](docs/PROCESS.md) | Start to end: the whole run, then each stage and what to check |
| [docs/decisions/](docs/decisions/README.md) | One guide per decision |
| [docs/SETTINGS.md](docs/SETTINGS.md) | Every setting and how to choose it |
| [docs/TRAPS.md](docs/TRAPS.md) | Things that cost real time |
| [docs/LESSONS.md](docs/LESSONS.md) | Findings from the three genomes, and my errors |
| [docs/IMPROVING.md](docs/IMPROVING.md) | Ideas not yet tried, and what is out of scope |
| [docs/SCRIPTS.md](docs/SCRIPTS.md) | Every script |
| [docs/TESTING.md](docs/TESTING.md) | Tests and how the regression was run |
| `examples/` | Hand decisions from each genome, as worked examples |
| `stages/` | One bash script per stage |
| `scripts/` | The tools |
| `legacy/` | The earlier reference-guided pipeline. Not part of this process |

## What gets automated, and what does not

| Automated | Needs judgement |
|---|---|
| Mapping, depth, chains, join verification | Haplotype depth and collapse (D01) |
| Hi-C evidence, allele table, end links | Chimera calls (D04) |
| Clustering and ordering searches | Group count, membership, shared chains (D05, D06, D14) |
| Bridge finding and ranking, connector restoration, arm search | Telomeric tips (D07), odd orders (D09, D10), phase (D15, D17) |
| Tail candidates and copy tests | Whether to place any tail (D16) |
| Read filling, join audit | Which read-derived joins stand (D11, D18) |
| Validation, and scans for forks, read tails, unplaced tangles and short homologs | SV claims (D12), and what those scans mean (D06, D11) |
| | Adoption, risk, Juicebox review: the user (D13) |

## Requirements

Python 3.10 or later with NumPy, PyYAML and matplotlib. minimap2 2.27 or later, samtools 1.12 or later, seqkit, bwa, samblaster, HapHiC's `filter_bam`. Optional: juicer_tools and compleasm. About 16 cores and 30 GB RAM were enough beyond the unitig alignments, which ran on an HPC with up to 50 GB.
