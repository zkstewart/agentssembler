# agentssembler

`agentssembler` provides a stable set of tools and information for an AI agent to autonomously assemble a haplotype-phased genome. Its aim is to work on any organism regardless of ploidy or biological complexity, although this has currently only been tested on diploid and triploid citrus plants. The reuse of tested tools minimises the risk of mistakes being made by a lone agent "reinventing the wheel". And the decision making guides ensure those tools get used in a consistent and repeatable fashion.

# Table of Contents
- [Getting started (for humans)](#getting-started-for-humans)
- [Getting started (for agents)](#getting-started-for-agents)
- [Background](#background)
  - [Why does `agentssembler` exist?](#why-does-agentssembler-exist)
  - [How is `agentssembler` implemented?](#how-is-agentssembler-implemented)
  - [What is `agentssembler` not meant to be?](#what-is-agentssembler-not-meant-to-be)
- [An agent's instructions for running `agentssembler`](#an-agents-instructions-for-running-agentssembler)
  - [Commands](#commands)
  - [Where things are](#where-things-are)
  - [What gets automated, and what does not](#what-gets-automated-and-what-does-not)
  - [Requirements](#requirements)
- [How to cite](#how-to-cite)

# Getting started (for humans)

Set up a computational environment for running `agentssembler`:
```bash
git clone https://github.com/zkstewart/agentssembler.git
cd agentssembler

conda env create -f environment.yml
conda activate agentssembler
export PATH=$PWD/bin:$PATH
```

Set up a config file or prompt an agent to do this for you:
```text
cp config/sample.example.yaml sample.yaml   # then fill it in

< --- OR --- >

Natural language prompt telling the agent where your files are, 
what your organism is, what the expected num. chromosomes and 
ploidy is, and for it to locate and edit the necessary program 
paths. You may also ask it to find the telomere motif relevant 
to your organism.
```

This config file will allow you to specify the data inputs including:
- A hifiasm unitig graph (`*.p_utg.gfa`) as well as the extracted FASTA from that file.
- Long corrected nanopore reads (FASTA format as output by `dorado correct` or equivalent).
- Hi-C reads (FASTQ with/without `gz` compression; trimmed or untrimmed).
- (Optional) raw nanopore data (BAM, CRAM or SAM; or FASTQ or FASTA, gzipped or not).
- (Optional) a related reference genome (single haplotype FASTA format).

In addition to biological information including:
- Number of chromosomes
- Ploidy number
- Telomere motif

Then provide a prompt to start the agent on the assembly process. An example is given below with `<your answer>` segments for you to fill your personalised information:

```text
Within the directory at `<my working directory>` you will find files to be used as input to the `agentssembler` genome assembly pipeline. Specifically, the input files can be found at:

- `<location of hifiasm utg gfa and fasta>` contains the unitigs GFA as output by hifiasm, as well as a FASTA extracted from that GFA.
- `<location of hi-c data>` contains the paired Hi-C reads.
- `<location of nanopore data>` contains corrected nanopore reads in FASTA format <OPTIONAL: and raw reads (BAM, CRAM, SAM, FASTQ or FASTA)>.
- `<location of reference genome>` which is a related reference genome.

The goal of this project is to assemble a genome of <my organism>. It is expected to be a <ploidy level> and have <number> chromosomes.

The `agentssembler` GitHub repository is located at `<cloned location>`. Contained within are Python scripts for performing various operations such as read alignment to the hifiasm unitigs. There is also an extensive amount of documentation detailing the process and pipeline to be followed; the `README.md` and `AGENTS.md` files are intended to be starting points for you to explore the repository and understand the intended methodology to be applied.

After exploring the repository, please begin by initialising a directory at `<my intended work directory>` with a suitable config file and beginning the pipeline. As much as possible, I would like you to adhere to the layout and process described by `agentssembler`. Of course, any unexpected bugs or issues should be dealt with appropriately and according to your best judgement. Otherwise, the pipeline should be followed through to its completion. Each decision point in the pipeline will be an opportunity to evaluate whether a deviation from the standard `agentssembler` process could result in a better assembly outcome. You can make a note of these ideas, and raise them as suggestions at the end of the pipeline if applicable. Otherwise, the goal is to adhere to the established pipeline as much as possible and follow it through until the assembly output has been generated.
```

# Getting started (for agents)
Before beginning an analysis, read [AGENTS.md](AGENTS.md), [docs/TRAPS.md](docs/TRAPS.md), [docs/decisions/README.md](docs/decisions/README.md) and "The whole run" in [docs/PROCESS.md](docs/PROCESS.md). Read the rest of `docs/` when the table below says to: each stage's section of PROCESS and each decision's guide when the run reaches them, and the reference documents when you need them. The development record (LESSONS, RECORD, IMPROVING, TESTING) is not needed to run an assembly.

Coordinate with the human user to set up the necessary conda environment and/or the `sample.yaml` configuration file.

Then start the pipeline, writing results to e.g., the `work/` directory:

```bash
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

Read the 'Report' markdown file and consult the indicated 'Guide' markdown (in this repository) to obtain relevant information for making a decision at this step. Then record your answer and carry on (example below):

```bash
agentssembler decide work/ D01 --by Claude --reason "Main peak 28x, second peak 56x in long unitigs of graph degree 2,2." --set hap_depth=28 --set repeat_depth=123 --set collapsed=true
agentssembler run work/
```

# Background
## Why does `agentssembler` exist?
Haplotype-phased genome assembly using Hi-C data is a process that typically requires a significant amount of manual intervention (e.g., using Juicebox) to produce the highest quality output. The need for manual intervention stems from the difficulty of crafting a software pipeline that can automatically adjust to *every possible data input*. The best pipelines make use of mathematical formulae and/or heuristic rules of thumb for assessing input data and determining numerical cut-offs to control **decisions** such as:

- Which haplotype do a set of contigs belong to?
- What order should those contigs be joined in to produce an accurate telomere-to-telomere chromosomal scaffold?
- If two contigs overlap end-to-end, is that a trustworthy join, or is a repeat region being collapsed?
- If two contigs do not overlap, does the Hi-C data unambiguously anchor them with a gap between them? Or can the raw long reads assist in supporting their union?

Tuning the algorithms and heuristic rules to hit the goldilocks zone of working *just right* for every organism (plant, animal, prokaryotic, and so on) and genetic background (simple or repeat-rich, low or high heterozygosity, diploid or polyploid) is nigh impossible. Software seeking to achieve this goal necessarily must become huge in scope (with hundreds of thousands of lines of code) to handle all edge cases, or invent complicated algorithms which attempt to make some sense out of the messiness of real biological data.

Almost all the difficulty of producing a good genome comes from making the right **decisions** of which contigs to join and in what order. Automatic methods rely on evidence from read alignments that must ultimately be transformed into some sort of numerical value that can be interpreted to make a choice: to join or not to join? An AI agent can answer that question through their "best judgement" by flexibly interpreting the data, without reliance upon hard-coded rules baked into an unchanging software tool.

## How is `agentssembler` implemented?
Assembly proceeds through stages, each involving the alternation between:
- Data generation, then
- AI review of stage output to weigh the evidence and determine appropriate numerical values to use in the next stage

Review may involve viewing a histogram to identify peak depth, or a Hi-C contact map to identify haplotype 
mosaicism. These and other processes are things that a human would normally perform, requiring time and expertise.

Moreover, it is rare that a genome is assembled perfectly in a single take. More often, some iteration is required. For example, a suitable numerical cut-off may be decided based on the best possible evidence, but the result may still *look wrong* despite everything working correctly. Iteration is then required to either: find a more suitable number, or find a different program or process better suited.

AI agents are now capable of identifying results that *look wrong* despite all else being correct. That ability enables a new approach to genome assembly.

## What is `agentssembler` *not* meant to be?

`agentssembler` is a proof of concept to answer the question "can we replace human decision making in the manually-intensive work of genome assembly and curation?" Obtaining a high-quality and useable genome is a benefit of this, but it is not the main point.

This software is *not* the result of painstaking coding effort; it is entirely vibe coded around a simple premise of *can an AI agent do this?* A human has written this README and, beyond that, everything has been prompted into existence.

Thus, `agentssembler` is *not* meant to be a tool that is used for decades to come. It *will* be deprecated by the work of a more advanced agent; at the current rate of improvement this will, in all likelihood, happen very soon.

# An agent's instructions for running `agentssembler`

## Commands

| Command | Does |
|---|---|
| `init CONFIG WORK` | Create a work directory |
| `run WORK [--to STAGE] [--only STAGE]` | Run until the next decision |
| `status WORK` | Stages done, decisions open |
| `decide WORK Dnn --by NAME --reason TEXT [--set k=v]` | Record a decision |
| `report WORK [STAGE]` | Print reports |
| `reset WORK STAGE [--remap] [--keepDecision]` | Undo a stage and everything after it, moving outputs to `WORK/archive/<timestamp>/` and setting aside the decisions they informed. Mappings stay unless `--remap`, or unless their target is rebuilt: a reset of `11_build` or earlier also drops stages 12-13's mappings to the scaffolds |
| `env WORK` | Print a stage's environment, to run it elsewhere |
| `compare WORK OTHER` | Side-by-side metrics of two finished builds, for choosing one at D13 |

## Where things are

| Path | Contents | Read when |
|---|---|---|
| [AGENTS.md](AGENTS.md) | **Start here if you are an agent.** Rules and the loop | Before starting |
| [docs/TRAPS.md](docs/TRAPS.md) | Things that cost real time | Before starting |
| [docs/decisions/](docs/decisions/README.md) | One guide per decision | The index before starting; each guide when the run stops at it |
| [docs/PROCESS.md](docs/PROCESS.md) | Start to end: the whole run, then each stage and what to check | "The whole run" before starting; each stage's section when the run reaches it |
| [docs/LAYOUT.md](docs/LAYOUT.md) | The work directory: which stage writes what, every file, and the files you write | When looking for a file |
| [docs/SETTINGS.md](docs/SETTINGS.md) | Every setting and how to choose it | When writing the config, or when a guide sets a value |
| [docs/HPC.md](docs/HPC.md) | Running a mapping on another machine and putting it in place | Only if mapping elsewhere |
| [docs/TANGLES.md](docs/TANGLES.md) | D16's tangle check: lost-route forks, and TTT for small tangles | At D16 |
| [docs/SCRIPTS.md](docs/SCRIPTS.md) | Every script | When a guide names a script |
| [docs/LESSONS.md](docs/LESSONS.md) | Principles the genomes taught, each linked to its evidence and its rule | Development record: not needed for a run |
| [docs/RECORD.md](docs/RECORD.md) | The evidence behind the process, genome by genome: experiments, numbers, errors | Development record: not needed for a run |
| [docs/IMPROVING.md](docs/IMPROVING.md) | What comes next: ideas not yet tried, open questions, what is out of scope or not worth repeating | Development record; at D13, before proposing a trial |
| [docs/TESTING.md](docs/TESTING.md) | Tests and how the regression was run | When changing code |
| `examples/` | Hand decisions from each genome, as worked examples | At a decision, for the form expected |
| `stages/` | One bash script per stage | |
| `scripts/` | The tools | |

## What gets automated, and what does not

| Automated | Needs judgement |
|---|---|
| Mapping, depth, chains, join verification | Haplotype depth and collapse (D01) |
| Hi-C evidence, allele table, end links | Chimera calls (D04) |
| Clustering and ordering searches | Group count, membership, shared chains (D05, D06, D14) |
| Bridge finding and ranking, connector restoration, arm search | Telomeric tips (D07), odd orders (D09, D10), phase (D15, D17) |
| Tail candidates and copy tests, lost-route fork scan | Tangle check and whether to place any tail (D16) |
| Read filling, join audit | Which read-derived joins stand (D11, D18) |
| Validation, and scans for forks, read tails, unplaced tangles and short homologs | SV claims (D12), and what those scans mean (D06, D11) |
| | Adoption, risk, Juicebox review: the user (D13) |

## Requirements

Python 3.10 or later with NumPy, PyYAML and matplotlib. minimap2 2.27 or later, samtools 1.12 or later, seqkit, bwa, samblaster. Optional: juicer_tools and compleasm; TTT (with PuLP 2.9 and GLPK) and GraphAligner, each in its own environment, if D16's tangle check calls for them ([TANGLES](docs/TANGLES.md)). About 16 cores and 30 GB RAM were enough beyond the unitig alignments, which ran on an HPC with up to 50 GB.

# How to cite
`agentssembler` is currently being evaluated. It is the author's hope that this will validate its use, and that reviewers will see its value as a potential publication.
