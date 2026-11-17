# Settings

Everything goes in the config (`config/sample.example.yaml`), except the values set by D01. The defaults are the values that worked on the genomes tried: treat them as starting points.

## From the sample

| Setting | How to set it |
|---|---|
| `chromosomes`, `ploidy` | Known biology, but see "What to enter for a polyploid or an inbred sample" below. Ask the user if unsure |
| `hap_depth` | Main histogram peak ([D01](decisions/D01_haplotype_depth.md), which lists the values used so far) |
| `repeat_depth` | Above anything the haplotypes can share: the larger of 4.4 and `ploidy` + 1.4, times `hap_depth` (D01). That is 4.4 x up to ploidy 3 and 5.4 x at ploidy 4 |
| `collapsed` | Long unitigs near 2x, with or without a visible peak (D01) |
| `hap_first` | Defaults to true at ploidy 2, and false above it. Can be set true at any ploidy ([D05](decisions/D05_group_count.md)), but only ploidy 2 has been run |
| `length_scale` | Left out, it is measured after stage 04 (below). Set a number to fix it: `1` gives the windows the citrus genomes were built with |
| `telomere_motif` | `CCCTAAA` is the plant default; `CCCTAA` is the vertebrate type. Confirm for other taxa |

## What to enter for a polyploid or an inbred sample

`ploidy` is the number of haplotypes that are **alleles of each other and distinct in the unitigs**, and `chromosomes` is the number of sets of them. That is not always the textbook ploidy.

| Sample | `chromosomes` | `ploidy` | Why |
|---|---|---|---|
| Diploid or triploid, heterozygous | Haploid number | 2 or 3 | The tested cases |
| Allopolyploid, heterozygous (subgenomes pair only within themselves) | Haploid number of the whole genome: 26 for cotton, not 13 | 2 | Homeologs are different chromosomes. They share no Hi-C, so entering 13 and 4 asks grouping to merge chromosomes that never touch |
| Allopolyploid, inbred | As above | 1 | Each subgenome's two copies are identical and collapsed |
| Diploid, inbred or doubled haploid | Haploid number | 1 | The main depth peak is both copies together. Entering 2 makes almost everything shared |
| Autopolyploid, heterozygous | Haploid number of one set | 4, 6, ... | Set `collapsed: true` at D01 if any long unitigs sit at 2x or more. Shared chains then get a dosage ([D14](decisions/D14_shared_chains.md)) |
| Autopolyploid from a recent doubling (two pairs of identical haplotypes) | Haploid number of one set | 2 | It is a diploid at twice the depth |

**None of the rows below the first has been run on real data.** `ploidy: 1` and ploidies above 3 are exercised only by unit tests on made-up data. Read the histogram at D01 before trusting the row you chose: the number of depth peaks in long unitigs says how many haplotypes the unitigs separate. If the peaks disagree with the config, tell the user, change `ploidy` in `WORK/config.yaml` and `reset 04_chains`.

An allopolyploid's subgenomes are not labelled. Chromosome names come from the reference, if it holds both subgenomes.

## Defaults

| Setting | Value | Meaning |
|---|---|---|
| `min_within` | 700 | Hi-C pairs/Mb within a chain below which it carries no placement information. Look at `10_groups/chain_stats.tsv` for a low-density class |
| `link_bonus` | 20000 | Order score for honouring a link. Set above the Hi-C differences between orders, so a confirmed link is effectively a hard constraint. This follows from the evidence order, not from tuning |
| `restore_connectors` | true | Restore dropped unitigs that are the graph path between chain ends (stage 06) |
| `hops` | 12 | Longest bridge path. Connectors add hops without much sequence: 3 and 8 have each missed a real path, and 16 found the same as 12 |

## Derived from `hap_depth`

Set automatically, so the scripts have no sample-specific defaults. Rounded to whole numbers. Bounds on shared sequence rise with `ploidy`, since up to that many haplotypes can carry one collapsed chain. The second column is the value at ploidy 2, which the diploid genomes were built with. The triploid was built with those values too, before the bounds followed ploidy: its tail candidates may differ on a rebuild.

| Value | Multiple at ploidy 2 | In general | Used by |
|---|---|---|---|
| 1x lower bound | 0.64 | 0.64 | `chain_class.py` |
| Shared range | 1.54 to 2.68 | 1.54 to `ploidy` + 0.68 | `chain_class.py`, `restore_arms.py` |
| 2x threshold | 1.5 | 1.5 | `assign_arms2.py`, `make_tails.py` |
| Copy target of an arm or tail chain | 1 or 2 | 1 below the 2x threshold, else depth over `hap_depth` rounded, from 2 to `ploidy` | `assign_arms2.py`, `make_tails.py` |
| Highest tail depth | 3 | `ploidy` + 1 | `make_tails.py` |
| Real allele | 0.36 | 0.36 | `assign_arms2.py` |
| Largest allele pair sum | 2.68 | 2.68, plus 1 for each copy of the target above 2 | `assign_arms2.py` |

## Derived from the length scale

Windows that are a fraction of a chromosome cannot be the same bp on yeast, citrus and wheat. They are multiplied by the **length scale: mean chromosome length over 40 Mb**, where mean chromosome length is the bp of the stage 04 chains over `chromosomes` x `ploidy`. Citrus is about 1 (13Q010: 38 Mb, so 0.95), and a diploid citrus build moves by under a fifth from the fixed values used before.

- **Held down by chain N50,** so that seeds are never longer than half the chain bases can reach: a 14 Gb genome with 2 Mb chains gets 2, not 17. The stage 04 report says when this applies.
- **Never below 0.02.**
- **Measured once,** when stage 04 finishes, and kept in `state/length_scale.json` for the rest of the build. `reset` to stage 04 or earlier measures it again.
- **Underestimated when most of the genome is collapsed,** since shared sequence is counted once: by up to half at ploidy 2.
- **Untested away from 1.** The stage 04 report raises an alarm below 0.67 or above 1.5.

| Value | At scale 1 | Used by |
|---|---|---|
| Seed length for grouping | 1 Mb | `hic_group.py --seedLen` |
| Chimera test: reach each side of a junction, flank for the contact profile, shortest side tested | 500 kb, 1 Mb, 100 kb | `chain_breaks.py` |
| Contact map bin | 250 kb | `contact_map.py` |
| Dosage bin against the reference; runs reported at 4 bins | 500 kb | `dosage.py` |

**Not scaled:** windows set by read length (10 kb depth windows, 5 kb overlaps, 20 kb allele blocks, the 50 kb allele conflict), counts of Hi-C pairs (`min_within`, link and pair minimums), and the fixed windows inside the review and validation scripts (`phase_scan.py`, `sv_candidates.py`, `inversion_test.py`). On a genome far from citrus scale, give those scripts their own values.

## Compute

| Setting | Meaning |
|---|---|
| `compute.threads` | Threads per job. About 16 cores and 30 GB RAM were enough, except for the stage 02 mappings, which took up to 50 GB |
| `compute.min_free_gb` | Mappings stop below this. Default 40 |
| `compute.guard_path` | Path whose free space counts. **Under WSL use the host drive**, such as `/mnt/d`, because `df /` overstates it |
