# Settings

Everything goes in the config (`config/sample.example.yaml`), except the values set by D01. The defaults are the values that worked on the genomes tried: treat them as starting points.

## From the sample

| Setting | How to set it |
|---|---|
| `chromosomes`, `ploidy` | Known biology. Ask the user if unsure |
| `hap_depth` | Main histogram peak ([D01](decisions/D01_haplotype_depth.md), which lists the values used so far) |
| `repeat_depth` | About 4.4 x `hap_depth` (D01) |
| `collapsed` | Long unitigs near 2x, with or without a visible peak (D01) |
| `hap_first` | Defaults to true at ploidy 2. Pairing of groups into chromosomes is implemented for ploidy 2 only |
| `telomere_motif` | `CCCTAAA` is the plant default; `CCCTAA` is the vertebrate type. Confirm for other taxa |

## Defaults

| Setting | Value | Meaning |
|---|---|---|
| `min_within` | 700 | Hi-C pairs/Mb within a chain below which it carries no placement information. Look at `10_groups/chain_stats.tsv` for a low-density class |
| `link_bonus` | 20000 | Order score for honouring a link. Set above the Hi-C differences between orders, so a confirmed link is effectively a hard constraint. This follows from the evidence order, not from tuning |
| `restore_connectors` | true | Restore dropped unitigs that are the graph path between chain ends (stage 06) |
| `hops` | 12 | Longest bridge path. Connectors add hops without much sequence: 3 and 8 have each missed a real path, and 16 found the same as 12 |

## Derived from `hap_depth`

Set automatically, so the scripts have no sample-specific defaults. Rounded to whole numbers.

| Value | Multiple | Used by |
|---|---|---|
| 1x lower bound | 0.64 | `chain_class.py` |
| Shared range | 1.54 to 2.68 | `chain_class.py`, `restore_arms.py` |
| 2x threshold | 1.5 | `assign_arms2.py`, `make_tails.py` |
| Highest tail depth | 3 | `make_tails.py` |
| Real allele | 0.36 | `assign_arms2.py` |
| Largest allele pair sum | 2.68 | `assign_arms2.py` |

## Compute

| Setting | Meaning |
|---|---|
| `compute.threads` | Threads per job. About 16 cores and 30 GB RAM were enough, except for the stage 02 mappings, which took up to 50 GB |
| `compute.min_free_gb` | Mappings stop below this. Default 40 |
| `compute.guard_path` | Path whose free space counts. **Under WSL use the host drive**, such as `/mnt/d`, because `df /` overstates it |
