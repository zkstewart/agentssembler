# Settings

Everything goes in the config (`config/sample.example.yaml`), except the values set by D01. These are the values that worked on the two genomes tried: treat them as starting points.

## From the sample

| Setting | 13Q010 | WatsonsChoice | How to set it |
|---|---|---|---|
| `chromosomes`, `ploidy` | 9, 3 | 9, 2 | Known biology. Ask the user if unsure |
| `hap_depth` | 27 | 28 | Main histogram peak (D01) |
| `repeat_depth` | 120 | 123 | About 4.4 x `hap_depth` (D01) |
| `collapsed` | false | true | Peak near 2x in long unitigs (D01) |
| `hap_first` | false | true | Defaults to true at ploidy 2. Pairing of groups into chromosomes is implemented for ploidy 2 only |
| `telomere_motif` | `CCCTAAA` | `CCCTAAA` | Plant default; `CCCTAA` is the vertebrate type. Confirm for other taxa |

## Defaults

| Setting | Value | Meaning |
|---|---|---|
| `min_within` | 700 | Hi-C pairs/Mb within a chain below which it carries no placement information. Look at `10_groups/chain_stats.tsv` for a low-density class |
| `link_bonus` | 20000 | Order score for honouring a link. Set above the Hi-C differences between orders, so a confirmed link is effectively a hard constraint. This follows from the evidence order, not from tuning: a collinearity audit by bonus (none 88.7%, 200 89.9%, 2000 91.7%, 20000 92.7%) only confirmed it |
| `hops` | 8 | Longest bridge path. 3 hops reversed a block on WatsonsChoice |

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
| `compute.threads` | Threads per job. About 16 cores and 30 GB RAM were enough for both genomes, excluding the HPC mappings of stage 02, which used 12 cores and up to 50 GB |
| `compute.min_free_gb` | Mappings stop below this. Default 40 |
| `compute.guard_path` | Path whose free space counts. **Under WSL use the host drive**, such as `/mnt/d`, because `df /` overstates it |
