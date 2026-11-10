# Running a mapping elsewhere

The costly steps are mappings. Any of them can be made on another machine, such as an HPC, and put in place: a step is skipped when its output exists. This is how stage 02 was run for all three genomes.

**Nothing here has to be done by the user.** Write the job script, hand it over or submit it if you have access, check what comes back, and carry on.

## Which mappings, and what to bring back

| Stage | Bring back | Needs on the other machine | Hours |
|---|---|---|---|
| `02_map_unitigs` | `02_ont_map/ont_vs_unitigs.paf.gz`, `03_unitig_align/unitigs_self.paf.gz`, `03_unitig_align/unitigs_vs_ref.paf.gz` (with a reference) | Unitig FASTA, corrected reads, reference | 2-6; up to 50 GB |
| `04_chains` | `06_chains/chains_self.paf.gz`, `06_chains/ont_vs_chains.paf.gz`, `07_chains2/ont_vs_chains.paf.gz` | Chains made earlier in the same stage, so only by running the stage there | 1-3 |
| `05_hic_chains` | `09_hic_map/pairs.tsv.gz` | `08_final_chains/`, `07_chains2/repeat_chains.txt`, unitig FASTA, Hi-C reads; and the tools bwa, samblaster, samtools and HapHiC's `filter_bam` | about 1.6 |
| `06_evidence` | `11_ont_chains/ont_vs_chains.paf.gz` | `12_broken/chains.fasta`, made earlier in the same stage, and corrected reads; so only by running the stage there, which needs most of the work directory | about 0.7 |

Stage 02 depends only on the inputs, so it is the one worth moving, and it can be started before anything else. The others depend on results of the run itself, and the stage 06 mapping is short enough that moving it rarely pays.

## Method

The stage scripts read everything from `AS_*` variables, so a stage runs anywhere those are set.

1. **On the other machine,** clone this repository, and make a directory to stand for the work directory. For stage 02 it needs only `00_inputs/` with the files under their fixed names (links are fine):

   ```text
   WORK/00_inputs/unitigs.fasta
   WORK/00_inputs/ont_corrected_reads.fasta
   WORK/00_inputs/reference.fasta        # only with a reference
   ```

   For stage 05, copy `08_final_chains/` and `07_chains2/repeat_chains.txt` there too, and add `00_inputs/hic_R1.fq.gz` and `hic_R2.fq.gz`.

2. **Get the environment** from the real work directory:

   ```bash
   agentssembler env WORK > env.sh
   ```

   Copy `env.sh` across and correct every path in it for that machine: `AS_WORK`, `AS_REPO`, `AS_GUARD_PATH`, the inputs, and the tools. Leave the settings alone.

3. **Run the stage script** in a job. For PBS:

   ```bash
   #!/bin/bash
   #PBS -N map_unitigs
   #PBS -l walltime=12:00:00
   #PBS -l mem=64G
   #PBS -l ncpus=24
   set -eo pipefail
   source /path/to/env.sh          # the corrected copy
   export AS_THREADS=24
   bash "$AS_REPO/stages/02_map_unitigs.sh"
   ```

   For SLURM, replace the header with `#SBATCH -c 24 --mem=64G -t 12:00:00`. A job that is killed leaves `.part` files, never a finished-looking output. Submit it again and it resumes at the first missing mapping.

4. **Check before copying back.** Each output must pass `gzip -t`. Record `md5sum` there, and compare after the copy.

5. **Put the files in place** under the paths in the table, creating the directory if needed, then `agentssembler run WORK`. The stage logs `skip, exists` for each one, writes its report, and the run goes on.

**Say in the work directory's `README.md`** which mappings were made elsewhere, with the minimap2 or bwa version used there. A mapping put in place has no `.log` beside it, so this is the only record.

## What must match

- **The same input files.** A mapping to a different unitig FASTA gives no error, only wrong depths.
- **minimap2 2.27 or later**, for `-x lr:hq`.
- **Stage 05 only: the same chains.** Compare `08_final_chains/chains.fasta.fai` on both sides before mapping, and again before using `pairs.tsv.gz` (T16).

## Reusing another build's mappings

Link them in before `run`, by the same paths. [PROCESS.md](PROCESS.md) says which are safe to reuse; `tests/specs/*.yaml` list the files.
