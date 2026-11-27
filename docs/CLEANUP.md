# Freeing disk in a finished work directory

A finished work directory holds 20-40 GB that no stage reads again: the Hi-C BAM and its bwa index, compleasm's own intermediates, and FASTAs that are plain concatenations or are written from an AGP. `scripts/cleanup/clean_work.py` deletes only those. It records each file with its checksum, and remakes any of them on demand.

**When:**
- Only after the run has finished (`state/15_validate.done`; the script refuses otherwise), and normally after D13.
- Only when the user asks: deleting is their call, as for the archive ([LAYOUT](LAYOUT.md), "Re-running").

**What it never touches:**
- Inputs and anything under `00_inputs/`.
- Every read mapping (`*.paf.gz` and its `.target`).
- `09_hic_map/pairs.tsv.gz` and `pairs.rescued.tsv.gz`, and all pair tables made from them.
- AGPs, chain FASTAs, the scaffold and final FASTAs, compleasm's tables and `summary.txt`.
- Reports, decisions, edits, `state/` and logs.
- Files the driver did not make.

So `tests/regress.py`, whose seeds are mappings and pairs, needs nothing it deletes ([TESTING](TESTING.md)). Neither does a `reset`: a stage that is re-run remakes its own intermediates.

## Use

Run from the repository (`S`), on the work directory (`WORK`):

```bash
python3 $S/scripts/cleanup/clean_work.py list WORK      # what would go, per rule, and what a rule covers but keeps, with why
python3 $S/scripts/cleanup/clean_work.py check WORK     # as delete, but deletes nothing
python3 $S/scripts/cleanup/clean_work.py delete WORK
python3 $S/scripts/cleanup/clean_work.py restore WORK   # remake what CLEANED.tsv lists and is missing
```

**Options:**
- `--rules a,b` limits the run to the named rules; `--skip a,b` leaves rules out. For example, `restore WORK --rules hic_map` remakes only the Hi-C mapping.
- `--noArchive` leaves `archive/` alone.
- `--tmp DIR` puts scratch files in `DIR` instead of `WORK`.

The tools come from the work directory's config, through `agentssembler env`.

**`delete` in steps:**
1. Every `exact` file (table below) is rebuilt in scratch and must match the original byte for byte. A file that does not match is kept, and the summary names it.
2. Every other file is checksummed. For the Hi-C BAM, the script records checksums of the pairs it yields.
3. The files are deleted.
4. Each deleted file gets a row in `WORK/CLEANED.tsv`: path, bytes, md5, rule, check, date. Deleting a file again after a restore replaces its row.
5. A dated line is added to `WORK/README.md`.

**`restore`** rebuilds each file and compares it with its row in `CLEANED.tsv`. It works in dependency order, so a concatenation is rebuilt before the subsets made from it. It exits 1 if an `exact` file or the Hi-C pairs differ; compleasm files are reported but never fail it. It first applies the stages' disk guard (`compute.guard_path`, `compute.min_free_gb`).

## The files

Paths are relative to `WORK`, and to each `WORK/archive/<time>/` unless the rule is archive only. `S` is the sample name; `hK` is a haplotype number. Every row was checked on Busk (diploid, collapsed) and 13Q010 (triploid). `exact` rows were rebuilt to an identical file; `restore` was run for the cheap rows and for one compleasm run. The Hi-C remap has not been run through `restore`.

| Rule | Files | Made by | Why nothing needs it | Remade by | Check |
|---|---|---|---|---|---|
| `hic_map` | `09_hic_map/HiC.bam`; `09_hic_map/targets.fasta.{amb,ann,bwt,pac,sa}` | stage 05 | Only `pairs.tsv.gz` and `pairs.rescued.tsv.gz` are read later; stage 05 maps again only when `pairs.tsv.gz` is missing | `bwa index`, then `bwa mem -5SP` → `samblaster` → `samtools view -F 3340`, as stage 05 runs them; 1-2 h | The BAM: pairs remade from it must match the recorded pairs. Its header names paths, so it has no md5. Keep `compute.threads` as it was, since bwa's output depends on `-t` |
| `hic_targets` | `09_hic_map/targets.fasta` | stage 05 | As above | `chains.fasta` plus the `decoys.txt` unitigs, named `decoy_*` | exact |
| `compleasm` | `01_genes/<lineage>/{miniprot_output.gff,translated_protein.fasta}`; `16_validate/compleasm/{hK,hK_unloc,unplaced}/<lineage>/` the same two | compleasm, in stages 01 and 15 | Stages read `full_table*.tsv` and `summary.txt`, which stay. A run with no known input is left alone | compleasm on the same input, in scratch; the two files are copied back. About 15 min a run | Recorded md5, reported only: compleasm is not deterministic, so a rerun can differ without anything being wrong. `restore` also compares the rerun's `summary.txt` with the kept one and says if it moved. It never fails on these files. The kept tables are what the stages used |
| `pre_target` | `12_shared/pre/target.fasta` | stage 10 | Read only in stage 10, which remakes it when missing | `cat pre/S.fasta pre/S.unplaced.fasta` | exact |
| `scaffold_target` | `13_scaffold/target.fasta`, **archive only** | stage 11 | Stages 12 and 13 read the live one | `cat S.graph.fasta S.graph.unplaced.fasta` | exact |
| `rawfill_in` | `17_raw_reads/fill/in.fasta` | stage 13 | Read only in stage 13 | `cat 15_final/S.final.fasta 15_final/S.final.unplaced.fasta` | exact |
| `rawfill_hap_in` | `17_raw_reads/fill/hK/in.fasta` (collapsed) | stage 13 | Read only in stage 13 | `seqkit grep -f hK/names.txt fill/in.fasta -w 0` | exact |
| `readfill_hap_in` | `14_readfill/hK/in.fasta` (collapsed), **archive only** | stage 12 | Stage 13 reads the live one | `seqkit grep -f hK/names.txt 13_scaffold/target.fasta -w 0` | exact |
| `readfill_masked` | `14_readfill/target.hK.fasta` (collapsed), **archive only** | stage 12 | Stage 13 maps raw reads to the live one | `mask_shared.py S.graph.agp 13_scaffold/target.fasta` | exact |
| `candidate` | `17_raw_reads/fill/candidate.fasta` | stage 13 | Stage 14 reads `candidate.agp`, not this | `agp_to_fasta.py candidate.agp` with the unitigs and both `fills.fasta` | exact |
| `bridge_reads` | `13_scaffold/bridges/reads.fasta` | stage 08 | Stage 08 deletes it itself; only archived or older runs have it | `seqkit grep -f ids.txt` on the corrected reads, then `fq2fa`; slow | exact |

`.fai` files stay: they are small, and some stages read them.

**A file stays, and `list` says why, when:**
- something it is remade from is missing, for example the Hi-C reads, `decoys.txt` or `alleles08.pairs.tsv`;
- `pairs.tsv.gz.target` no longer matches `09_hic_map/chains.fasta`;
- it is a symlink;
- compleasm is not configured.

## Not covered

- **`HiC.filtered.bam`** in builds made before stage 05 wrote pairs directly. It can be remade from `HiC.bam` with HapHiC's `filter_bam` ([LAYOUT](LAYOUT.md)); deleting it is the user's call.
- **Whole `archive/<time>/` folders.** They are the record of what a re-run replaced. The script removes only rebuildable files inside them; the folder itself goes only when the user agrees ([LAYOUT](LAYOUT.md), "Re-running").
- **Trial work directories and anything not made by the driver**, such as scratch FASTAs in `16_validate/compleasm/tmp/`. List them for the user and let them decide.

## Sizes seen

- **13Q010** (triploid, 2026-11): 32 files, 39.0 GB of 56 GB, including the first build's archive.
  - Hi-C BAM 20.1 GB
  - compleasm 9.6 GB
  - concatenations 5.0 GB
  - candidate FASTAs 1.9 GB
  - Hi-C target and index 2.8 GB
- **Busk** (diploid, collapsed): 25 files, 29.4 GB of 41 GB. `check` rebuilt every exact file identically.
