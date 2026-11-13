Note to agents: this file is not needed for running the `agentssembler` pipeline, and can be ignored. It is used solely for downstream comparison of scaffolders' results.


# `hic_crosstalk.py`

This script runs a pipeline involving mapping of Hi-C data to one or more scaffolding outputs, to measure their objective quality from the perspective of the Hi-C mappings. The output is a score which allows intuitive comparison of different assemblies of the same organism and from the same underlying data. The main intent is to assess whether the Hi-C supports the scaffolding choices. Practically, this means looking for signs of haplotype switches/mosaicism and misassemblies. Detailed contiguity (N50, auNG and the like) is left to `scaffold_metrics.py`; here, sequence left unscaffolded is only counted as not delivered (category O below), so that a tool cannot score well by leaving difficult sequence out of its chromosomes. All efforts are taken for this comparison to be unbiased, involving several curation steps including:
- Labelling every assembly the same way. A reference is used only to give each sequence a common chromosome name and position in the maps; haplotype labels from any tool (including the reference) are never used.
- Finding homologous copies from each assembly's own sequences, irrespective of their chromosome names, so that real translocations in one haplotype are not penalised.
- Removing read pairs that cannot be placed uniquely (MAPQ < 20), such as those from repeats or from regions where the haplotypes are near-identical.
- Masking regions where the Hi-C data could not reveal a phasing error even if one existed, such as collapsed or hemizygous sequence, or sequence whose other copies are not on chromosome-scale scaffolds. Regions with more copies than the ploidy allows are also kept out of the phasing score, and reported separately.
- Never masking contacts between chromosomes, so that misjoins stay visible.

## The short summary line

The table written as output by `hic_crosstalk.py compare` has a row titled `hic_breakdown`, which contains a BUSCO-inspired short summary of the Hi-C quality results. Every bp of the assembly (or of `--genome_size`, if given) is placed into exactly one of five categories, which sum to 100%:

| | Category | Meaning |
|---|---|---|
| C | consistent | On a backbone (chromosome-scale scaffold), scorable for phasing, and Hi-C raises no conflict. **`hic_score` = C%** |
| X | conflict | In a 100 kb bin flagged as a haplotype switch or misjoin, shown as `[switch, misjoin]` |
| D | duplicated | On a backbone, in a region with more backbone copies than the ploidy allows |
| U | untestable | On a backbone, but no complete set of homologous copies on other backbones: collapsed, hemizygous, or the other copy is not scaffolded |
| O | off-backbone | Fragments and unplaced sequence (and missing sequence, when `--genome_size` is larger than the assembly) |

## Tips for interpretation

With respect to interpretation, this is a useful crib sheet:

| | A higher value means… |
|---|---|
| **C** consistent | **Better:** More of the assembly is chromosome-scale, has all its haplotype copies scaffolded, and Hi-C agrees with the layout. |
| **X** conflict | **Worse:** More sequence where Hi-C contradicts the assembly: a haplotype switch or a join between chromosomes. Each copy involved is counted. |
| **D** duplicated | **Worse:** This may indicate unpurged duplicates or false haplotypes that were built into the chromosomes. Duplicates left unplaced count as O instead. |
| **U** untestable | **Worse†:** More of the scaffolded sequence lacks a scaffolded partner copy, so its phasing is unable to be reliably checked. †See below for further interpretation notes. |
| **O** off-backbone | **Worse:** More sequence left as fragments or unplaced contigs, or missing when you use `--genome_size`. |

X is an approximate size, not an exact count of wrongly placed bp. A switch flags both haplotypes, and the flagged area extends up to ~2 Mb (`--far`) either side of the switch point. To count events, use the `switch_runs` and `misjoin_runs` rows, or the `switch_runs.tsv` and `misjoin_runs.tsv` files, which say where each one is.

Two contact rates are reported alongside the breakdown, for context:
- `homolog_fraction`: the share of within-chromosome contacts that land on another haplotype's copy, over the phase-testable regions. Lower is better.
- `interchrom_fraction`: the share of contacts between backbones that join different chromosomes. Most of this is background noise from the Hi-C library, so only differences between assemblies of the same data mean anything.

## What U means in practice

A region can be classified as **U**(ntestable) for one of these reasons:

1. **The other haplotype's copy is missing from the assembly.** It was collapsed into one copy or purged. The assembly is short on haplotype sequence. This may indicate a flaw in the assembly.
2. **The other copy exists but isn't on a backbone.** It's a fragment or an unplaced contig, so it's unscaffolded. This may indicate a flaw in the assembly.
3. **The organism really has only one copy there, or the copies are too different to align.** Examples are sex chromosomes, large insertions or deletions present on one haplotype only, and highly divergent or repeat-rich regions such as centromeres. This does not necessarily indicate a flaw in the assembly.

A translocation is not a cause of U. Homologous copies are found wherever they are in the assembly, whatever chromosome they are named after.

For assemblies of the same organism, cause 3 is roughly the same in all of them, so **a difference in U mostly reflects causes 1 and 2.** In plain terms:

> A higher U means more of this assembly's chromosome-scale sequence lacks a scaffolded partner copy. Either the other haplotype was lost or collapsed, or it was left unscaffolded. That region isn't shown to be wrong; the assembly just doesn't deliver it as a complete, checkable phased set.

- **U and O are linked.** If haplotype 2's copy of a region is left as a fragment, it counts as O, and haplotype 1's copy becomes U. One scaffolding gap therefore costs bp in two categories. This should be a fair penalisation, since neither copy ends up as a verified phased pair, but it's why U and O often move in opposite directions between tools.
- **To tell cause 1 from cause 2:** if `all_bp` is well below the expected size of all haplotypes, sequence was probably lost, which is collapse. If O is high, the sequence is present but unscaffolded.
- **U sequence is still checked for misjoins.** The misjoin scan covers every backbone bin, so U sequence is not placed on the wrong chromosome as far as Hi-C can tell. What cannot be checked is its phasing. Order and orientation within a chromosome are not scored by any category; view the contact maps for those.

## Some final notes

- **Off-backbone sequence counts against the score.** If it didn't, a tool could score well by leaving difficult sequence unscaffolded.
- **Use `--genome_size`** in `score` or `compare` (the expected size of all haplotypes combined, from an independent estimate such as k-mers) to give every assembly the same denominator. Otherwise a tool that throws sequence away gets a smaller denominator and a higher score, and one that outputs lots of debris is penalised.
- **Switch and misjoin bins are flagged relative to each assembly's own typical bin.** A bin is flagged when its contacts with one other sequence, relative to its own long-range contacts, are more than 4× (`--fold`) the assembly's median bin. This adjusts for the Hi-C library, ploidy and scaffold length, which would otherwise shift the flags.
- **Many misjoin flags are probably biological.** On some data tested while developing this script, several of the strongest flagged bins sat at chromosome ends, which tend to contact each other inside the nucleus. They should appear in both assemblies and cancel out when you compare them. Thus, these values are solely for ***interpretation in comparison to another scaffolder***, rather than being a truly objective indication of quality that would stand up to comparison of different organisms.

# `compleasm_copies.py`

This script counts BUSCO gene copies with the ploidy as the expected number, using the sequence labels from `hic_crosstalk.py prepare` so that every tool is judged the same way, without haplotype labels.

Run compleasm once per tool, on everything the tool output. The `hic_crosstalk.py` mapping target with its decoys removed is exactly that, and its sequence names match the labels:

```bash
awk '/^>/{keep = ($1 !~ /^>decoy_/)} keep' xt/<tool>/prepare/mapping_target.fa > <tool>.fa
compleasm run -a <tool>.fa -o compleasm/<tool> -l eudicotyledons_odb12 -t 12
python compleasm_copies.py --ploidy 2 --run <tool> compleasm/<tool>/eudicotyledons_odb12/full_table.tsv xt/<tool>/prepare/labels.v2.tsv [--run ...] --detail copies.tsv
```

compleasm's own Single/Duplicated split doesn't mean much for a pooled phased assembly, where every gene should be "duplicated". The script ignores that split and counts complete copies of each gene. The categories sum to 100%, and chromosomal + unscaffolded = C:

| | Meaning | Matching Hi-C category |
|---|---|---|
| **C** | Exactly ploidy copies, the expected state | |
| ↳ chromosomal | All copies are on backbones (chromosome-scale scaffolds), whatever chromosome they are named after | C |
| ↳ unscaffolded | One or more copies are on a fragment or unplaced sequence | O, with the backbone copies as U |
| **D-** | Fewer copies than the ploidy (a haplotype copy was collapsed or lost) | U |
| **D+** | More copies than the ploidy (false duplication or unpurged haplotigs) | D, or O if the extra copies are unplaced |
| **F / M** | Fragmented or missing | |

The output is a BUSCO-style line, e.g. `C:95.3%[chromosomal:95.1%,unscaffolded:0.1%],D-:2.0%,D+:2.5%,F:0.1%,M:0.2%,n:2805`. The `--detail` file lists every gene's copies, with the sequence, class and chromosome of each.

Pooling every sequence in one run means a weaker copy of a gene could be outcompeted by a stronger one and go uncounted. Every tool is treated the same way, so this does not favour any tool, and haplotype copies of a gene are usually nearly identical. If D- looks suspiciously high, check it against a per-haplotype run.

# `scaffold_metrics.py`

This script measures contiguity, which `hic_crosstalk.py` deliberately leaves alone. It has two subcommands, one for each question.

## `whole`: contiguity of everything a tool output

All of a tool's output files (every haplotype, plus unplaced sequence) are pooled, so tools that split their output differently are measured the same way. The input is a directory of annotarium `fasta stats` tables named `<genome>.[<part>.]<tool>.stats.tsv`. For example, `Busk.hap1.yahs`, `Busk.hap2.yahs` and `Busk.unplaced.yahs` are pooled as YaHS. The input contig sets (`p_utg`, `p_ctg`) are listed alongside as baselines.

| Value | Meaning |
|---|---|
| **auNG** | **The main measure.** The area under the NGx curve: `sum(L²) / G`, where G is the total length of the shared input unitigs (`<genome>.p_utg`) unless `--genome_size` is given. Unlike N50, it uses every sequence length rather than one threshold, so small changes do not make it jump. G is the same for every tool of a genome, so the auNG ranking does not depend on the choice of G. Higher is better. |
| auN | The same, divided by the tool's own total length. It rewards a tool that discards short sequence, so prefer auNG for ranking. |
| N50 / NG50, L50 / LG50, N90 / NG90 | The usual threshold measures, for familiarity. NG rankings can depend on G. |
| total_bp, pct_of_G | How much sequence was output, against G. |
| pct_of_input | Output length as a % of the contig set that tool was given (`--input_map`, default `p_utg`). Well below 100% means sequence was discarded; above 100% means sequence was duplicated or added (gap Ns also add a little). |
| n_gaps | Number of scaffolding gaps, assuming each is one run of exactly `--gap_size` Ns (default 100). The script warns when that assumption cannot hold. |

## `placed`: contiguity of the chromosome-scale sequence

This uses the `labels.v2.tsv` file of `hic_crosstalk.py prepare`, so every sequence has the same chromosome naming and class (backbone, fragment or unplaced) for every tool. For each chromosome, the ploidy longest backbones are taken as its copies ("anchored"). Any further backbones and fragments are pieces that were not joined into those copies.

| Value | Meaning |
|---|---|
| **anchored_pct_output / _input / _G** | **The main measure.** How much sequence is in the chromosome copies, as a % of the tool's output, of its input contigs, or of G. Higher is better. |
| copies vs expected_copies | Chromosome copies found vs chromosomes × ploidy. Fewer means a haplotype copy of some chromosome is missing or too fragmented to reach backbone size. |
| pieces_per_copy | All pieces (backbones and fragments) per expected copy. 1.0 means every chromosome copy is in one piece. |
| extra_backbones, fragments | Chromosome-assigned pieces that were not joined into the copies. |
| anchored_gaps, contigs | With `--stats`: gaps within the copies, and contigs per chromosome. More gaps for the same length means more joins were made, each of which could be wrong. |
| busco_* | With `--compleasm`: gene copies counted within the anchored copies only. `busco_all_copies` is genes with at least ploidy copies, `busco_excess` more than ploidy, and `busco_copy_ratio` the mean of min(copies, ploidy) / ploidy. `compleasm_copies.py` gives the simpler, whole-assembly summary of these. |
| possible_chimeras, anchored_dom_bp | `possible_chimeras` counts copies with less than 90% of their aligned sequence on one reference chromosome. `anchored_dom_bp` is the anchored length with each copy's share aligning to other chromosomes subtracted. **Caution:** a real translocation produces the same signal as a misjoin, so these values are not a measure of error. Check the Hi-C maps and the `hic_crosstalk.py` misjoin runs before treating them as errors. |

Per-chromosome values are in `placed.per_chromosome.tsv` and per-assembly values in `placed.summary.tsv`.
