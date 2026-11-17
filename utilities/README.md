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
| D | duplicated | On a backbone, in a region with more backbone copies than the ploidy allows. Every copy of such a region counts, not only the extra one |
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
- **Switch and misjoin bins are flagged relative to each assembly's own typical bin.** A bin is flagged when two things are each more than 4× (`--fold`) the assembly's median bin: its contacts with one other sequence relative to its own long-range contacts, and those contacts as a share of all the bin's contacts. This adjusts for the Hi-C library, ploidy and scaffold length, which would otherwise shift the flags. The second condition was added after a spike-in test (`tests/crosstalk_spikein.py`, results in `docs/TESTING.md`) showed that the first alone also flags bins that have only lost their own long-range contacts.
- **Scores made before that change are kept.** They are in `xtalk_q<mapq>/`; scores under the current rule are written to `xtalk_q<mapq>.v2/`, and `compare` reads only those. Rerunning `score` reuses the mapping, so it takes minutes.
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

# `agp_compare.py`

This script compares two AGP files that scaffold the same set of input contigs or unitigs, such as the agentssembler AGP against the AGP of HapHiC or YaHS. It reports how far the two agree on which contigs belong together, their order and orientation, and where the gaps are. It needs nothing but the two AGPs: no reads, FASTA or reference genome. In its reference mode it instead compares one AGP with a reference genome, such as an earlier assembly of the same genome (see "Reference mode" below).

**It measures agreement, not correctness.** Neither AGP is treated as the truth. The first one (A) only supplies the scaffold names, so every result is a place where the two tools made different choices. To decide which tool is right at such a place, look at the Hi-C contact maps of `hic_crosstalk.py`, which does not score order and orientation within a chromosome itself.

```bash
python agp_compare.py agentssembler.agp haphic.agp --labels agentssembler haphic -o agp_compare -p Busk
```

Give the AGP with meaningful scaffold names (e.g. `chr1_h1`) first. The scaffolds of the second AGP (B, e.g. `group1`) are then each matched to the A scaffold they share the most sequence with. Only the Python standard library is required; `scipy` (for the one-to-one assignment) and `matplotlib` (for the dot plot) are used when installed.

## How the comparison is made

1. **Markers.** Every contig is cut at each boundary that either AGP uses. A piece used by both AGPs is a shared marker, weighted by its length. This makes the AGPs comparable when one tool uses whole contigs and the other trims or splits them, as agentssembler does when it joins unitigs by their overlaps.
2. **Scaffold matching.** A table of the bp shared between every A scaffold and every B scaffold gives each B scaffold its best A scaffold.
3. **Junctions.** Each AGP is reduced to the markers that are scaffolded in both. Two neighbouring markers form a junction, and it is conserved when the same two contig ends are joined in the other AGP.
4. **Blocks.** A run of conserved junctions is a block: a stretch with identical order and orientation in both AGPs. Each block is then classed by where it sits in the matched B scaffold.
5. **Gaps.** Each gap is looked up in the other AGP through the junction it lies in.

Sequence that only one tool scaffolded is counted in section 1 of the summary and then set aside, so a contig left out by one tool is never counted as an order or orientation difference.

## Output files

All files are written to `--outdir` with the `--prefix` in front.

| File | Content |
|---|---|
| `summary.txt` | The full report, also printed to the screen. Read this first. |
| `dotplot.png` | Every block drawn as a line, A scaffolds on the x axis and their matched B scaffolds on the y axis, coloured by block class. |
| `matches_B_to_A.tsv` | One row per B scaffold: its best A scaffold and the agreement values for that pair. |
| `matches_A_to_B.tsv` | One row per A scaffold: how much of it B scaffolded, and in which B scaffolds. |
| `contingency.tsv` | bp shared between every A scaffold and B scaffold, with unscaffolded sequence as a final row and column. |
| `blocks.tsv` | Every block with its coordinates in both AGPs and its class. Use this to find where a difference is. |
| `markers.tsv` | Every shared marker with its position and strand in both AGPs. |
| `junctions.tsv` | Every junction of A and of B, and what the other AGP does with the same two markers. `within_component` marks a junction that lies inside one contig, not at a join. |
| `gaps.tsv` | Every gap of A and of B, and its verdict in the other AGP. |

## Tips for interpretation

The summary has five sections.

**1. Component content** says how much of the sequence scaffolded by one tool is scaffolded, unscaffolded or absent in the other. "Part not used" means the contig is in the other AGP but not that part of it, such as a trimmed overlap. A large unscaffolded value shows that one tool placed sequence the other left as fragments.

**2. Scaffold matching** says whether the two tools formed the same groups.

| Value | Meaning |
|---|---|
| Reciprocal best pairs | B scaffolds whose best A scaffold also has them as its best B scaffold. When this equals the number of scaffolds, the two tools found the same chromosomes and haplotypes. |
| Purity of B scaffolds | The share of shared bp that lies in each B scaffold's best A scaffold. Below 100% means some B scaffolds combine sequence from more than one A scaffold. |
| Purity of A scaffolds | The same in the other direction. Below 100% means some A scaffolds are divided among several B scaffolds. |
| Chromosome-level purity | The same as B purity, but counting all haplotypes of a chromosome as one. It is reported when the A names have the form `<chromosome>_h<N>` (`--chrom_regex`). If it is higher than the B purity, the difference is sequence placed on another haplotype of the right chromosome. |

**3. Order and orientation** counts junctions and blocks. Each junction of one AGP gets one of these statuses in the other:

| Status | Meaning |
|---|---|
| `adjacent` | The same two contig ends are joined. The tools agree. |
| `adjacent_orientation_differs` | The two markers are neighbours in both, but one of them is flipped. |
| `same_scaffold_not_adjacent` | Both markers are in one scaffold of the other AGP, but not next to each other: an order difference. |
| `different_scaffolds` | The markers are in different scaffolds of the other AGP: a grouping difference. |

The breakpoint count is the number of junctions that are not `adjacent`. One flipped contig inside a scaffold breaks two junctions, and one contig moved elsewhere breaks three, so breakpoints overstate the number of events. The block classes are the better guide to how much sequence is affected:

| Class | Meaning |
|---|---|
| `collinear` | In the same place and orientation in both. |
| `inverted` | In the same place, but in the opposite orientation. |
| `relocated` | In the matched scaffold, but at a different position along it. |
| `other_scaffold` | In a B scaffold whose best match is a different A scaffold. The `detail` column of `blocks.tsv` says whether that is another haplotype of the same chromosome or another chromosome. |

**4. Gaps** gives each gap one of these verdicts in the other AGP:

| Verdict | Meaning |
|---|---|
| `gap_in_both` | The other tool makes the same join and also has a gap there. |
| `join_gapless_in_other` | The other tool makes the same join without a gap, either directly or through contigs that this tool did not scaffold. In `gaps.tsv`, `other_components_between` gives the number of such contigs, and `implied_gap_size` the amount of sequence the other tool has there beyond what this tool already has. A negative size means the two sides of the gap overlap in the other tool's layout. |
| `join_made_orientation_differs` | The same two markers are neighbours in the other AGP, but one is flipped. |
| `join_not_made_in_other` | The markers on either side of the gap are not neighbours in the other AGP. |
| `not_assessable` | The gap lies beyond the outermost shared marker of its scaffold. |

Agentssembler joins most unitigs through their graph overlaps, while Hi-C scaffolders put a gap at every join. A large `join_gapless_in_other` count for the other tool's gaps is therefore expected, and shows how many of its joins agentssembler confirms with sequence.

**5. Per-scaffold matches** lists one row per B scaffold:

| Column | Meaning |
|---|---|
| purity | The share of the B scaffold's shared bp that lies in its best A scaffold. |
| compl. | The share of the A scaffold's sequence that this B scaffold holds. It is lowered both by sequence B placed elsewhere and by sequence B left unscaffolded. |
| ori | Whether the B scaffold runs in the same direction as the A scaffold (`+`) or is reversed (`-`). This carries no meaning on its own, as the direction of a scaffold is arbitrary. |
| blocks | The number of blocks in the pair. One block means the two scaffolds are identical over their shared sequence. |
| collin. | The share of the pair's bp on its collinear backbone, i.e. in the same order in both. |
| orient. | The share of the pair's bp whose strand agrees with the direction of the pair. |
| inv, reloc | The numbers of inverted and relocated blocks. |
| other hits | Other A scaffolds with sequence in this B scaffold, in Mb. Any entry here marks a B scaffold that combines sequence of several A scaffolds. |

## Reference mode: comparing an AGP with a reference genome

With `--ref_paf`, one AGP is compared with a reference genome instead of a second AGP. This is meant for cases such as reassembling a genome that already has a chromosome-level assembly. The comparison stays at the level of the AGP's own contigs, so gapless joins are tested as well as gapped ones.

First align the contigs that the AGP is made of to the reference. Then give the alignments to the script along with the AGP:

```bash
minimap2 -x asm5 -t 16 reference.fa unitigs.fa > unitigs_vs_ref.paf
python agp_compare.py agentssembler.agp --ref_paf unitigs_vs_ref.paf --labels reference agentssembler -o agp_compare_ref -p Busk
```

The queries of the PAF must be the contigs named in the AGP (for agentssembler, the input unitigs), not the scaffolds. One PAF therefore serves every AGP built from the same contigs. `asm5` suits a reference of the same individual; use `asm10` or `asm20` for a more distant one.

The alignments say where the reference has each piece of each contig. That layout takes the role of A and the AGP is B, so every table and value described above keeps its meaning, with "A scaffold" read as "reference sequence". The alignments are used as follows:

| Rule | Reason |
|---|---|
| Secondary alignments are dropped | Only the best placement of each piece is used. |
| Alignments covering less than `--min_aln` bp of the contig are ignored (default 10,000) | Short alignments are mostly repeats. A contig shorter than this is never placed. |
| Where two alignments of one contig overlap, the better one keeps the shared part | Each piece of a contig is placed once. |
| Alignments with a MAPQ below `--min_mapq` are left out as ambiguous (default 10) | The piece fits two places almost equally well, as where two haplotypes are identical. Its placement would be a guess, so it counts neither for nor against the AGP. |

Points that differ from a comparison of two AGPs:

| Where | How to read it |
|---|---|
| Section 1, "placed" | How much of the AGP's scaffolded sequence has a unique place on the reference. The rest is split into ambiguous sequence and sequence with no alignment, and neither is compared. |
| Section 1, "Reference covered" | The share of the reference that the AGP's scaffolds account for. "Covered more than once" is reference sequence claimed by two scaffolded pieces, which points to duplicated sequence in the AGP. |
| Section 2 | Each AGP scaffold is matched to a reference sequence. With a phased reference whose names fit `--chrom_regex`, a chromosome-level purity above the B purity is sequence on another haplotype of the right chromosome. |
| Section 3, "inside one component" | Junctions that lie within a single contig, not at a join between contigs. A difference there was not made by scaffolding: the contig and the reference disagree. |
| Blocks, `other_scaffold` | `other_haplotype` in the `detail` column is a phase switch relative to the reference. `other_chromosome` is a join between chromosomes. |
| Section 4, `join_gapless_in_other` | The reference has the two sides of the gap next to each other, so it supports the join. `implied_gap_size` in `gaps.tsv` is the reference sequence between the two sides, less the AGP's own sequence there: an estimate of the true size of the gap. A negative size means the two sides overlap on the reference. |
| Section 4, `join_not_made_in_other` | The reference does not have the two sides of the gap next to each other. |
| Section 5, `compl.` | The shared bp as a share of the length of the reference sequence. |

Things to keep in mind:

- **The reference is not the truth unless you know it to be.** A difference means the AGP and the reference disagree. If the reference is another individual, real structural differences are reported in the same way as errors.
- **A phase switch can be the reference's.** Where the AGP follows one reference haplotype and then another, either assembly may have the switch. The Hi-C maps of `hic_crosstalk.py` can say which.
- **Ambiguous sequence is not tested.** In a genome with long stretches where the haplotypes are identical, much sequence may be left out. Section 1 says how much.
- **The reference is known only through the alignments.** Reference sequences with no alignment at all are not counted, and runs of N in the reference are not read.
- **Contigs missing from the FASTA you aligned have no placement.** For agentssembler these are the short `rf_*` fill pieces, which are not unitigs.
- **Unplaced reference sequences count as scaffolds.** Every reference sequence with a placed contig is treated as a scaffold. To restrict the comparison to the chromosomes, select them by name with `--scaffold_regex_a`.

## Some final notes

- **What counts as a scaffold.** An object with two or more components is a scaffold (`--min_components`), and anything smaller is unscaffolded. This suits tools like HapHiC that list every leftover contig as its own object. If a tool names all its objects alike, or a real chromosome is a single contig, select the scaffolds by name with `--scaffold_regex_a` or `--scaffold_regex_b`.
- **Only the direction of the whole B scaffold is adjusted, and it is chosen by order, not by strand.** A scaffold and its reverse complement are the same scaffold, so before a pair is measured the script decides whether to read the B scaffold forwards or backwards (the `ori` column). That is the only adjustment made. Nothing inside the scaffold is moved or flipped, so every difference between the two AGPs remains.
  - `collin.` and `orient.` measure two separate properties of each block. `collin.` asks whether the block is at the same *position* in the sequence of blocks. `orient.` asks whether the block itself is on the same *strand*. A block can be in the right position and on the wrong strand: this is an inversion in place.
  - Reading the B scaffold backwards reverses the order of all its blocks and flips the strand of all of them, together. The single choice of direction therefore cannot repair a block that is flipped while its neighbours are not.
  - The direction is the one that puts the most bp in the same order as A. `orient.` is then measured against that direction, whatever its value turns out to be. If the direction were chosen by strand instead, `orient.` could never fall below 0.5, but an inverted block would then make the order of the whole scaffold appear reversed.
  - A high `collin.` with a low `orient.` therefore means: the blocks follow one another in the same order in both AGPs, but much of the sequence lies in blocks that are flipped where they stand. A value below 0.5 means that more than half of the pair's bp is in such blocks, which usually comes down to one or two large inverted blocks. The `inverted` rows of `blocks.tsv` say which.
  - As an example, take a pair of five blocks. Read backwards, all but one small block are in the same order as A, but two blocks holding 59% of the bp are flipped in place. Read backwards, `collin.` is 0.99 and `orient.` is 0.41. Read forwards, `orient.` would be 0.59, but only 56% of the bp would still be in order. The script reports the first reading: the order agrees, and two blocks are inverted.
- **Two descriptions can fit the same difference.** A scaffold made of two large blocks in swapped order is equally well described as both blocks being inverted in a reversed scaffold. The script reports whichever keeps more bp in order, so check `blocks.tsv` or the dot plot before describing a large difference as an inversion.
- **Swapping the two AGPs** changes the names, and can change which blocks are called `other_scaffold`, since that depends on the best match of each B scaffold. The junction, breakpoint and block counts stay the same.
- **Small contigs dominate the counts and large ones the bp.** Most differences between Hi-C scaffolders are short contigs flipped or shifted by a few places, where Hi-C signal is weakest. Compare the block counts with their bp before concluding that two scaffoldings differ substantially.
