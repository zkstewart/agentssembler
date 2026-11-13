# hic_crosstalk.py

Note to agents: this file is not needed for running the `agentssembler` pipeline, and can be ignored. It is used solely for downstream comparison of scaffolders' results.

This script runs a pipeline involving mapping of Hi-C data to one or more scaffolding outputs, to measure their objective quality from the perspective of the Hi-C mappings. The output is a score which allows intuitive comparison of different assemblies of the same organism and from the same underlying data. The main intent is to assess whether the Hi-C supports the scaffolding choices. Practically, this means looking for signs of haplotype switches/mosaicism and misassemblies. Contiguity, in terms of fragmentation, is something which this script intentionally seeks to *not* evaluate; this is the domain of auNG calculation rather than Hi-C mapping. All efforts are taken for this comparison to be unbiased, involving several curation steps including:
- Masking of regions which are expected to unfairly penalise Hi-C mapping data for its inability to resolve the assembly structure, such as repetitive, collapsed, or homozygous sequence shared by multiple haplotypes.
- Penalising misjoins but ignoring fragmentation.

## The short summary line

`hic_crosstalk.py` produces a `comparison.tsv` file which encapsulates the most important details from the analysis. Within that file is a row titled 'hic_breakdown' which contains a BUSCO-inspired short summary of the Hi-C quality results. Every bp of the assembly is placed into one of five categories as below:

| | Category | Meaning |
|---|---|---|
| C | consistent | On a backbone, scorable for phasing, and Hi-C raises no conflict. **`hic_score` = C%** |
| X | conflict | In a flagged switch or misjoin bin, shown as `[switch, misjoin]` |
| D | duplicated | More copies than the ploidy allows |
| U | untestable | Collapsed, hemizygous, or no partner backbone |
| O | off-backbone | Fragments and unplaced sequence |

## Tips for interpretation

With respect to interpretation, this is a useful crib sheet:

| | A higher value means… |
|---|---|
| **C** consistent | **Better:** More of the assembly is chromosome-scale, has all its haplotype copies scaffolded, and Hi-C agrees with the layout. |
| **X** conflict | **Worse:** More sequence where Hi-C contradicts the assembly: a haplotype switch or a join between chromosomes. Each copy involved is counted. |
| **D** duplicated | **Worse:** This may indicate unpurged duplicates or false haplotypes. |
| **U** untestable | **Worse†:** More of the scaffolded sequence lacks a scaffolded partner copy, so its phasing is unable to be reliably checked. †See below for further interpretation noes. |
| **O** off-backbone | **Worse:** More sequence left as fragments or unplaced contigs, or missing when you use `--genome_size`. |

## What U means in practice

A region can be classified as **U**(ntestable) for one of these reasons:

1. **The other haplotype's copy is missing from the assembly.** It was collapsed into one copy or purged. The assembly is short on haplotype sequence. This may indicate a flaw in the assembly.
2. **The other copy exists but isn't on a backbone.** It's a fragment or an unplaced contig, so it's unscaffolded. This may indicate a flaw in the assembly.
3. **The organism really has only one copy there, or the copies are too different to align.** Examples are sex chromosomes, large insertions or deletions present on one haplotype only, and highly divergent or repeat-rich regions such as centromeres. This does not necessarily indicate a flaw in the assembly.

For assemblies of the same organism, cause 3 is roughly the same in all of them, so **a difference in U mostly reflects causes 1 and 2.** In plain terms:

> A higher U means more of this assembly's chromosome-scale sequence lacks a scaffolded partner copy. Either the other haplotype was lost or collapsed, or it was left unscaffolded. That region isn't shown to be wrong; the assembly just doesn't deliver it as a complete, checkable phased set.

- **U and O are linked.** If haplotype 2's copy of a region is left as a fragment, it counts as O, and haplotype 1's copy becomes U. One scaffolding gap therefore costs bp in two categories. This should be a fair penalisation, since neither copy ends up as a verified phased pair, but it's why U and O often move in opposite directions between tools.
- **To tell cause 1 from cause 2:** if `all_bp` is well below the expected size of all haplotypes, sequence was probably lost, which is collapse. If O is high, the sequence is present but unscaffolded.

## Some final notes

- **Off-backbone sequence counts against the score.** If it didn't, a tool could score well by leaving difficult sequence unscaffolded.
- **Use `--genome_size`** in `score` or `compare` (the expected size of all haplotypes combined) to give every assembly the same denominator. Otherwise a tool that throws sequence away gets a smaller denominator and a higher score, and one that outputs lots of debris is penalised.
- **Many misjoin flags are probably biological.** On some data tested while developing this script, several of the strongest flagged bins sat at chromosome ends, which tend to contact each other inside the nucleus. They should appear in both assemblies and cancel out when you compare them. Thus, these values are solely for ***interpretation in comparison to another scaffolder***, rather than being a truly objective indication of quality that would stand up to comparison of different organisms.
