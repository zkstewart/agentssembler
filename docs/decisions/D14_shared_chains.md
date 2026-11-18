# D14. Is a chain shared between haplotypes?

**Only for `collapsed: true`.** A shared chain is collapsed sequence that belongs in more than one haplotype of its chromosome: in both at ploidy 2, and in 2 to `ploidy` of them above that ("Dosage" below). `make_shared_copies.py` writes a copy `<chain>.<group>` for each haplotype that carries it.

**Read:** `11_ont_chains/chain_class.tsv`, `10_groups/shared.candidates.txt`, `10_groups/groups.shared.tsv`, and `scripts/review/hub_depth.py`, `scripts/review/reused_depth.py`.

**Yes when** interior depth (outside graph overlaps) is near 2x and Hi-C is tied between homologs.

**Not shared when** the chain is high copy: over about `ploidy` + 1 times haplotype depth (3x on a diploid), with Hi-C across all groups. Unassign it (WatsonsChoice chain0184, 262x, and chain0267, 100x).

**The score ratio alone is too strict.** Hub chains at 2x with a ratio of 1.5 were left unassigned or forced into one haplotype on WatsonsChoice. Read the call with depth and graph position. 11 chains were promoted by hand, for example:

`chain0175  shared  c8  Interior depth 64.6x; Hi-C to chr3 5.23x next chromosome (32,356 pairs), no separation between g17 and g18; 13 complete genes; graph-linked hub of the chr3 gap`

On Busk five short hubs (190-340 kb) had medians of 71-84x, depth class `specific`, with 50-130 kb overlaps at each end and 125-145x in between. The auto-candidates missed them, because class comes from the median. The report lists such chains from the middle third of their windows, **which only nominates them**.

**Then run `scripts/review/hub_depth.py`, and require real interior.** I first promoted seven on window depth, then reversed two: chain0171 and chain0206 had 4.7 and 4.4 kb outside graph overlaps, and about 360 Hi-C pairs. A depth over that little sequence says nothing. The five kept had 27-119 kb. Unassigned chains can still serve as bridge arms, so leaving one out loses little. Read degree too: chain0108 is tied in Hi-C but has degree 1,1 and depth about 60-88x, and was left alone.

**Why this matters:** promotion and 8-hop ordering, with tails, took genes complete only in unplaced sequence from 52 to 0. They were tested together, so their separate contributions are not known.

**Candidate list:** by default it comes from depth class. To replace it, add rows to `edits/shared.txt` (one chain per line, under its header) before `07_group` runs, or reset that stage. With no rows, the automatic list is used.

**Edit:** `edits/group_edits.tsv`, action `shared`, value `with:<a chain on that chromosome>`.

## Dosage: which haplotypes, above ploidy 2

At ploidy 2 a shared chain is in both haplotypes and there is nothing to choose. Above it, a chain at 2x of 4 is in two haplotypes, and copying it into all four would put sequence where it does not exist. **Untested on real data.**

`hic_group.py` sets each shared chain's dosage to its median depth over haplotype depth, from 2 to the number of haplotypes of its chromosome. A chain below that number is listed in `10_groups/groups.dosage.tsv` with the groups chosen to carry it. The `call` column says what separated them from the groups left out, in the evidence order:

| Call | Meaning | Trust |
|---|---|---|
| `graph` | Every chosen group has a member with a graph link to the chain, and no group left out has one | Highest |
| `alleles` | Every group left out holds an allele of the chain, and no chosen group does: a haplotype cannot carry both | High, if the allele is not a repeat hit |
| `hic` | Links per group depth to the weakest chosen group are at least twice those to the best group left out (`hic_ratio`) | Read the ratio |
| `arbitrary` | None of those. The groups are the best guess and no more | None: the report raises an alarm |

**Check before accepting:**

- **The dosage itself.** It comes from the median, which long overlaps pull down (the Busk hubs above). Run `scripts/review/hub_depth.py` on chains whose `median_depth` sits between two whole multiples.
- **Every `arbitrary` row.** Either find evidence and name the haplotypes, or unassign the chain. An unassigned chain can still serve as a bridge arm, at up to its depth's worth of haplotypes.
- **Warnings in `12_shared/copies.log`** after stage 08: a graph link from a group's member into a shared chain that has no copy in that group contradicts the call, and the graph outranks Hi-C.

**Edit:** `edits/group_edits.tsv`, action `shared_in`, value a comma list of the groups that carry the chain, each written `with:<a chain in that group>`:

`chain0175  shared_in  with:chain0012,with:chain0040  3.1x outside overlaps is wrong: 2.1x by hub_depth; graph links to both`

`shared` still means every haplotype, and replaces a dosage call. `10_groups/final.dosage.tsv` holds the calls after edits, and is what stage 08 reads.
