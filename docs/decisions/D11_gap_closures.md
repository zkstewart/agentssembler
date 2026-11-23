# D11. Which gaps need more than the automatic fill?

**Read:** the per-gap report (`14_readfill/*.report.tsv`) and `reports/12_fill.md`.

| Outcome | Meaning | Action |
|---|---|---|
| `reads` | Filled from reads | Accepted; audited in D18 if raw reads exist |
| `reads_overlap` | Flanks overlap; one trimmed | Accepted; audited in D18. **Can be wrong in tandem arrays** |
| `no_span` | No spanning read | Leave. On 13Q010, WatsonsChoice and Busk most gaps were this |
| `too_few_reads` | One read | Leave |
| `inconsistent` | Reads disagree | Leave |
| `local_realign_failed` | Reads leave a flank more than `--maxTrim` from its end | Leave |
| `read_conflict` | Reads continue elsewhere | Report as a misjoin candidate. On 13Q010 chr3_h2 they fit repeat anchoring |
| `trim_too_long` | A trim would remove a whole component | Examine by hand |

**Action available:** with a shared flank on a collapsed sample, a read can get MAPQ 0 in spite of masking. If a gap reports many low-MAPQ spanning reads, list its scaffold in `edits/fill_retry.tsv` (`chr6_h2<TAB>1<TAB>reason`) and `reset 12_fill --keepDecision`. It closed one gap on WatsonsChoice with 28 of 28 reads.

**Do not loosen `--minIdentity` or `--maxTrim` globally.** Loosening gained 1 gap and lost 2 on WatsonsChoice.

**Direction check:** if a scaffold came out reversed against its homolog as a whole, list it in `edits/flip.txt`. The `11_build` report raises this from `13_scaffold/homolog_direction.tsv`; `15_validate` repeats the check on the final files. It changes writing direction only, and `12_fill` applies it: listed before stage 12, it costs nothing; after, it costs `reset 12_fill`.

**Look behind the gap, above all where a homolog is much longer.** `14_readfill/gap_tails.tsv` lists where reads go past each gap edge (`17_raw_reads/gap_tails.tsv` with raw reads, after D11). `scripts/review/unplaced_components.py` lists connected unplaced tangles and the scaffolds they touch; `15_validate` writes the same table. The reports raise an alarm for a homolog at least 1 Mb and 5% shorter, and for tails that mostly enter one unplaced sequence. Tails entering a scaffold are only counted: on 13Q010 most led to homologs at divergent ends, or repeats. So are tails into unplaced sequence that duplicates placed sequence: the `also_placed` column counts reads whose tail aligns nearly as well to a placed scaffold. On 292, chain0285 took 4 of 4 such tails at MAPQ 0, a copy of placed chain0124. A low MAPQ with `also_placed` 0 means something else: a near-identical unplaced allele pair, like chain0091 and chain0092. A short homolog with no gap can be real: check its telomeres. On Busk chr7_h2, 14 of 14 q20 tails entered unplaced chain0099 at about 6.4 kb, along a graph path through dropped unitigs, in a 7.24 Mb tangle touching only chr7_h2.

**What to do with it:** report it. Nothing there closes a gap, since a tangle cannot be ordered by this process (a small one with a boundary on each side belongs to D16's tangle check, [TANGLES](../TANGLES.md); if one is found only now, the fix is a `reset` back through D16), and placing the first chain trades risk for completeness: put it to the user at [D13](D13_user_decisions.md). Exclude the largest Hi-C contributor before reading a component's Hi-C: one chain with a repeat fraction of 0.90 carried 77% of the links and pointed at the wrong chromosome.
