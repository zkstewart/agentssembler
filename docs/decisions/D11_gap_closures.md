# D11. Which gaps need more than the automatic fill?

**Read:** the per-gap report (`14_readfill/*.report.tsv`) and `reports/12_fill.md`.

| Outcome | Meaning | Action |
|---|---|---|
| `reads` | Filled from reads | Accepted; audited in D18 if raw reads exist |
| `reads_overlap` | Flanks overlap; one trimmed | Accepted; audited in D18. **Can be wrong in tandem arrays** |
| `no_span` | No spanning read | Leave. On both genomes most gaps are this |
| `too_few_reads` | One read | Leave |
| `inconsistent` | Reads disagree | Leave |
| `local_realign_failed` | Reads leave a flank more than `--maxTrim` from its end | Leave. Raising the limit gained 1 gap and lost 2 on WatsonsChoice |
| `read_conflict` | Reads continue elsewhere | Report as a misjoin candidate. On 13Q010 chr3_h2 they fit repeat anchoring |
| `trim_too_long` | A trim would remove a whole component | Examine by hand |

**Action available:** with a shared flank on a collapsed sample, a read can get MAPQ 0 in spite of masking. If a gap reports many low-MAPQ spanning reads, list its scaffold in `edits/fill_retry.tsv` (`chr6_h2<TAB>1<TAB>reason`) and `reset 12_fill --keepDecision`. It closed one gap on WatsonsChoice with 28 of 28 reads.

**Do not loosen `--minIdentity` or `--maxTrim` globally.** Loosening gained 1 and lost 2.

**Direction check:** if a scaffold came out reversed against its homolog as a whole (seen in `15_validate`), list it in `edits/flip.txt`. It changes writing direction only.

**Look behind the gap, above all where a homolog is much longer.** `scripts/review/gap_tails.py` lists where reads go past each gap edge. `scripts/review/unplaced_components.py` lists connected unplaced tangles and the scaffolds they touch. On Busk chr7_h2, 14 of 14 q20 tails entered unplaced chain0099 at about 6.4 kb, along a graph path through dropped unitigs, in a 7.24 Mb tangle touching only chr7_h2.

**What to do with it:** report it. Nothing there closes a gap, since a tangle cannot be ordered, and placing the first chain trades risk for completeness: put it to the user at [D13](D13_user_decisions.md). Exclude the largest Hi-C contributor before reading a component's Hi-C: one chain with a repeat fraction of 0.90 carried 77% of the links and pointed at the wrong chromosome.
