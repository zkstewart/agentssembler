# D11. Which gaps need more than the automatic fill?

**Read:** the per-gap report (`14_readfill/*.report.tsv`) and `reports/12_fill.md`.

| Outcome | Meaning | Action |
|---|---|---|
| `reads` | Filled from reads | Accepted; audited in D18 if raw reads exist |
| `reads_overlap` | Flanks overlap; one trimmed | Accepted; audited in D18. **Can be wrong in tandem arrays** |
| `no_span` | No spanning read | Leave. On all three genomes most gaps are this |
| `too_few_reads` | One read | Leave |
| `inconsistent` | Reads disagree | Leave |
| `local_realign_failed` | Reads leave a flank more than `--maxTrim` from its end | Leave |
| `read_conflict` | Reads continue elsewhere | Report as a misjoin candidate. Once, they fitted repeat anchoring |
| `trim_too_long` | A trim would remove a whole component | Examine by hand |

**Action available:** with a shared flank on a collapsed sample, a read can get MAPQ 0 in spite of masking. If a gap reports many low-MAPQ spanning reads, list its scaffold in `edits/fill_retry.tsv` (`<scaffold><TAB>1<TAB>reason`) and `reset 12_fill --keepDecision`. It once closed a gap with 28 of 28 reads.

**Do not loosen `--minIdentity` or `--maxTrim` globally.** Loosening gained 1 gap and lost 2.

**Direction check:** if a scaffold came out reversed against its homolog as a whole (seen in `15_validate`), list it in `edits/flip.txt`. It changes writing direction only.

**Look behind the gap, above all where a homolog is much longer.** `14_readfill/gap_tails.tsv` lists where reads go past each gap edge (`17_raw_reads/gap_tails.tsv` with raw reads, after D11). `scripts/review/unplaced_components.py` lists connected unplaced tangles and the scaffolds they touch; `15_validate` writes the same table. The reports raise an alarm for a homolog at least 1 Mb and 5% shorter, and for tails that mostly enter one unplaced sequence. Tails entering a scaffold are only counted: on one genome most led to homologs at divergent ends, or repeats. A short homolog with no gap can be real: check its telomeres. On one scaffold, 14 of 14 q20 tails entered the same unplaced chain a few kb past the edge, along a graph path through dropped unitigs, in a tangle of several megabases touching only that scaffold.

**What to do with it:** report it. Nothing there closes a gap, since a tangle cannot be ordered, and placing the first chain trades risk for completeness: put it to the user at [D13](D13_user_decisions.md). Exclude the largest Hi-C contributor before reading a component's Hi-C: one chain with a repeat fraction of 0.90 carried 77% of the links and pointed at the wrong chromosome.
