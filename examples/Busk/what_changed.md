# Where Busk disagreed with the process, and what changed

Busk was the first genome run by following the guides and making each decision, not replaying earlier ones. Each row is a point where following the repo as it stood would have given a worse result, or where it failed.

| # | What happened | Result if followed as written | Change |
|---|---|---|---|
| 1 | `build_chains.py` crashed on an overlap longer than a unitig | Stopped at stage 04 | Overlap carries into the next unitig. WatsonsChoice had one in the safe direction by luck, and 4 more unwalked |
| 2 | That crash left a partial AGP (21 of 266 chains), and resume treated it as finished | A 70-minute mapping on wrong input, wasted | `step` moves failed output aside; `r_chains` compares chain counts between stages |
| 3 | **"Redundant" unitigs were the graph path between chain ends.** 40 of 60 gaps had one | 60 gaps. With them: 8 before filling | `restore_connectors.py`, on by default |
| 4 | `hops: 8` missed a 9-chain allele path once connectors were present: one false duplicate, two chains unplaced | Those errors | Default 12; 16 gave the same |
| 5 | Chains ran through collapsed forks, so 2x hubs went to one haplotype only | 13 of 14 h2-only genes in one chr2 region | `review/chain_forks.py`, `review/embedded_runs.py`, D04 fork section. h1 missing genes went from 19 to 7 |
| 6 | D14 accepted nominations with too little measured interior | Two wrong shared chains (under 5 kb interior), which I recorded first and reversed | D14: require interior; `r_group` lists nominations as candidates |
| 7 | D17 was asked before its input existed | An uninformed decision | Moved to the stage 10 gate; guide states the exact command |
| 8 | Group edits named labels, which are numbered by size and shift after a break | Old edits would have misplaced chain0082 and the chr2 hub. Caught by comparing with archived groups | `with:<chain>` edits |
| 9 | Subtelomeric satellites grouped by Hi-C with the wrong chromosome | Wrong telomere count per group | D06 telomere ledger; `r_group` alarm |
| 10 | D01 expected a visible 2x peak; Busk has none, and a 4-40x shelf instead | Misreading as a low-coverage class | D01, D05 wording; report gives histograms |
| 11 | compleasm's tools were not on PATH and the stage carried on | Validation finishing without gene content | Tools checked; failure stops the stage |
| 12 | The report printed `Specific-pair homolog fraction: 0.0000`. The parser skipped every row starting with `c`, so every `chr` row. **The same bug would have hit every genome** | A perfect-looking false number. Real values are 3.2-4.8% | Fixed |
| 13 | Merged joins were never read-tested at final coordinates | New joins taken on trust | `validate/join_reads.py`: 671 of 741 have q20 reads; connector joins 116 of 132, others 555 of 609 |
| 14 | Nothing asked where reads go past a gap edge. 14 of 14 q20 tails at chr7_h2 enter chain0099 | A 7.24 Mb tangle unnoticed behind a 4.85 Mb short haplotype | `review/gap_tails.py`, `review/unplaced_components.py`; D11 |
| 15 | A low-confidence chain, chain0042, stayed low-confidence by Hi-C alone | Unresolved 4.6 Mb | D06: walk the graph from the chain's gap end first. It settled the question |
| 16 | A pooled "raw 40, corrected 3 q20 spanners" looked like a strong correction effect. 37 of the 40 are at one array | A wrong claim about correction, which I made in passing | Per-site reading, already a rule. Followed late |

## Where the process held

- **The 5 kb overlap rule.** I doubted it for two genuine overlaps of 2.6 and 3.5 kb that were left as gaps. Raw reads closed one at exactly the overlap length, and the other has one read. But a 2.6 kb overlap also links chain0019 and chain0038, on homologs, with Hi-C 8x against. Keep the rule; curate with evidence.
- **Refusing array closures.** The chr6_h2 gap had 306 raw spanners, and length disagreed by 8-14 kb.
- **Evidence order.** Hi-C pointed the tangle at chr4 until one 95% repeat chain was excluded; the graph had pointed at chr7_h2 throughout.

## Decisions

| Decision | Outcome |
|---|---|
| D01 | `hap_depth` 62, `repeat_depth` 273, `collapsed` true |
| D03 | 2 of 523 joins broken by rule, 3 weak joins kept |
| D04 | 2 fork breaks, `chain_breaks.txt` |
| D06, D14 | `group_edits.tsv` |
| D09, D10 | `order_links.tsv`, `order_edits.tsv` |
| D16 | 4 candidates, none placed |
| D11, D18 | 0 fills from corrected reads; 1 raw overlap closure, passing audit |
| D12 | No inverted block over 400 kb; 2 translocation candidates and 1 dosage run not supported, `sv_calls.tsv` |
