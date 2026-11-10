# D13. What belongs to the user?

An agent stops here and asks. Do not record D13 yourself until the user has answered, and quote their answer in the reason.

**Hand over:** `16_validate/juicebox/out_JBAT.hic` and `.assembly`, `16_validate/sv/sv_calls.tsv`, `16_validate/curation_targets.tsv` (you write this one, under the header stage 15 puts there), the work directory's `README.md` filled in from its template, and the reports.

**The user decides:**

- Which assembly to adopt, when a change trades contiguity against placed sequence or risk.
- How much risk to accept for completeness. Tails ([D16](D16_tails.md)) and chains behind a gap ([D11](D11_gap_closures.md)) are the clearest cases. On all three genomes so far the user chose to err toward minimising assembly risk, including leaving Busk's chain0099 unplaced.
- Whether medium-confidence SVs stay as assembled after inspection.
- Juicebox edits. **This route has not been exercised:** bring edits back as breaks, links or group edits with reasons, and re-run from the earliest affected stage. Do not paste a Juicebox export over the AGP.

**Declare this limit in any report:** the reference was used as an alarm during development. It shaped where people looked, so "reference-free" is true of the evidence behind each decision, not of the whole process.
