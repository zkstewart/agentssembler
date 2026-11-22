# D13. What belongs to the user?

An agent stops here and asks. Do not record D13 yourself until the user has answered, and quote their answer in the reason.

**Hand over:** `16_validate/juicebox/out_JBAT.hic` and `.assembly`, `16_validate/sv/sv_calls.tsv`, `16_validate/curation_targets.tsv` (you write this one, under the header stage 15 puts there), the work directory's `README.md` filled in from its template, and the reports.

**The user decides:**

- Which assembly to adopt, when a change trades contiguity against placed sequence or risk.
- **Whether to trial anything outside the process.** That means a different method or setting, a tool the guides do not name, or a fix specific to this assembly. Propose it here, after the run is complete, with the evidence for it. Do not try it first. If the user agrees, build it in a second work directory and give `agentssembler compare`'s table, where each gain lies, and the weakest change with how to revert it alone. The tangle check is not such a trial: it is part of D16.
- How much risk to accept for completeness. Tails ([D16](D16_tails.md)) and chains behind a gap ([D11](D11_gap_closures.md)) are the clearest cases. On all three genomes so far the user chose to err toward minimising assembly risk, including leaving Busk's chain0099 unplaced.
- Whether medium-confidence SVs stay as assembled after inspection.
- **Unlocalized scaffolds are part of the assembly** (`16_validate/unloc/`; `unloc_purge` 0.2): hand them over with the scaffolds, and `unplaced_rest.fasta` in place of the unplaced file. Give the user compleasm with and without them (`compleasm/hK_unloc/`), the Mb by phase evidence (most are `arbitrary_*`: phase unknown), and the duplication cost; leaving them out, or another `unloc_purge`, is the user's call. They add genes without touching the scaffolds; they do not replace placing a chain when the evidence allows it.
- With `hic_allele_rescue` on: the rescued regions for Juicebox review (orientation, and sequence on both haplotypes), from `10_groups/groups.rescued.tsv` and `sv/inv/merged.tsv`.
- Juicebox edits. **This route has not been exercised:** bring edits back as breaks, links or group edits with reasons, and re-run from the earliest affected stage. Do not paste a Juicebox export over the AGP.

**Declare this limit in any report:** the reference was used as an alarm during development. It shaped where people looked, so "reference-free" is true of the evidence behind each decision, not of the whole process.
