# Decision guides

The driver stops at each decision below and points to its guide. Read the stage report and the guide, look at the evidence, make any edits, and record the decision with a reason:

```bash
agentssembler decide WORK D04 --by NAME --reason "No chain has two signals. The one flagged chain has cross_ratio alone."
```

Decisions are recorded in `WORK/decisions/`, with a checksum of every edit file at that moment. An edit always needs its own reason as well.

| Decision | Question | Stops after | Judgement needed |
|---|---|---|---|
| [D01](D01_haplotype_depth.md) | What is the haplotype depth? Is there collapsed sequence? | `03_depth` | Reading a histogram |
| [D03](D03_weak_joins.md) | Which weak joins are real? | `04_chains` | Usually none beyond the rule |
| [D04](D04_chimeric_chains.md) | Is a chain chimeric? | `05_hic_chains` | High |
| [D05](D05_group_count.md) | Why not the expected number of groups? | `07_group` | Medium |
| [D06](D06_wrong_group.md) | Is a chain in the wrong group? | `07_group` | Medium |
| [D14](D14_shared_chains.md) | Is a chain shared between haplotypes? | `07_group` | High, collapsed samples only |
| [D07](D07_telomere_tips.md) | Should a telomeric tip be attached? | `08_bridges` | Medium |
| [D09](D09_order.md) | Is an odd order real or a search artefact? | `09_order` | High |
| [D10](D10_curation.md) | Should a junction be set by hand? | `09_order` | High, and rarely yes |
| [D17](D17_phase_switches.md) | Is a gap end a phase switch? | `10_arms_tails` | Medium, collapsed samples only |
| [D15](D15_arm_phase.md) | Which arm goes where? | `10_arms_tails` | Reporting, mostly |
| [D16](D16_tails.md) | Should a tail be placed? | `10_arms_tails` | High |
| [D11](D11_gap_closures.md) | Which gaps need more than the automatic fill? | `12_fill` | Low |
| [D18](D18_raw_read_joins.md) | Which read-derived joins survive raw reads? | `13_rawreads` | Medium, needs raw BAMs |
| [D12](D12_sv_calls.md) | Is an SV candidate real? | `15_validate` | High |
| [D13](D13_user_decisions.md) | What belongs to the user? | `15_validate` | The user's |

The numbers are gaps in the D-series left by merging: D02 (a failed self-alignment) is now a stage check, and D08 (bonus strength) is a setting explained in [SETTINGS](../SETTINGS.md).

To change an earlier answer: edit the file, then `agentssembler reset WORK STAGE` from the earliest stage it affects. Later stages and the decisions made from them are set aside, because they rest on the old result.
