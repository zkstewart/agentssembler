# Traps

Each of these cost real time.

| # | Trap | What to do |
|---|---|---|
| T1 | `pkill -f pattern` from a command containing that pattern kills your own shell, and a `while pgrep -f pattern` loop in that command never ends, because it finds itself | Bracket a letter, as in `pgrep -f "[r]egress.py"`, then kill by PID |
| T2 | A failed command can leave a stale output that you read as new. This happened with crosstalk and with an SV table | Check timestamps, and never suppress errors on a metric step |
| T3 | `contact_map.py` groups file is one line per chromosome: `chr1<TAB>chr1_h1,chr1_h2`. Other layouts crash or overwrite | `15_validate` writes it for you. `--out` must exist |
| T4 | A chain in `-` orientation reverses coordinates: position x sits at scaffold `start + length - x` | Lift through the AGP; do not do it mentally |
| T5 | PAFs against `target.fasta` use pre-fill coordinates, which differ from final ones | Use the AGP the PAF was mapped to. Stage 15 remaps to final coordinates |
| T6 | A glob matching two PAFs passes an extra argument | Name files exactly |
| T7 | `minimap2 -c` costs about 50% more CPU. Depth and end links do not need it | Omit it unless you need identities |
| T8 | minimap2's PAF identity without `-c` is not an identity (0.958 for raw reads that are about 0.99 by CIGAR) | Do not report it as one |
| T9 | Haplotype numbers do not link across chromosomes | Never interpret `h1` across chromosomes |
| T10 | A process can look complete when its outputs come from an earlier order | `agentssembler reset` handles this; do not run a stage's scripts by hand on outdated inputs |
| T11 | A raw FASTA can fill the disk and crash the machine, more so when the work disk is a file on another drive | Stage 13 streams. Set `guard_path` |
| T12 | Corrected read IDs carry `:N` when a read was split | Strip it to match raw reads, and key spanning by ID and length so that two fragments do not count as one spanning read |
| T13 | Pooled depth ratios are dominated by repeat arrays | Report medians per site |
| T14 | A direction flag can leave a scaffold reversed against its homolog, which reads as a huge inversion | The validation report flags it. Use `edits/flip.txt` |
| T15 | `chain_breaks.tsv`'s `cross_ratio` is 0 at satellites | Require two signals (D04) |
| T16 | `build_chains.py` numbering changes with any change upstream | Chain names apply to one run only. Never copy an edit file between samples, or across a change upstream |
| T17 | A step that crashed could leave a partial output, and the next run skipped it as finished. On Busk a crash in `build_chains.py` left 21 of 266 chains, and join verification ran on them, wasting a 70-minute mapping | `step` now moves it to `.failed`. If a run was interrupted by something other than a failed command, such as a kill or a crash of the machine, check the newest outputs before resuming |
| T18 | compleasm needs `miniprot` and `hmmsearch` on PATH, and a pipeline that only prints a message on failure hides it | The driver adds compleasm's directory to PATH, `check_tools` tests both, and stage 15 now stops on failure |
| T19 | Group and chromosome labels are numbered by size. A change upstream renumbers them, and an edit written against an old label silently moves a chain to another chromosome | Write `with:<chain>` in `group_edits.tsv`. After any reset above stage 07, compare groups by chain with the archived run |
| T20 | Inputs can vanish between stages, when someone frees disk space. An empty stand-in for missing reads ran without error up to stage 11, and changed the scaffolds, because stage 08 chooses bridges with reads | `run` now stops if an input is missing and warns if one is empty. Tell the user which read sets must stay. Treat a surprise in the results as a question about the inputs first |
