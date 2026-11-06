# Legacy: reference-guided pipeline

Not part of the workflow. This is the earlier pipeline that came before the graph-first process, kept for reference and for comparing with a baseline. It scaffolds with the reference as a guide, which the rules in `AGENTS.md` forbid for a final assembly.

These files are unchanged from their last use on 13Q010, except that `hic_phase.py` now needs `ANNOTARIUM_DIR` set, for `github.com/zkstewart/annotarium`. The shell scripts read a config, with paths resolved against it. Sample-specific configs were left out.

`run_pipeline.sh` is the entry point. Nothing here has been re-tested in this repository.
