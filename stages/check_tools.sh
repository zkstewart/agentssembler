#!/bin/bash
# Report which tools can be found. Missing required tools stop the run; missing optional ones disable the steps named.
fail=0
for v in AS_MINIMAP2 AS_SAMTOOLS AS_SEQKIT AS_BWA AS_SAMBLASTER AS_FILTER_BAM; do
  t="${!v:-}"; if command -v "$t" >/dev/null 2>&1; then echo "ok       $v=$t"; else echo "MISSING  $v=$t"; fail=1; fi
done
for v in AS_JUICER_TOOLS AS_COMPLEASM AS_PLOT_PYTHON; do
  t="${!v:-}"; if [ -n "$t" ] && { command -v "$t" >/dev/null 2>&1 || [ -e "$t" ]; }; then echo "ok       $v=$t"; else echo "optional $v not set: the steps using it are skipped"; fi
done
if [ -n "${AS_COMPLEASM:-}" ]; then
  d=$(dirname "$AS_COMPLEASM")
  for t in miniprot hmmsearch; do if PATH="$PATH:$d" command -v "$t" >/dev/null 2>&1; then echo "ok       $t"; else echo "MISSING  $t, needed by compleasm"; fail=1; fi; done
fi
python3 -c "import numpy" || fail=1
exit $fail
