# Sourced by every stage. The driver exports AS_* from the config and from recorded decisions.
set -euo pipefail
cd "$AS_WORK"
S="$AS_REPO/scripts"
export PATH="$(dirname "$AS_MINIMAP2"):$(dirname "$AS_SAMTOOLS"):$(dirname "$AS_SEQKIT"):$PATH"
log() { echo "[$(date +%H:%M:%S)] $*"; }
need() { for f in "$@"; do [ -s "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }; done; }
# Stops the stage when free space where the work lives falls below the limit. A full disk can crash a machine.
guard() { local f; f=$(df -B1G "$AS_GUARD_PATH" | awk 'NR==2{print $4}'); [ "$f" -ge "$AS_MIN_FREE_GB" ] || { echo "ERROR: only $f GB free on $AS_GUARD_PATH (limit $AS_MIN_FREE_GB)" >&2; exit 3; }; }
# step OUTPUT command...: skip when OUTPUT exists. Long mappings write to .part and rename, so a partial file never counts.
# To use outputs mapped elsewhere (an HPC, say), put them in place and the step is skipped.
step() { local o="$1"; shift; if [ -s "$o" ]; then log "skip, exists: $o"; else log "run: $o"; "$@"; need "$o"; fi; }
# map_paf OUT.paf.gz minimap2-args...: a guarded mapping that writes OUT only when complete
map_paf() {
  local o="$1"; shift
  [ -s "$o" ] && { log "skip, exists: $o"; return 0; }
  guard; log "map: $o"
  "$AS_MINIMAP2" "$@" 2> "$o.log" | gzip -1 > "$o.part" &
  local p=$!
  while kill -0 $p 2>/dev/null; do
    local f; f=$(df -B1G "$AS_GUARD_PATH" | awk 'NR==2{print $4}')
    if [ "$f" -lt "$AS_MIN_FREE_GB" ]; then kill $p; pkill -P $$ minimap2 || true; echo "ERROR: disk guard stopped $o" >&2; exit 3; fi
    sleep 30
  done
  wait $p; gzip -t "$o.part"; mv "$o.part" "$o"
}
agp2fa() { python3 "$S/agp_to_fasta.py" "$1" "$2" "$3"; "$AS_SAMTOOLS" faidx "$3"; }
# has_rows FILE: true when an edit file holds more than comments
has_rows() { [ -s "$1" ] && grep -qv '^#' "$1"; }
collapsed() { [ "$AS_COLLAPSED" = "true" ]; }
has_ref() { [ -n "${AS_REFERENCE:-}" ]; }
has_raw() { [ -n "${AS_RAW_BAMS:-}" ]; }
# per-haplotype list, from the scaffold names in an AGP
haps_of() { awk -F'\t' '$1!~/^#/{n=split($1,a,"_"); if(a[n]~/^h[0-9]+$/)print a[n]}' "$1" | sort -uV; }
