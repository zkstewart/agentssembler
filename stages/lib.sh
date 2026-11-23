# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Sourced by every stage. The driver exports AS_* from the config and from recorded decisions.
set -euo pipefail
cd "$AS_WORK"
S="$AS_REPO/scripts"
export PATH="$(dirname "$AS_MINIMAP2"):$(dirname "$AS_SAMTOOLS"):$(dirname "$AS_SEQKIT"):$PATH"
# compleasm finds miniprot and hmmsearch on PATH, so its own environment has to be there too.
if [ -n "${AS_COMPLEASM:-}" ]; then export PATH="$PATH:$(dirname "$AS_COMPLEASM")"; fi
log() { echo "[$(date +%H:%M:%S)] $*"; }
need() { for f in "$@"; do [ -s "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }; done; }
# Stops the stage when free space where the work lives falls below the limit. A full disk can crash a machine.
guard() { local f; f=$(df -B1G "$AS_GUARD_PATH" | awk 'NR==2{print $4}'); [ "$f" -ge "$AS_MIN_FREE_GB" ] || { echo "ERROR: only $f GB free on $AS_GUARD_PATH (limit $AS_MIN_FREE_GB)" >&2; exit 3; }; }
# step OUTPUT command...: skip when OUTPUT exists. Long mappings write to .part and rename, so a partial file never counts.
# To use outputs mapped elsewhere (an HPC, say), put them in place and the step is skipped.
# A failed command's partial output is moved to OUTPUT.failed: left in place, the next run would skip the step and build on it.
step() {
  local o="$1"; shift
  if [ -s "$o" ]; then log "skip, exists: $o"; return 0; fi
  log "run: $o"; local rc=0; "$@" || rc=$?
  if [ "$rc" -ne 0 ]; then [ -e "$o" ] && mv "$o" "$o.failed"; echo "ERROR: step failed with exit $rc: $o" >&2; exit "$rc"; fi
  need "$o"
}
# map_paf OUT.paf.gz minimap2-args...: a guarded mapping that writes OUT only when complete
# A kept mapping (`keep` in stages.py) survives a reset without --remap. OUT.target holds the checksum of the sequence it was mapped to, so
# a kept mapping whose target has since changed stops the stage instead of being reused. Mappings made before this have no OUT.target.
target_sum() { md5sum < "$1" | cut -d' ' -f1; }
check_target() { # OUT TARGET
  [ -s "$1.target" ] && [ -f "$2" ] && [ "$(target_sum "$2")" != "$(cat "$1.target")" ] && {
    echo "ERROR: $1 was mapped to an earlier $2. Reset with --remap so it is made again." >&2; exit 4; }
  return 0
}
map_paf() {
  local o="$1"; shift; local t="${@: -2:1}"
  [ -s "$o" ] && { check_target "$o" "$t"; log "skip, exists: $o"; return 0; }
  guard; log "map: $o"
  "$AS_MINIMAP2" "$@" 2> "$o.log" | gzip -1 > "$o.part" &
  local p=$!
  while kill -0 $p 2>/dev/null; do
    local f; f=$(df -B1G "$AS_GUARD_PATH" | awk 'NR==2{print $4}')
    if [ "$f" -lt "$AS_MIN_FREE_GB" ]; then kill $p; pkill -P $$ minimap2 || true; echo "ERROR: disk guard stopped $o" >&2; exit 3; fi
    sleep 30
  done
  wait $p; gzip -t "$o.part"; mv "$o.part" "$o"
  if [ -f "$t" ]; then target_sum "$t" > "$o.target"; fi
}
agp2fa() { python3 "$S/agp_to_fasta.py" "$1" "$2" "$3.part"; mv "$3.part" "$3"; "$AS_SAMTOOLS" faidx "$3"; }
# has_rows FILE: true when an edit file holds more than comments
has_rows() { [ -s "$1" ] && grep -qv '^#' "$1"; }
collapsed() { [ "$AS_COLLAPSED" = "true" ]; }
has_ref() { [ -n "${AS_REFERENCE:-}" ]; }
hic_rescue() { [ "${AS_HIC_RESCUE:-false}" = "true" ]; }
# compleasm gene copies on the unitigs (stage 01), when compleasm is set: evidence for how many genes a decision gains or loses
GENES="01_genes/${AS_LINEAGE:-none}/full_table_busco_format.tsv"
has_genes() { [ -s "$AS_WORK/$GENES" ]; }
has_raw() { [ -n "${AS_RAW_READS:-}" ]; }
# raw_kind FILE: sam for BAM, CRAM or SAM (by extension); otherwise fastq or fasta, from the first character, gzipped or not
raw_kind() {
  case "$1" in
    *.bam|*.cram|*.sam) echo sam ;;
    *) local c; c=$(gzip -cdf -- "$1" 2>/dev/null | head -c1 || true); if [ "$c" = "@" ]; then echo fastq; else echo fasta; fi ;;
  esac
}
# raw_records FILE [IDS]: id, sequence and mean Q per read, tab-separated; only the reads named in IDS when given.
# Q is the qs:f tag of a BAM, the mean error probability of a FASTQ as a Phred score, and NA when unknown.
_raw_pick() { if [ -n "$2" ]; then "$AS_SEQKIT" grep -j 4 -f "$2" "$1"; else "$AS_SEQKIT" seq -j 4 "$1"; fi; }
raw_records() {
  local f="$1" ids="${2:-}"
  case "$(raw_kind "$f")" in
    sam) "$AS_SAMTOOLS" view -@ 2 ${ids:+-N "$ids"} "$f" | awk -F'\t' 'BEGIN{OFS="\t"} {q="NA"; for(i=12;i<=NF;i++) if($i~/^qs:f:/){q=substr($i,6); break}; print $1,$10,q}' ;;
    fastq) _raw_pick "$f" "$ids" | "$AS_SEQKIT" fx2tab -i -Q -q ;;
    fasta) _raw_pick "$f" "$ids" | "$AS_SEQKIT" fx2tab -i -Q | awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$2,"NA"}' ;;
  esac
}
# per-haplotype list, from the scaffold names in an AGP
haps_of() { awk -F'\t' '$1!~/^#/{n=split($1,a,"_"); if(a[n]~/^h[0-9]+$/)print a[n]}' "$1" | sort -uV; }
