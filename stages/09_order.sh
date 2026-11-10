#!/bin/bash
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
# Order each chromosome's groups in parallel by Hi-C with graph links and telomeres as near-hard constraints. Reference check is an alarm. Writes into 12_shared. Stops for D09, D10.
source "$AS_REPO/stages/lib.sh"
D=12_shared
step $D/order_links.base.tsv python3 "$S/make_order_links.py" $D/groups.tsv $D/end_links.tsv $D/graph.bridges.tsv $D/order_links.base.tsv
{ cat $D/order_links.base.tsv; if has_rows edits/order_links.tsv; then grep -v '^#' edits/order_links.tsv | cut -f1-3; fi; } > $D/order_links.tsv
mkdir -p $D/by_chrom; CH=10_groups/final.chroms.tsv
if [ ! -s $D/order.raw.tsv ]; then
  pids=()
  for c in $(cut -f2 $CH | sort -u); do
    if [ -s $D/by_chrom/$c.order.tsv ]; then continue; fi
    awk -v c="$c" 'NR==FNR{if($2==c)k[$1]=1;next} /^#/{print;next} ($3 in k)' $CH $D/groups.tsv > $D/by_chrom/$c.groups.tsv
    python3 "$S/hic_order.py" $D/by_chrom/$c.groups.tsv $D/pairs.tsv.gz $D/chains.fasta.fai $D/by_chrom/$c.order.tsv --endLinks $D/order_links.tsv --linkBonus "$AS_LINK_BONUS" --telomeres $D/telomeres.tsv --stats $D/chain_stats.tsv --minWithin "$AS_MIN_WITHIN" > $D/by_chrom/$c.log 2>&1 &
    pids+=($!)
  done
  for p in "${pids[@]:-}"; do [ -n "$p" ] && wait "$p"; done
  { head -n 1 $D/by_chrom/$(cut -f2 $CH | sort -u | head -1).order.tsv; for c in $(cut -f2 $CH | sort -u); do grep -v '^#' $D/by_chrom/$c.order.tsv; done; } | sort -s -t$'\t' -k1,1 > $D/order.raw.tsv
fi
if has_rows edits/order_edits.tsv; then python3 "$S/apply_order_edits.py" edits/order_edits.tsv $D/order.raw.tsv $D/groups.tsv $D/order.tsv $D/groups.edited.tsv; else cp $D/order.raw.tsv $D/order.tsv; fi
has_ref && python3 "$S/eval_order.py" $D/order.tsv $D/chains.agp 03_unitig_align/unitigs_vs_ref.paf.gz > $D/eval_order.txt || true
python3 "$S/review/unhonoured.py" $D/order.tsv $D/end_links.tsv $D/graph.bridges.tsv $D/telomeres.tsv 10_groups/final.chroms.tsv > $D/unhonoured.txt 2>&1 || echo "unhonoured.py failed: see $D/unhonoured.txt"
