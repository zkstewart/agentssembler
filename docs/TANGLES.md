# Tangles: the D16 tangle check, and TTT

At [D16](decisions/D16_tails.md), before choosing any tails, check whether a tangle in the unitig graph is hiding a haplotype's route. There are two tools:

- **Lost-route forks** (`12_shared/chain_forks.tsv`, rows `LOST`): written by stage 10 for collapsed samples. Each row is a chain that runs straight through a fork. The other exit leads to an unplaced, about 1x chain. The fix is a chain break.
- **[TTT](https://github.com/marbl/TTT)**: for a small tangle that bridges and tails cannot cross. TTT solves the tangle as a whole, in two steps:
  1. It gives every node an integer multiplicity that fits coverage, using a MILP.
  2. It finds the traversal that best fits read-to-graph alignments.

  The pipeline's own bridges work one junction at a time, with reads that map linearly, and they ignore links under 5 kb of overlap that have fewer than 3 reads. So TTT can reach sequence they cannot.

Both results enter as `edits/` rows, each with its reason. Then `reset` to the earliest stage they affect, and run on. **This is part of the process, not a second build:** the run still ends in one assembly. Record what the check found and did in D16's reason, including "nothing to do".

The evidence comes from one genome, GammaCitrus (diploid, 9 chromosomes). Its numbers are in [RECORD](RECORD.md), under "TTT on graph tangles" and "Lost-route forks". On GammaCitrus, TTT ran after the pipeline had finished. The signals it acted on were already visible at D16, which is why the check now sits there.

## Lost-route forks

**Background.** Stage 04 drops a unitig whose sequence is already held by kept copies. Beside a collapsed hub, the dropped unitig may have been one haplotype's only route through it. The hub is then left with one exit on each side, so `build_chains.py` runs a chain straight through the fork. The other haplotype's arm is lost. It stays unplaced, or the arm on the chain gets used by both haplotypes. Stage 06 restores dropped connectors only between chain ends, so it cannot see a fork in the middle of a chain. The stage 07 scan (`10_groups/chain_forks.txt`) can't separate these forks from ordinary bubbles until placement is known, which is why stage 10 repeats it.

**A `LOST` row means all of the following:**

- the hub unitig is 2x;
- the chain is placed, or lies on a tail candidate;
- the other exit reaches a chain that is neither placed nor on a tail;
- that chain is about 1x and at least 50 kb long.

Rows are sorted by the length of that chain. `break_at` is the D04 break that separates the hub from the chain's arm on the fork side.

**Act on a row when the lost chain is worth having:**

- **The chain's arm is reused.** It is `over` in `12_shared/arm_usage.tsv`, and the lost chain is that arm's likely allele. GammaCitrus's chain0073 → chain0081 (779 kb) was this case, and so were chain0093, chain0138 and chain0185.
- **The lost chain leads on** to further unplaced 1x or telomeric chains. GammaCitrus's chain0107 → chain0173 was this case: it was h2's route to the chr5 telomere.
- **Genes:** compleasm on the unplaced sequence later shows genes there that are missing from both haplotypes. You can't see that yet at D16, so judge from the first two signals.

**Leave a row alone** when the lost chain is a bubble arm that rejoins the same chain within a few hundred kb. Breaking there only trades a contiguous join for a gap. On GammaCitrus, the chr2c rows chain0083 → chain0193 and chain0089 → chain0240 were this case. Grouped as TTT proposed, those chains were misplaced by Hi-C.

**A hub with forks on both sides** may need a break on each side, so read the `.` rows for the same hub too. chain0107's hub utg000467l had exits on both sides.

**Edits, with each one's reason:**

1. `edits/chain_breaks.txt`: the row's `chain` and `break_at` ([D04](decisions/D04_chimeric_chains.md), the fork case). The part to the right of the break becomes `<chain>b`. Several breaks in one chain are allowed (parts `<chain>b`, `<chain>c`, ...).
2. Then `agentssembler reset WORK 06_evidence --remap`. The chains changed, so stage 06's read mapping is redone. The Hi-C mapping is kept, because breaks are applied after it.
3. Stage 06 now finds the dropped route between the new chain ends and restores it as connectors.
4. Re-decide D05-D10 on the new reports:
   - make the hub shared ([D14](decisions/D14_shared_chains.md));
   - group the lost chain with its haplotype if Hi-C did not ([D06](decisions/D06_wrong_group.md));
   - add an order link ([D10](decisions/D10_curation.md)) where the restored link is below the bridge threshold. That is where TTT helps.

## When to use TTT

Use it on a tangle that is **small** (tens of unitigs) and **has a boundary on every side**. A boundary is a long, about 1x, heterozygous unitig, and together they must cut the tangle off from the rest of the graph. Use it when one of these holds:

| Signal at D16 | Where to read it | What TTT can add | GammaCitrus |
|---|---|---|---|
| A group end with no telomere whose tails stop inside a small tangle (a 2-2 tangle on a diploid). Beyond it, unplaced ~1x chains lead to telomeric chains | `12_shared/tails.tsv`, `12_shared/telomeres.tsv` (unplaced chains with a telomere), `scripts/ttt/tangle_explore.py` from the tail's last unitig | Both haplotypes carried through to the telomere | **chr5.** Both haplotypes got about 1.7 Mb and a telomere, and 46 genes per haplotype became complete, 40-44 of them on these arms |
| A `LOST` row you act on, where the restored route is below the bridge threshold, or runs through a tangle rather than a single path | `chain_forks.tsv`, then `12_shared/graph.bridges.tsv` after the reset | Which exit the reads favour, and the junctions to link | **chr2a.** TTT preferred chain0081, but its reads lay on a segment with an exact copy. A weak result: test it as below |

**Do not use it for:**

- **Large tangles, or tangles with dead ends.** At a dead end, every path needs connections TTT bans ("Z connections"), and it reports "No valid path found". GammaCitrus's main chr2 tangle had 254 unitigs, and its MILP grew by about 340 MB a minute before stalling.
- **Phase:** a D07 tip choice, or phase between bubbles. TTT scores the swap the same, and lists it as a "swappable interval".
- **Reuse the graph forces:** TTT can only confirm it (GammaCitrus chr4 chain0177).
- **Bubble alleles at a chain end that already has a direct graph join.** Adding one replaces a join with a gap.

**Tools.** TTT and GraphAligner are not part of the main environment. Put each in its own conda environment:

- `graphaligner`: GraphAligner 1.0.20, from bioconda.
- `ttt`: TTT's repository, PuLP 2.9 and GLPK. PuLP 4 removed an API that TTT uses.

Installing them is the user's call. Ask the first time a signal calls for TTT. If they decline, record that in D16 and go on with tails.

## TTT procedure

Work in `WORK/ttt/`, the one place in the work directory for this ([LAYOUT](LAYOUT.md)). `W` is the work directory, `S` the repository, and `HAP` is D01's haplotype depth.

1. **Find the tangle and its boundaries.** Start from the unitig at the problem: the last unitig of a stopped tail, or the hub of a `LOST` row.

   ```bash
   python3 $S/scripts/ttt/tangle_explore.py $W/00_inputs/unitigs.gfa $W/04_copynumber/depth.unitigs.tsv $W/12_broken/chains.agp \
     $W/12_shared/pre/SAMPLE.placement.tsv utg000509l,utg000124l --hap HAP --out $W/ttt/tangles/chr5
   ```

   The output lists the anchors it met, which are the boundary candidates, the side of each that faces the tangle, and where each sits in the pre-build. Write a boundary file with one `in<TAB>out` line per haplotype. TTT does not decide which entrance pairs with which exit, so for a 2-2 tangle write both pairings.

2. **Coverage.** Use corrected-read depth, not hifiasm's `rd:i`, which set real alleles to multiplicity 0 on GammaCitrus.

   ```bash
   awk -F'\t' '$1!~/^#/{print $1"\t"$3}' $W/04_copynumber/depth.unitigs.tsv > $W/ttt/coverage.tsv
   ```

3. **Reads.** Take only the reads that touch the tangle, and align them to the whole graph, so that look-alike copies elsewhere still compete for them. On GammaCitrus, 186,250 reads took 1 h 54 min on 16 threads and used 5.7 GB.

   ```bash
   python3 $S/scripts/ttt/tangle_reads.py $W/02_ont_map/ont_vs_unitigs.paf.gz $W/ttt/tangles/chr5 $W/ttt/ids.txt
   seqkit grep -f $W/ttt/ids.txt $W/00_inputs/ont_corrected_reads.* | seqkit fq2fa -o $W/ttt/reads.fasta
   GraphAligner -g $W/00_inputs/unitigs.gfa -f $W/ttt/reads.fasta -a $W/ttt/reads.gaf -x vg -t 16
   ```

   Delete `reads.fasta` afterwards: the GAF is what you keep.

4. **Run TTT,** once per pairing:

   ```bash
   python TTT.py --graph $W/00_inputs/unitigs.gfa --alignment $W/ttt/reads.gaf --coverage $W/ttt/coverage.tsv --median-unique HAP \
     --boundary-nodes $W/ttt/chr5_A.tsv --outdir $W/ttt/runs/chr5_A --milp-time-limit 1500
   ```

   Rerun with `--num-initial-paths 30`. Paths that change between runs are guesses; on GammaCitrus they were identical.

5. **Accept or reject each run:**

   ```bash
   python3 $S/scripts/ttt/check_traversal.py $W/00_inputs/unitigs.gfa $W/ttt/coverage.tsv $W/12_broken/chains.agp \
     $W/12_shared/pre/SAMPLE.placement.tsv $W/ttt/runs/chr5_A --hap HAP
   ```

   - **Every adjacency must be a GFA `L` line.** TTT reads hifiasm graphs heuristically.
   - **Multiplicity often disagrees with depth.** That alone does not reject a run: where the graph gives both haplotypes no other way through a 0.5-1x node, TTT gives it multiplicity 2. Reject the run if a disagreement is the change you would make.
   - **Tied pairings mean the pairing is unknown.** Break the tie with `scripts/ttt/hic_phase.py` and report it as weak. On GammaCitrus it was 41 pairs against 24.
   - **Check where the supporting reads lie.**
     - Look for exact copies among the exit unitigs in `03_unitig_align/unitigs_self.paf.gz`.
     - Count the reads in `02_ont_map/ont_vs_unitigs.paf.gz` that run from the hub past any shared segment into each exit.
     - Reads on a segment with an exact copy support both copies. On GammaCitrus the read counts disagreed with TTT at chr2a: 157 raw q20 reads ran past the duplicate into the arm, and none into chain0081.

6. **Translate to edits.** The script only proposes rows. Copy the ones you accept, each with its reason.

   ```bash
   python3 $S/scripts/ttt/ttt_to_edits.py $W/ttt/runs/chr5_B --work $W --gaf $W/ttt/reads.gaf --depth $W/ttt/coverage.tsv --hap HAP --sample SAMPLE --tag chr5
   ```

   It proposes groups only for chains that are ungrouped and unplaced in the pre-build. It proposes order links only for junctions with a graph link. It reports what chains cannot express: an "inside-chain switch".

   **TTT rows through a tangle replace the tails through it. Do not accept both.**

   | TTT result | Edit | Guide | Reset from |
   |---|---|---|---|
   | An unplaced chain on one haplotype's traversal | `group_edits.tsv`: `group with:<anchor>` | D06 | `08_bridges` |
   | Consecutive chains with a graph link | `order_links.tsv`, type `ttt` | D09 | `09_order` |
   | A link with under 5 kb overlap and fewer than 3 reads, which bridges ignore | `order_links.tsv`, type `ttt`, with the GraphAligner read count in the reason. Expect a gap there | D10 | `09_order` |
   | The path leaves a chain partway along (an inside-chain switch) | `chain_breaks.txt`, as a lost-route fork | D04 | `06_evidence --remap` |
   | The path runs through unitigs that stage 04 dropped, between ends that are not chain ends | Cannot be expressed. It stays a gap | | |

7. **Reset and run on.**
   - Reset from the earliest stage in the table. With chain breaks, that means `06_evidence --remap`, then run `ttt_to_edits.py` again on the new chains before writing group and order rows.
   - Later decisions are set aside. Re-decide each on its new report.
   - At D16 again, choose tails on the new pre-build. Record what TTT changed, and its weakest change, in D16.
   - List any weak result as a curation target at D13, as for any other judgement.

## GammaCitrus: what TTT gave

Of six tangles, one gain: the 17-unitig chr5 telomeric tangle, which took both haplotypes telomere to telomere. The larger tangles gave no path, or the pipeline's own answer. The D16 check would have listed both of GammaCitrus's lost routes before TTT ran. Each tangle's result and the net change are in [RECORD](RECORD.md), "TTT on graph tangles" and "Lost-route forks".
