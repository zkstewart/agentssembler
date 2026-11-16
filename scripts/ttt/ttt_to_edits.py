#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Translate accepted TTT traversals into agentssembler chain-level edits (docs/TANGLES.md). Proposes rows only; a person
(or the agent at the decision point) copies them into edits/ with the reasons.

Per traversal (one haplotype path), unitig runs are collapsed into chain steps (12_broken/chains.agp) with orientation,
entry end and exit end (chainN+ = 3' end, chainN- = 5' end). The traversal's haplotype group is the group of its first
chain (the entrance anchor).

Proposed group_edits: a chain >= --minLen used by exactly one traversal of the tangle and unassigned in the
grouping (10_groups/final.tsv) and not placed in the stage-10 pre-build (12_shared/pre placement; bridge chains are placed
without a group and stay so) joins that traversal's group (`group with:<anchor>`). Chains already grouped are left as
they are: TTT marks allele swaps between bubbles as equally good ("swappable intervals"), so moving them is not a gain.
Chains used by both traversals are reported, with depth, and not made shared here (acceptance: multiplicity must agree
with depth; that is a separate call).

Proposed order_links: consecutive chain steps whose exit and entry unitigs are the chains' terminal unitigs at those ends
and whose ends have a graph link in 10_groups/end_links.tsv, when both chains end up in the traversal's group. Each with
the number of GraphAligner reads (mapq >= --minq) that contain the junction's two unitigs adjacently.

Run at D16, so tail chains are not yet placed: TTT rows through a tangle replace the tails through it; do not accept both.
Usage: ttt_to_edits.py run_dir[,run_dir...] --work agentssembler_workdir --gaf aln.gaf --depth coverage.tsv --hap HAP_DEPTH --sample SAMPLE [--tag chr5]'''
import argparse, collections, re, os
p = argparse.ArgumentParser()
p.add_argument("runs"); p.add_argument("--work", required=True); p.add_argument("--gaf", required=True)
p.add_argument("--depth", required=True); p.add_argument("--tag", default="ttt"); p.add_argument("--minLen", type=int, default=20000)
p.add_argument("--minq", type=int, default=20); p.add_argument("--hap", type=float, required=True); p.add_argument("--sample", required=True)
a = p.parse_args(); W = a.work
AGP = collections.defaultdict(list)
for l in open(f"{W}/12_broken/chains.agp"):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if s[4] == "W": AGP[s[0]].append((s[5], s[8]))
U2C = {u: (c, i, o) for c, us in AGP.items() for i, (u, o) in enumerate(us)}
CL = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(f"{W}/12_broken/chains.fasta.fai")}
GR = {}
for l in open(f"{W}/10_groups/final.tsv"):
    if l[0] != "#": s = l.rstrip("\n").split("\t"); GR[s[0]] = s[2]
SH = {l.split("\t")[0] for l in open(f"{W}/10_groups/final.shared.tsv") if l[0] != "#"} if os.path.exists(f"{W}/10_groups/final.shared.tsv") else set()
EL = {}
for l in open(f"{W}/10_groups/end_links.tsv"):
    if l[0] == "#": continue
    x, y, r, q, ov = l.rstrip("\n").split("\t"); EL[frozenset((x, y))] = (int(r), int(ov))
PLC = {}
for l in open(f"{W}/12_shared/pre/{a.sample}.placement.tsv"):
    if l[0] != "#": x = l.rstrip("\n").split("\t"); PLC[x[0]] = x[4]
D = {l.split("\t")[0]: float(l.split("\t")[1]) for l in open(a.depth)}

def steps_of(run):
    for l in open(os.path.join(run, "traversal.gaf")):
        s = l.rstrip("\n").split("\t"); path = next(x for x in s if x[:1] in "<>")
        yield s[0], re.findall(r"([<>])([^<>]+)", path)

def chain_steps(st):
    '''[(chain, orient, first_unitig_step, last_unitig_step, bp)] for runs of one chain; '-' for no chain.'''
    out = []
    for o, u in st:
        c, i, ao = U2C.get(u, ("-", -1, "+"))
        orient = "+" if (o == ">") == (ao == "+") else "-"
        if out and out[-1][0] == c and c != "-" and out[-1][1] == orient: out[-1][3] = (o, u); out[-1][4] += 1
        else: out.append([c, orient, (o, u), (o, u), 1])
    return out

def terminal(c, u, end):
    us = AGP[c]; return (us[-1][0] if end == "+" else us[0][0]) == u

# reads: adjacent unitig pairs per read path, mapq filtered
need = set(); runs = a.runs.split(","); trav = []
for run in runs:
    for name, st in steps_of(run):
        trav.append((run, name, st))
        for (o1, u1), (o2, u2) in zip(st, st[1:]): need.add(frozenset((u1, u2)))
RS = collections.Counter()
for l in open(a.gaf):
    s = l.split("\t")
    if int(s[11]) < a.minq: continue
    us = re.findall(r"[<>]([^<>]+)", s[5])
    for x in {frozenset((u1, u2)) for u1, u2 in zip(us, us[1:])} & need: RS[x] += 1

use = collections.defaultdict(set)
for run, name, st in trav:
    for c, *_ in chain_steps(st):
        if c != "-": use[c].add((run, name))
G, O, R = [], [], []
for run, name, st in trav:
    cs = chain_steps(st); anchor = cs[0][0]; grp = GR.get(anchor, "?")
    R.append(f"== {os.path.basename(run)} {name}: anchor {anchor} ({grp}); " + " ".join(f"{c}{o}" for c, o, *_ in cs if c != "-" and CL.get(c, 0) >= a.minLen))
    newgrp = {}
    for c, o, f, l, n in cs:
        if c == "-" or CL.get(c, 0) < a.minLen: continue
        if len(use[c]) > 1:
            R.append(f"   both traversals: {c} {CL[c]/1e3:.0f} kb, grouped as {GR.get(c, 'shared' if c in SH else 'none')}")
            continue
        if GR.get(c, "unassigned") == "unassigned" and c not in SH and c != anchor and not PLC.get(c, "dropped").startswith("placed"):
            newgrp[c] = grp
            dp = sum(D.get(u, 0) for u, _ in AGP[c]) / max(1, len(AGP[c])) / a.hap
            G.append(f"{c}\tgroup\twith:{anchor}\tD06 (TTT {a.tag}): on TTT's traversal from {anchor} (run {os.path.basename(run)}), used once; {CL[c]/1e3:.0f} kb, mean unitig depth {dp:.2f}x; unassigned and {PLC.get(c, 'dropped')} in the pre-build")
    for (c1, o1, f1, l1, n1), (c2, o2, f2, l2, n2) in zip(cs, cs[1:]):
        if "-" in (c1, c2): continue
        e1 = "+" if o1 == "+" else "-"; e2 = "-" if o2 == "+" else "+"
        if not (terminal(c1, l1[1], e1) and terminal(c2, f2[1], e2)):
            R.append(f"   inside-chain switch {c1}{o1}->{c2}{o2} at {l1[1]}->{f2[1]} (allele choice differs from the chain's)"); continue
        g1 = newgrp.get(c1, GR.get(c1)); g2 = newgrp.get(c2, GR.get(c2))
        ev = EL.get(frozenset((c1 + e1, c2 + e2)))
        rs = RS[frozenset((l1[1], f2[1]))]
        tag = f"{c1}{e1} {c2}{e2} reads(q>={a.minq}) {rs} graph {ev}"
        if g1 == g2 == grp and ev: O.append(f"{c1}{e1}\t{c2}{e2}\tttt\tD09 (TTT {a.tag}): consecutive on TTT's traversal from {anchor}; graph link, {ev[1]} bp overlap; {rs} GraphAligner reads q>={a.minq} span the junction")
        else: R.append(f"   junction not proposed: {tag}; groups {g1}/{g2} (traversal {grp})")
print("\n".join(R)); print("## group_edits rows"); print("\n".join(G)); print("## order_links rows"); print("\n".join(O))
