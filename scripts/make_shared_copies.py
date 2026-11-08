#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Give each shared (collapsed) chain one copy per haplotype, so later steps can keep their one-group-per-chain model.
A copy is named <chain>.<group>. Every per-chain input for ordering and scaffolding is rewritten for the copies:

  groups    each copy in its group; other rows unchanged
  agp, fai  copies replace the chain; the same unitig range appears once per haplotype
  pairs     shared-specific pair: goes to the copy in the specific chain's group, so each copy sees its own haplotype
            shared-shared and within-chain pair: goes to every copy, since nothing says which haplotype it came from
            shared-unassigned pair: dropped
  links     shared-specific link: to the copy in the specific chain's group
            shared-shared link, same chromosome: between same-group copies
            shared-unassigned link: from every copy, so bridges can be found from each haplotype
  stats, telomeres, repeat: each copy repeats the chain's row

A shared chain with no Hi-C assignment takes its chromosome from graph links (>= --minOverlap) to assigned shared chains or
group members, if they agree. Otherwise it stays as a single unassigned chain.
Usage: make_shared_copies.py groups_prefix chains.agp chains.fai pairs.tsv.gz end_links.tsv chain_stats.tsv telomeres.tsv out_dir
'''
import argparse, collections, gzip, os
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in "prefix agp fai pairs links stats telomeres out".split(): p.add_argument(x)
p.add_argument("--minOverlap", type=int, default=5000)
a = p.parse_args(); os.makedirs(a.out, exist_ok=True)
L = {l.split("\t")[0]: int(l.split("\t")[1]) for l in open(a.fai)}
gr = {}
for l in open(a.prefix + ".tsv"):
    if l[0] != "#":
        c, n, g = l.rstrip().split("\t"); gr[c] = g
chrom_of_group = dict(l.split() for l in open(a.prefix + ".chroms.tsv"))
groups_of_chrom = collections.defaultdict(list)
for g, c in chrom_of_group.items(): groups_of_chrom[c].append(g)
shared = {}                                     # chain -> chromosome
sh_all = set()
for l in open(a.prefix + ".shared.tsv"):
    if l[0] == "#": continue
    r = l.rstrip().split("\t"); sh_all.add(r[0])
    if r[2] != "unassigned": shared[r[0]] = r[2]
LK = []
for l in open(a.links):
    if l[0] == "#": continue
    r = l.rstrip("\n").split("\t"); LK.append(r)
# graph-assign the rest, repeating because neighbours can be assigned in turn
def chrom_of(c):
    if c in shared: return shared[c]
    return chrom_of_group.get(gr.get(c, ""), None)
changed = True
while changed:
    changed = False
    votes = collections.defaultdict(collections.Counter)
    for x, y, rd, q, ov in LK:
        if int(ov) < a.minOverlap: continue
        u, v = x[:-1], y[:-1]
        for s, t in ((u, v), (v, u)):
            if s in sh_all and s not in shared and chrom_of(t): votes[s][chrom_of(t)] += 1
    for c, v in votes.items():
        if len(v) == 1:
            shared[c] = next(iter(v)); changed = True; print(f"# {c} assigned to {shared[c]} by graph links")
for c in sh_all - set(shared): print(f"# {c} ({L[c]} bp) stays unassigned")
copies = {c: [f"{c}.{g}" for g in groups_of_chrom[k]] for c, k in shared.items()}
def copy_in(c, g): return f"{c}.{g}"

with open(f"{a.out}/groups.tsv", "w") as o:
    o.write("#chain\tlength\tgroup\n")
    for c, g in gr.items():
        if c in shared: continue
        o.write(f"{c}\t{L[c]}\t{'unassigned' if g == 'shared' else g}\n")
    for c, k in shared.items():
        for g in groups_of_chrom[k]: o.write(f"{copy_in(c, g)}\t{L[c]}\t{g}\n")
with open(f"{a.out}/chains.agp", "w") as o:
    for l in open(a.agp):
        if l[0] == "#": o.write(l); continue
        r = l.split("\t", 1)
        if r[0] in shared:
            for n in copies[r[0]]: o.write(f"{n}\t{r[1]}")
        else: o.write(l)
with open(f"{a.out}/chains.fasta.fai", "w") as o:
    for c, n in L.items():
        for x in copies.get(c, [c]): o.write(f"{x}\t{n}\t0\t0\t0\n")
def group_of(c): return gr[c] if gr.get(c) in chrom_of_group else None
n_out = collections.Counter()
with gzip.open(a.pairs, "rt") as fh, gzip.open(f"{a.out}/pairs.tsv.gz", "wt", compresslevel=1) as o:
    for l in fh:
        x, u, y, v = l.rstrip("\n").split("\t"); sx, sy = x in shared, y in shared
        if not sx and not sy: o.write(l); n_out["plain"] += 1; continue
        if sx and sy:
            if x == y or shared[x] == shared[y]:
                for g in groups_of_chrom[shared[x]]: o.write(f"{copy_in(x, g)}\t{u}\t{copy_in(y, g)}\t{v}\n")
                n_out["shared_both"] += 1
            else: n_out["shared_other_chrom"] += 1
            continue
        s, t = (x, y) if sx else (y, x); g = group_of(t)
        if g is None or chrom_of_group[g] != shared[s]: n_out["dropped"] += 1; continue
        if sx: o.write(f"{copy_in(x, g)}\t{u}\t{y}\t{v}\n")
        else: o.write(f"{x}\t{u}\t{copy_in(y, g)}\t{v}\n")
        n_out["shared_to_copy"] += 1
print("# pairs", dict(n_out))
with open(f"{a.out}/end_links.tsv", "w") as o:
    o.write("#end_a\tend_b\treads\tq_gap_median\tgraph_overlap\n")
    for x, y, rd, q, ov in LK:
        u, eu, v, ev = x[:-1], x[-1], y[:-1], y[-1]; su, sv = u in shared, v in shared
        if not su and not sv: o.write(f"{x}\t{y}\t{rd}\t{q}\t{ov}\n"); continue
        if su and sv:
            if shared[u] != shared[v]: continue
            for g in groups_of_chrom[shared[u]]: o.write(f"{copy_in(u, g)}{eu}\t{copy_in(v, g)}{ev}\t{rd}\t{q}\t{ov}\n")
            continue
        s, es, t, et = (u, eu, v, ev) if su else (v, ev, u, eu); g = group_of(t)
        if g is None:
            for g2 in groups_of_chrom[shared[s]]: o.write(f"{copy_in(s, g2)}{es}\t{t}{et}\t{rd}\t{q}\t{ov}\n")
        elif chrom_of_group[g] == shared[s]: o.write(f"{copy_in(s, g)}{es}\t{t}{et}\t{rd}\t{q}\t{ov}\n")
for src, dst in ((a.stats, "chain_stats.tsv"), (a.telomeres, "telomeres.tsv")):
    with open(f"{a.out}/{dst}", "w") as o:
        for l in open(src):
            if l[0] == "#": o.write(l); continue
            r = l.split("\t", 1)
            for n in copies.get(r[0], [r[0]]): o.write(f"{n}\t{r[1]}")
print(f"shared chains copied: {len(shared)}, {sum(L[c] for c in shared)/1e6:.1f} Mb; left unassigned: {len(sh_all - set(shared))}")
