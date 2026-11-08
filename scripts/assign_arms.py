#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Diploid bridge choice. chain_bridges.py finds paths per haplotype copy; a bubble between shared chains then appears once
per haplotype, and choose_bridges.py alone would give both haplotypes the arm with more reads, losing the other allele.

  prepare: collapse copy-level bridges to base chains (same junction, same paths) for choose_bridges.py
  assign:  read its ranking and write one choice per copy-level junction. Parallel junctions (the same base junction in
           two groups) take different paths when a second exists. Which haplotype takes which is settled by Hi-C links
           from the arm chains to each group's specific chains when >= --minLinks and >= --ratio; otherwise it is recorded
           as arbitrary. Other junctions take their best path whose unassigned chains are not already used by the other
           haplotype; forced re-use is reported.

Usage: assign_arms.py prepare graph.bridges.tsv base.bridges.tsv
       assign_arms.py assign graph.bridges.tsv choices_base.tsv groups.tsv links.tsv choices.tsv report.tsv
'''
import argparse, collections
fl = {"+": "-", "-": "+"}
base = lambda c: c.split(".g")[0]
bend = lambda e: base(e[:-1]) + e[-1]
def rev(path): return [z[:-1] + fl[z[-1]] for z in reversed(path)]
def norm(x, y, path):
    '''junction key and path written from the lower base end'''
    bx, by = bend(x), bend(y)
    return (bx, by, path) if bx <= by else (by, bx, rev(path))
def read_bridges(f):
    for l in open(f):
        if l[0] == "#": continue
        r = l.rstrip("\n").split("\t")
        if r[2] == "yes": yield r[0], r[1], r[3].split(","), int(r[4]), r[5]
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("mode"); ap.add_argument("files", nargs="+")
ap.add_argument("--minLinks", type=int, default=5); ap.add_argument("--ratio", type=float, default=2.0)
a = ap.parse_args()
if a.mode == "prepare":
    seen = {}
    for x, y, path, bp, mo in read_bridges(a.files[0]):
        bx, by, pth = norm(x, y, path); seen[(bx, by, tuple(pth))] = (bp, mo)
    with open(a.files[1], "w") as o:
        o.write("#end_a\tend_b\tsame_group\tpath\tpath_bp\tmin_overlap\n")
        for (bx, by, pth), (bp, mo) in sorted(seen.items()):
            o.write(f"{bx}\t{by}\tyes\t{','.join(pth)}\t{bp}\t{mo}\n")
    print(len(seen), "base candidate paths,", len({(k[0], k[1]) for k in seen}), "base junctions"); raise SystemExit
bridges, choices, groups, links, out, report = a.files
rank = {}                                       # base junction -> [(path tuple, left, right, full)], best first
for l in open(choices):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    k = (s[0], s[1]); v = [(tuple(s[3].split(",")), int(s[4]), int(s[5]), int(s[6]))]
    if s[8] != ".":
        for alt in s[8].split(";"):
            pth, c = alt.rsplit(":", 1); c = c.split("/"); v.append((tuple(pth.split(",")), int(c[0]), int(c[1]), int(c[2])))
    rank[k] = v
gr = {}; L = {}
for l in open(groups):
    if l[0] != "#":
        c, n, g = l.rstrip().split("\t"); gr[c] = g; L[c] = int(n)
members = collections.defaultdict(set)
for c, g in gr.items():
    if g != "unassigned" and ".g" not in c: members[g].add(c)
lk = collections.defaultdict(dict)
for l in open(links):
    x, y, n = l.rstrip().split("\t")
    if y != ".": lk[x][y] = int(n); lk[y][x] = int(n)
J = collections.defaultdict(list)               # base junction -> copy junctions
for x, y, path, bp, mo in read_bridges(bridges):
    bx, by, pth = norm(x, y, path); k = (bx, by)
    if (x, y) not in [q[:2] for q in J[k]]: J[k].append((x, y, gr[x[:-1]]))
def support(k, pth):
    for r in rank[k]:
        if r[0] == pth: return r
def hic(path, g):
    return sum(lk[z[:-1]].get(m, 0) for z in path for m in members[g])
used = collections.defaultdict(set); rows = []; rep = []; stats = collections.Counter()
order = sorted(J, key=lambda k: (-len(J[k]), -rank[k][0][3]))
for k in order:
    copies = J[k]; R = rank[k]
    if len(copies) == 2 and len(R) >= 2:
        (x1, y1, g1), (x2, y2, g2) = copies
        p, q = R[0][0], R[1][0]
        a1 = hic(p, g1) + hic(q, g2); a2 = hic(p, g2) + hic(q, g1)
        if max(a1, a2) >= a.minLinks and max(a1, a2) >= a.ratio * max(min(a1, a2), 1): how = "hic"
        else: how = "arbitrary"
        if a2 > a1 and how == "hic": p, q = q, p
        pick = [(copies[0], p), (copies[1], q)]; stats["parallel_" + how] += 1
    elif len(copies) == 2:
        pick = [(copies[0], R[0][0]), (copies[1], R[0][0])]; how = "one_path"; stats["parallel_one_path"] += 1
    else:
        c = copies[0]; other_g = None; choice = None
        for r in R:
            if not any(used[z[:-1]] - {c[2]} for z in r[0]): choice = r[0]; break
        if choice is None: choice = R[0][0]; how = "reused"; stats["single_reused"] += 1
        else: how = "best"; stats["single_best"] += 1
        pick = [(c, choice)]
    for (x, y, g), pth in pick:
        for z in pth: used[z[:-1]].add(g)
        r = support(k, pth); bx = bend(x)
        P = list(pth) if bx == k[0] else rev(list(pth))
        rows.append(f"{x}\t{y}\tyes\t{','.join(P)}\t{r[1]}\t{r[2]}\t{r[3]}\t{sum(L[z[:-1]] for z in pth)}\t{how}:" + (";".join(f"{','.join(t[0])}:{t[1]}/{t[2]}/{t[3]}" for t in R if t[0] != pth) or ".") + "\n")
        rep.append(f"{g}\t{x}\t{y}\t{how}\t{','.join(P)}\t{r[1]}\t{r[2]}\t{r[3]}\n")
with open(out, "w") as o:
    o.write("#end_a\tend_b\tchosen\tpath\tleft_reads\tright_reads\tfull_reads\tpath_bp\talternatives\n"); o.writelines(rows)
with open(report, "w") as o:
    o.write("#group\tend_a\tend_b\thow\tpath\tleft_reads\tright_reads\tfull_reads\n"); o.writelines(rep)
print(dict(stats))
