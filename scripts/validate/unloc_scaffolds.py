#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Unlocalized scaffolds: unplaced chains whose chromosome is known but whose position is not, written per haplotype as <scaffold>_unloc_NN
(the AGP/INSDC "unlocalized" class). Nothing in the scaffolds changes; the files sit beside the final ones, for the user to adopt or not (D13).
Ploidy 2 only.
Chromosome: every graph contact of the chain's unplaced component (--components) is to one chromosome's scaffolds.
Haplotype, within a chromosome: allele pairs (--alleles; >= --minAllele bp aligned and >= --minFrac of the shorter) must split.
  opposite_placed   an allele placed on one haplotype's scaffold sends the chain to the other
  arbitrary_pair    other allele sets are split heaviest pairs first, then the larger side goes to the haplotype with less unloc sequence:
                    the phase of the set is not known
  copy_2x           no allele and about 2x haplotype depth: written to both haplotypes
  arbitrary_single  no allele: the lighter haplotype
Left unplaced (unloc.skipped.tsv): repeats (> --maxRepeat repeated 31-mers), chains with alleles placed on both haplotypes, and, with --purge f,
the shorter of two chains on one haplotype that align over >= f of its length (overlapping fragments of one haplotype: these duplicate genes).
Measure the result with compleasm on each haplotype plus its unloc sequences (stage 15 does): adding the per-haplotype and unplaced
tables does not work, because compleasm drops weaker copies only when a better one is in the same run (on one 292 build that sum gave
D 3.4% where the run gave 2.0%).
Usage: unloc_scaffolds.py AGP UNPLACED_FASTA out_dir --hap DEPTH --components unplaced_components.txt --placement placement.tsv
  --repeat chain_repeat.tsv --depth depth.unitigs.tsv --chainsAgp chains.agp --alleles alleles.pairs.tsv [--purge F]'''
import argparse, ast, collections, os, re
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("agp"); p.add_argument("unplaced"); p.add_argument("out"); p.add_argument("--hap", type=float, required=True)
p.add_argument("--minAllele", type=int, default=50000); p.add_argument("--minFrac", type=float, default=0.1)
p.add_argument("--maxRepeat", type=float, default=0.5); p.add_argument("--copy2x", type=float, default=1.6)
p.add_argument("--purge", type=float, default=0.2)
for x in ("components", "placement", "repeat", "depth", "chainsAgp", "alleles"): p.add_argument("--" + x, required=True)
a = p.parse_args(); os.makedirs(a.out, exist_ok=True)
sample = os.path.basename(a.agp).split(".final")[0]
key = lambda s: re.sub(r"_h\d+$", "", s)                                 # scaffold -> chromosome key
def fasta(path):
    name, seq = None, []
    for l in open(path):
        if l[0] == ">":
            if name: yield name, "".join(seq)
            name, seq = l[1:].split()[0], []
        else: seq.append(l.strip())
    if name: yield name, "".join(seq)
def rows(path):
    for l in open(path):
        if l[0] != "#" and l.strip(): yield l.rstrip("\n").split("\t")
unpl = dict(fasta(a.unplaced))
U = {n[len("unplaced_"):]: n for n in unpl if n.startswith("unplaced_")}   # chain -> record name
L = {c: len(unpl[n]) for c, n in U.items()}
scaf = sorted({r[0] for r in rows(a.agp)})
haps_of = collections.defaultdict(list)
for s in scaf: haps_of[key(s)].append(s)
place = {}                                                                # placed chain (base name) -> scaffolds
for s in rows(a.placement):
    if s[4] == "placed": place.setdefault(s[0].split(".")[0], set()).add(s[3])
rep = {s[0]: float(s[2]) for s in rows(a.repeat)}
ud = {s[0]: float(s[2]) for s in rows(a.depth)}
dsum = collections.Counter(); dlen = collections.Counter()               # depth per chain, length-weighted over its unitigs
for s in rows(a.chainsAgp):
    if len(s) > 8 and s[4] not in "NU" and s[5] in ud:
        n = int(s[7]) - int(s[6]) + 1; dsum[s[0]] += ud[s[5]] * n; dlen[s[0]] += n
depth = {c: dsum[c] / dlen[c] for c in dlen}
chrom = {}; comp = None
for l in open(a.components):
    m = re.match(r"^\d+ chains, [\d.]+ Mb; graph contacts: (.*)$", l.strip())
    if m: comp = m.group(1); continue
    if l.strip().startswith("chains:") and comp and comp != "none":
        ks = {key(s) for s in ast.literal_eval(comp)}
        if len(ks) == 1:
            k = ks.pop()
            for c in l.split()[1:]:
                if c in U: chrom[c] = k
        comp = None
skip = {}
for c in list(chrom):
    if rep.get(c, 0) > a.maxRepeat: skip[c] = "repeat"; del chrom[c]
    elif len(haps_of[chrom[c]]) != 2: skip[c] = f"chromosome {chrom[c]} has {len(haps_of[chrom[c]])} scaffolds"; del chrom[c]
al = collections.defaultdict(dict); alAll = collections.defaultdict(dict)  # allele pairs: kept for phasing; any length for --purge
for s in rows(a.alleles):
    x, y, b = s[0], s[1], min(int(s[2]), int(s[3])); alAll[x][y] = alAll[y][x] = b
    lx = L.get(x) or dlen.get(x); ly = L.get(y) or dlen.get(y)
    if b >= a.minAllele and lx and ly and b >= a.minFrac * min(lx, ly): al[x][y] = b; al[y][x] = b
par, pty = {}, {}                                                         # signed union-find: parent, parity to parent (1 = other haplotype)
def find(x):
    if par.setdefault(x, x) == x: pty.setdefault(x, 0); return x
    p0 = par[x]; r = find(p0)
    if p0 != r: pty[x] ^= pty[p0]; par[x] = r
    return r
def parity(x): find(x); return 0 if par[x] == x else pty[x]
def union(x, y, diff):
    rx, ry = find(x), find(y); px, py = parity(x), parity(y)
    if rx == ry: return (px ^ py) == diff
    par[ry] = rx; pty[ry] = px ^ py ^ diff; return True
out_rows, conflicts = [], collections.Counter()
for k in sorted(set(chrom.values())):
    h1, h2 = sorted(haps_of[k]); H = {h1: "H1", h2: "H2"}
    par.clear(); pty.clear(); union("H1", "H2", 1)
    cs = sorted((c for c in chrom if chrom[c] == k), key=lambda c: -L[c])
    edges = []
    for c in cs:
        for o, b in al[c].items():
            if chrom.get(o) == k: edges.append((b, c, o))
            elif o in place: edges += [(10**12 + b, c, H[sc]) for sc in place[o] if sc in H]   # placed alleles bind first
    for b, x, y in sorted(edges, reverse=True):
        if not union(x, y, 1): conflicts[k] += 1
    for c in cs:
        if len({sc for o in al[c] if o in place for sc in place[o] if sc in H}) == 2: skip[c] = "allele placed on both haplotypes"
    cs = [c for c in cs if c not in skip]
    hbp = collections.Counter(); tied = find("H1"); sets = collections.defaultdict(list); side = {}
    for c in cs: sets[find(c)].append(c)
    for r, mem in sorted(sets.items(), key=lambda t: (-sum(L[c] for c in t[1]), t[1][0])):
        if r == tied:
            for c in mem: side[c] = (h1 if parity(c) == parity("H1") else h2, "opposite_placed")
            continue
        if len(mem) == 1 and not any(chrom.get(o) == k for o in al[mem[0]]):
            c = mem[0]
            if depth.get(c, 0) >= a.copy2x * a.hap: side[c] = ("both", "copy_2x"); continue
            t = min((h1, h2), key=lambda h: (hbp[h], h)); side[c] = (t, "arbitrary_single"); hbp[t] += L[c]; continue
        a0 = sum(L[c] for c in mem if parity(c) == 0); a1 = sum(L[c] for c in mem if parity(c) == 1)
        light, heavy = sorted((h1, h2), key=lambda h: (hbp[h], h)); big = 0 if a0 >= a1 else 1
        for c in mem:
            h = light if parity(c) == big else heavy; side[c] = (h, "arbitrary_pair"); hbp[h] += L[c]
    if a.purge > 0:                                                       # overlapping fragments of one haplotype: keep the longer
        for c in sorted(cs, key=lambda c: (L[c], c)):
            if side[c][0] == "both": continue
            for o, b in sorted(alAll[c].items()):
                if o != c and o in side and o not in skip and side[o][0] == side[c][0] and L[o] >= L[c] and b >= a.purge * L[c]:
                    skip[c] = f"overlaps {o} on the same haplotype ({b} bp aligned)"; break
        cs = [c for c in cs if c not in skip]
    for c in cs:
        h, ph = side[c]
        for hh in ((h1, h2) if h == "both" else (h,)): out_rows.append((hh, c, L[c], ph, depth.get(c, 0) / a.hap))
n_by = collections.Counter()
with open(f"{a.out}/{sample}.unloc.fasta", "w") as fo, open(f"{a.out}/{sample}.unloc.tsv", "w") as ft:
    ft.write("#name\tchain\tlength\tscaffold\tphase\tdepth_x_hap\n")
    for h, c, n, ph, dx in sorted(out_rows, key=lambda r: (r[0], -r[2], r[1])):
        n_by[h] += 1; nm = f"{h}_unloc_{n_by[h]:02d}"; s = unpl[U[c]]
        fo.write(f">{nm} {c}\n" + "".join(s[i:i + 80] + "\n" for i in range(0, len(s), 80)))
        ft.write(f"{nm}\t{c}\t{n}\t{h}\t{ph}\t{dx:.2f}\n")
used = {r[1] for r in out_rows}
with open(f"{a.out}/{sample}.unplaced_rest.fasta", "w") as fo:
    for nm, s in unpl.items():
        if nm[len("unplaced_"):] not in used: fo.write(f">{nm}\n" + "".join(s[i:i + 80] + "\n" for i in range(0, len(s), 80)))
with open(f"{a.out}/{sample}.unloc.skipped.tsv", "w") as fo:
    fo.write("#chain\tlength\treason\n")
    for c, r in sorted(skip.items()): fo.write(f"{c}\t{L.get(c, 0)}\t{r}\n")
ph = collections.Counter(); bp = collections.Counter()
for r in out_rows: ph[r[3]] += 1; bp[r[3]] += r[2]
print(f"Unlocalized: {len(used)} chains, {sum(L[c] for c in used) / 1e6:.2f} Mb as {len(out_rows)} sequences ({a.out}/{sample}.unloc.fasta); "
      f"by phase: {', '.join(f'{k} {ph[k]} ({bp[k] / 1e6:.1f} Mb)' for k in sorted(ph))}; left unplaced {len(skip)} (unloc.skipped.tsv); "
      f"allele conflicts {dict(conflicts) or 0}. Phase is not known for arbitrary_* sequences.")
for h in sorted(n_by): print(f"  {h}: {n_by[h]} sequences, {sum(r[2] for r in out_rows if r[0] == h) / 1e6:.2f} Mb")
