#!/usr/bin/env python3
'''
Build pseudochromosomes from ordered chains. Sequence comes from unitigs only, so the AGP is at unitig level.

Placement: a chain stays in its pseudochromosome if it has informative Hi-C (within-chain density >= --minWithin)
or, having little, Hi-C still places it decisively without a link (position LLR >= --minLLR), or its telomere faces
outwards at a scaffold end (--telomeres), or a graph/read link connects it to such a chain. Other chains become labelled unplaced scaffolds: being in a group says little when
there are few pairs, and a default position would be a false claim.

Junctions:
  direct   graph link between facing chain ends: merged, overlap removed, no gap
  bridge   path through short chains chosen by choose_bridges.py: merged, no gap
  read     reads join the ends but no graph overlap: gap, left to read filling
  hic      Hi-C order only: gap of --gap Ns, type U

Names and direction are conventions only: names come from the table given with --names, and a pseudochromosome is
written in the direction of its name's reference when --flip lists it.

--tails adds one-sided graph paths from make_tails.py beside the chain end they link to, at a Hi-C gap or a scaffold end.

Usage: build_scaffolds.py --order order.tsv --agp chains.agp --links end_links.tsv --choices choices.tsv
         --stats chain_stats.tsv --names names.tsv --unitigs unitigs.fasta --out prefix [--extraUnplaced fasta]
'''
import argparse, collections
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in ("order", "agp", "links", "choices", "stats", "names", "unitigs", "out"): p.add_argument("--" + x, required=True)
p.add_argument("--minWithin", type=float, default=700); p.add_argument("--minLLR", type=float, default=500)
p.add_argument("--minOverlap", type=int, default=5000); p.add_argument("--minReads", type=int, default=3)
p.add_argument("--gap", type=int, default=100)
p.add_argument("--telomeres", help="telomeres.py output for the chains")
p.add_argument("--extraUnplaced", help="FASTA of further unplaced sequence, copied as is")
p.add_argument("--tails", help="make_tails.py output: end, group, path written outward from the end")
a = p.parse_args()
fl = {"+": "-", "-": "+"}

rows = collections.defaultdict(list); clen = {}
for l in open(a.agp):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    if s[4] == "W": rows[s[0]].append([s[5], int(s[6]), int(s[7]), s[8]])
    clen[s[0]] = max(clen.get(s[0], 0), int(s[2]))
within = {l.split("\t")[0]: float(l.split("\t")[2]) for l in open(a.stats) if l[0] != "#"}
names = {}; flips = set()
for l in open(a.names):
    if l[0] == "#": continue
    g, n, f = l.rstrip().split("\t"); names[g] = n
    if f == "flip": flips.add(g)

ov = {}; rd = {}
for l in open(a.links):
    if l[0] == "#": continue
    x, y, r, q, o = l.rstrip().split("\t"); k = frozenset((x, y)); ov[k] = int(o); rd[k] = int(r)
bridge = {}; alts = collections.Counter(); bridge_info = {}
for l in open(a.choices):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t")
    k = frozenset((s[0], s[1])); bridge[k] = (s[0], s[1], s[3].split(",")); bridge_info[k] = s
def junction(c1, o1, c2, o2):
    '''kind, and the path (oriented chains) between two chains in the scaffold'''
    e1 = c1 + o1; e2 = c2 + fl[o2]; k = frozenset((e1, e2))
    if ov.get(k, -1) >= a.minOverlap: return "direct", []
    if k in bridge:
        x, y, pth = bridge[k]
        if x == e1: return "bridge", [(z[:-1], z[-1]) for z in pth]
        return "bridge", [(z[:-1], fl[z[-1]]) for z in reversed(pth)]
    if rd.get(k, 0) >= a.minReads: return "read", []
    return "hic", []

order = collections.defaultdict(list); group_of = {}; status = {}
for l in open(a.order):
    if l[0] == "#": continue
    s = l.rstrip("\n").split("\t"); group_of[s[2]] = s[0]
    if s[1] == ".": status[s[2]] = "unplaced:too_few_pairs"; continue
    order[s[0]].append({"c": s[2], "o": s[3], "pllr": float(s[7])})

tel = {}
if a.telomeres:
    for l in open(a.telomeres):
        if l[0] != "#":
            r = l.rstrip().split("\t"); tel[r[0]] = r[4]
def tel_out(x, xs):
    call = tel.get(x["c"], "--")
    if x is xs[0] and (call[0] if x["o"] == "+" else call[1]) == "T": return True
    return x is xs[-1] and (call[1] if x["o"] == "+" else call[0]) == "T"
# Keep linked components that contain a chain placed with confidence
plan = {}
for g, xs in order.items():
    comp = []; cur = [xs[0]]
    for u, v in zip(xs, xs[1:]):
        if junction(u["c"], u["o"], v["c"], v["o"])[0] in ("direct", "bridge", "read"): cur.append(v)
        else: comp.append(cur); cur = [v]
    comp.append(cur)
    keep = []
    for cp in comp:
        ok = any(within.get(x["c"], 0) >= a.minWithin or (len(cp) == 1 and x["pllr"] >= a.minLLR) or tel_out(x, xs) for x in cp)
        for x in cp:
            status[x["c"]] = "placed" if ok else "unplaced:hic_not_decisive"
            if ok: keep.append(x)
    plan[g] = keep

U = {}; n = None
for l in open(a.unitigs):
    if l[0] == ">": n = l[1:].split()[0]; U[n] = []
    else: U[n].append(l.strip())
for k in U: U[k] = "".join(U[k])
tr = str.maketrans("ACGTNacgtn", "TGCANtgcan"); rc = lambda s: s.translate(tr)[::-1]

def oriented(c, o):
    if o == "+": return [r[:] for r in rows[c]]
    return [[r[0], r[1], r[2], fl[r[3]]] for r in reversed(rows[c])]
def trim(rs, k, what):
    '''remove k bp from the start of an oriented list of components'''
    while k > 0 and rs:
        r = rs[0]; n = r[2] - r[1] + 1
        if k >= n: rs.pop(0); k -= n; continue
        if r[3] == "+": r[1] += k
        else: r[2] -= k
        k = 0
    return rs

tails = {}
if a.tails:
    for l in open(a.tails):
        if l[0] != "#":
            r = l.rstrip("\n").split("\t"); tails[r[0]] = r[2].split(",")
tail_used = collections.Counter()
agp = open(a.out + ".agp", "w"); agp.write("##agp-version 2.1\n"); fa = open(a.out + ".fasta", "w")
jt = open(a.out + ".junctions.tsv", "w"); jt.write("#scaffold\tpos\tleft\tright\tkind\toverlap\tspanning_reads\tpath\tleft_join_reads\tright_join_reads\talternatives\n")
used_in_bridge = {}
out = {}
for g in sorted(plan, key=lambda g: names[g]):
    xs = plan[g]; name = names[g]; els = []
    seq = [(x["c"], x["o"]) for x in xs]
    if g in flips: seq = [(c, fl[o]) for c, o in reversed(seq)]
    comps = []; pos = 0; pieces = []; part = 0
    def add(rs):
        global pos, part
        for r in rs:
            m = r[2] - r[1] + 1; part += 1
            agp.write(f"{name}\t{pos+1}\t{pos+m}\t{part}\tW\t{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\n"); pos += m
            x = U[r[0]][r[1] - 1:r[2]]; pieces.append(x if r[3] == "+" else rc(x))
    def tail_out(end):
        '''scaffold runs out of this end: add its tail in the path's own direction'''
        prev = end
        for z in tails.get(end, []):
            c_, s_ = z[:-1], z[-1]; o_ = max(ov.get(frozenset((prev, c_ + fl[s_])), 0), 0)
            add(trim(oriented(c_, s_), o_, None)); prev = c_ + s_; used_in_bridge[c_] = name; tail_used[end] += 1
    def tail_in(end, c, o):
        '''scaffold runs into this end of (c, o): add its tail reversed, then the chain with the last overlap removed'''
        path = tails.get(end, [])
        if not path: add(oriented(c, o)); return
        rp = [(z[:-1], fl[z[-1]]) for z in reversed(path)]
        add(oriented(*rp[0])); used_in_bridge[rp[0][0]] = name
        for (c0, s0), (c_, s_) in zip(rp, rp[1:]):
            o_ = max(ov.get(frozenset((c0 + s0, c_ + fl[s_])), 0), 0)
            add(trim(oriented(c_, s_), o_, None)); used_in_bridge[c_] = name
        c0, s0 = rp[-1]; o_ = max(ov.get(frozenset((c0 + s0, end)), 0), 0)
        add(trim(oriented(c, o), o_, None)); tail_used[end] += len(path)
    tail_in(seq[0][0] + fl[seq[0][1]], *seq[0])
    for (c1, o1), (c2, o2) in zip(seq, seq[1:]):
        kind, pth = junction(c1, o1, c2, o2); e1 = c1 + o1; e2 = c2 + fl[o2]; k = frozenset((e1, e2))
        if kind == "direct":
            jt.write(f"{name}\t{pos}\t{c1}{o1}\t{c2}{o2}\tdirect\t{ov[k]}\t{rd.get(k, 0)}\t.\t.\t.\t.\n")
            add(trim(oriented(c2, o2), ov[k], k))
        elif kind == "bridge":
            prev = e1; info = bridge_info[k]
            jt.write(f"{name}\t{pos}\t{c1}{o1}\t{c2}{o2}\tbridge\t.\t.\t{','.join(z + s for z, s in pth)}\t{info[4]}\t{info[5]}\t{info[8]}\n")
            for z, s in pth:
                entry = z + fl[s]; o_ = max(ov.get(frozenset((prev, entry)), 0), 0)
                add(trim(oriented(z, s), o_, None)); prev = z + s; used_in_bridge[z] = name
            o_ = max(ov.get(frozenset((prev, e2)), 0), 0)
            add(trim(oriented(c2, o2), o_, None))
        else:
            if kind == "hic": tail_out(e1)
            part += 1
            agp.write(f"{name}\t{pos+1}\t{pos+a.gap}\t{part}\tU\t{a.gap}\tscaffold\tyes\tproximity_ligation\n")
            jt.write(f"{name}\t{pos}\t{c1}{o1}\t{c2}{o2}\t{kind}\t.\t{rd.get(k, 0)}\t.\t.\t.\t.\n")
            pos += a.gap; pieces.append("N" * a.gap)
            if kind == "hic": tail_in(e2, c2, o2)
            else: add(oriented(c2, o2))
    tail_out(seq[-1][0] + seq[-1][1])
    fa.write(f">{name}\n{''.join(pieces)}\n")
    out[name] = pos
fa.close(); agp.close(); jt.close()

placed = {x["c"] for xs in plan.values() for x in xs}
with open(a.out + ".unplaced.fasta", "w") as o, open(a.out + ".placement.tsv", "w") as t:
    t.write("#chain\tlength\tgroup\tname\tstatus\n")
    for c in clen:
        g = group_of.get(c, "unassigned")
        if c in placed: st = "placed"
        elif c in used_in_bridge: st = "placed:bridge"
        else: st = status.get(c, "unplaced:no_group")
        nm = names.get(g, ".") if st.startswith("placed") else (names.get(g, ".") + "?" if g in names else ".")
        if c in used_in_bridge: nm = used_in_bridge[c]
        t.write(f"{c}\t{clen[c]}\t{g}\t{nm}\t{st}\n")
        if not st.startswith("placed"):
            q = []
            for r in rows[c]:
                x = U[r[0]][r[1] - 1:r[2]]; q.append(x if r[3] == "+" else rc(x))
            o.write(f">unplaced_{c}\n{''.join(q)}\n")
    if a.extraUnplaced:
        for l in open(a.extraUnplaced): o.write(l)
if a.tails:
    miss = [e for e in tails if not tail_used[e]]
    print("tails placed:", len(tails) - len(miss), "not placed:", miss)
print("pseudochromosomes:", len(out), f"{sum(out.values())/1e6:.1f} Mb")
for n_ in sorted(out): print(f"  {n_}\t{out[n_]/1e6:.2f}")
