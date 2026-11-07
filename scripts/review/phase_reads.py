#!/usr/bin/env python3
'''Do raw reads phase arms labelled arbitrary? No new mapping. For each such junction, finds the arm on the scaffold and the nearest haplotype-specific chain on each side, then counts reads linking them.
cis: one alignment (MAPQ >= 20) overlaps the arm and runs at least minAnchor into specific sequence of this scaffold. trans: an alignment overlaps the arm here, and another alignment of the same read lies in specific sequence of a homolog.
confirmed: cis >= 2 and at least 5 x trans. conflict: trans >= 2 and over cis. mixed: other cases with at least 2 of either. thin or no_reads: reads cannot say. Most arms can lie too far from specific sequence for reads to reach.
Hi-C cannot do this (see the D15 guide). A conflict is a prompt, and no arm has been swapped on this evidence.
Usage: phase_reads.py graph.agp chains.agp junctions.tsv groups.tsv out.tsv --paf a.paf.gz [b.paf.gz] [--suffix _h1 _h2]'''
import argparse, collections, gzip, re
p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
for x in "agp chains junctions groups out".split(): p.add_argument(x)
p.add_argument("--paf", nargs="+", required=True); p.add_argument("--suffix", nargs="*", default=[]); p.add_argument("--minAnchor", type=int, default=2000); p.add_argument("--minQ", type=int, default=20); a = p.parse_args()
U2C = {}
for l in open(a.chains):
    if l[0] == "#": continue
    f = l.rstrip("\n").split("\t")
    if f[4] == "W": U2C.setdefault(f[5], []).append(f[0])
G = {l.split("\t")[0]: l.rstrip("\n").split("\t")[2] for l in open(a.groups) if l[0] != "#"}
# Specific: an assigned chain with no .gNN suffix; a unitig shared with other chains still counts only if it belongs to such a chain
def spec(u):
    return any(re.match(r"g\d+$", G.get(c, "")) and ".g" not in c and not c.startswith(("chain8", "chain9")) for c in U2C.get(u, []))
S = collections.defaultdict(list)
for l in open(a.agp):
    if l[0] == "#": continue
    f = l.rstrip("\n").split("\t")
    if f[4] == "W": S[f[0]].append((int(f[1]) - 1, int(f[2]), f[5], spec(f[5])))
for k in S: S[k].sort()
J = []
for l in open(a.junctions):
    if l[0] == "#": continue
    f = l.rstrip("\n").split("\t")
    if len(f) > 10 and f[10].startswith("arbitrary") and f[4] == "bridge":
        J.append(dict(s=f[0], pos=int(f[1]), left=f[2], right=f[3], path=f[7]))
def intervals(s):
    '''merged specific intervals of a scaffold'''
    out = []
    for x0, x1, u, sp in S[s]:
        if sp:
            if out and x0 <= out[-1][1] + 1000: out[-1][1] = max(out[-1][1], x1)
            else: out.append([x0, x1])
    return out
I = {s: intervals(s) for s in S}
stem = lambda t: t.rsplit("_h", 1)[0]
def inspec(t, x0, x1, min_bp):
    return any(min(x1, b) - max(x0, a_) >= min_bp for a_, b in I.get(t, []))
for j in J:
    v = S[j["s"]]; ch = {c.rstrip("+-") for c in j["path"].split(",")}
    us = [x for x in v if set(U2C.get(x[2], [])) & ch and abs((x[0] + x[1]) / 2 - j["pos"]) < 2e6]
    # arm: path unitigs nearest the junction
    if us:
        us.sort(key=lambda x: abs((x[0] + x[1]) // 2 - j["pos"])); arm0 = min(x[0] for x in us); arm1 = max(x[1] for x in us)
    else: arm0, arm1 = j["pos"], j["pos"] + 1
    j.update(arm0=arm0, arm1=arm1)
    L = [b for a_, b in I[j["s"]] if b <= arm0 + 1000]; R = [a_ for a_, b in I[j["s"]] if a_ >= arm1 - 1000]
    j["lspec"] = max(L) if L else None; j["rspec"] = min(R) if R else None
    j.update(cisL=set(), cisR=set(), through=set(), trans=set(), arm_reads=set())
by = collections.defaultdict(list)
for j in J: by[j["s"]].append(j)
def scan(f, suf):
    reads = collections.defaultdict(list)
    for l in gzip.open(f, "rt"):
        x = l.split("\t", 12)
        if int(x[11]) >= a.minQ: reads[x[0]].append((x[5], int(x[7]), int(x[8])))
    for r, al in reads.items():
        for t, x0, x1 in al:
            if suf and not t.endswith(suf): continue
            for j in by.get(t, []):
                if min(x1, j["arm1"]) - max(x0, j["arm0"]) < a.minAnchor: continue   # must overlap the arm, not cover it all
                j["arm_reads"].add(r)
                cl = j["lspec"] is not None and x0 <= j["lspec"] - a.minAnchor and x1 >= j["arm0"] + a.minAnchor
                cr = j["rspec"] is not None and x1 >= j["rspec"] + a.minAnchor and x0 <= j["arm1"] - a.minAnchor
                if cl: j["cisL"].add(r)
                if cr: j["cisR"].add(r)
                if cl and cr: j["through"].add(r)
                for t2, y0, y1 in al:
                    if t2 != t and stem(t2) == stem(t) and inspec(t2, y0, y1, a.minAnchor): j["trans"].add(r)
if a.suffix:
    for f, s in zip(a.paf, a.suffix): scan(f, s)
else: scan(a.paf[0], "")
n = len(J); cl = collections.Counter()
with open(a.out, "w") as o:
    o.write("scaffold\tpos\tarm_start\tarm_end\tarm_bp\tdist_left_specific\tdist_right_specific\tarm_reads\tcis_left\tcis_right\tthrough\ttrans\tclass\n")
    for j in J:
        d0 = j["arm0"] - j["lspec"] if j["lspec"] is not None else -1; d1 = j["rspec"] - j["arm1"] if j["rspec"] is not None else -1
        c, t = len(j["cisL"] | j["cisR"]), len(j["trans"])
        k = "no_reads" if not c and not t else "thin" if c < 2 and t < 2 else "confirmed" if c >= 2 and c >= 5 * t else "conflict" if t >= 2 and t > c else "mixed"
        cl[k] += 1
        o.write(f"{j['s']}\t{j['pos']}\t{j['arm0']}\t{j['arm1']}\t{j['arm1'] - j['arm0']}\t{d0}\t{d1}\t{len(j['arm_reads'])}\t{len(j['cisL'])}\t{len(j['cisR'])}\t{len(j['through'])}\t{t}\t{k}\n")
print(f"Arbitrary bridge junctions {n}: " + ", ".join(f"{k} {v}" for k, v in cl.most_common()) + ".")
