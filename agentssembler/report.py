"""One short report per stage: the numbers to look at, alarms, and what to decide. Reports describe; they never decide."""
import collections, glob, gzip, os, re, statistics as st
from .stages import BY_ID, GUIDES

def rows(p, comment="#"):
    if not os.path.exists(p): return
    for l in open(p):
        if l.strip() and not l.startswith(comment): yield l.rstrip("\n").split("\t")
def Mb(x): return f"{x / 1e6:.1f} Mb"
def fai(p): return {r[0]: int(r[1]) for r in rows(p)}
def n50(v):
    v = sorted(v, reverse=True); t = sum(v) / 2; c = 0
    for x in v:
        c += x
        if c >= t: return x
    return 0

def r_inputs(w, c, s):
    L = fai(f"{w}/00_inputs/unitigs.fasta.fai"); return [f"Unitigs: {len(L)}, {Mb(sum(L.values()))}, N50 {Mb(n50(L.values()))}", f"Chromosomes {c['chromosomes']}, ploidy {c['ploidy']}."], []
def r_depth(w, c, s):
    L = fai(f"{w}/00_inputs/unitigs.fasta.fai")
    def hist(title, pairs):
        H = collections.Counter()
        for d, n in pairs: H[min(int(d // 4) * 4, 200)] += n
        t = sum(H.values()) or 1; o = [f"{title} ({Mb(t)}), percentage per 4x bin:", "```"]
        o += [f"{k:4d}x {100 * H[k] / t:5.1f}% {'#' * int(100 * H[k] / t)}" for k in sorted(H) if H[k] / t >= 0.002]
        return o + ["```"], max(H, key=H.get) if H else 0
    C = [(r[0], float(r[2]), float(r[6])) for r in rows(f"{w}/04_copynumber/class_depth.tsv")]
    out = []
    o, _ = hist("Primary depth, all unitigs", [(p, 10000) for _, p, _c in C]); out += o
    o, peak = hist("Primary depth, unitigs of 1 Mb or more", [(p, 10000) for u, p, _c in C if L.get(u, 0) >= 1000000]); out += o
    o, _ = hist("Class depth, adding near-identical copies, all unitigs", [(k, 10000) for _, _p, k in C]); out += o
    # Everything below uses the main peak in long unitigs, and only describes. The peak is a candidate for hap_depth, not a decision.
    h = peak + 2 if peak else 0
    by = collections.defaultdict(list)
    for u, p, k in C: by[u].append((p, k))
    M = {u: (st.median(x[0] for x in v), st.median(x[1] for x in v)) for u, v in by.items()}
    deg = collections.Counter()
    if os.path.exists(f"{w}/00_inputs/unitigs.gfa"):
        for l in open(f"{w}/00_inputs/unitigs.gfa"):
            if l[0] == "L": x = l.split("\t"); deg[(x[1], x[2])] += 1
    if h:
        twin = sum(L[u] for u, (p, k) in M.items() if p < 0.64 * h and 0.64 * h <= k < 1.54 * h)
        low = sum(L[u] for u, (p, k) in M.items() if k < 0.64 * h)
        two = sorted(((L[u], u, p) for u, (p, k) in M.items() if 1.54 * h <= p <= 2.68 * h and k <= 2.68 * h), reverse=True)
        out += [f"Taking {h}x as a candidate haplotype depth (main bin in long unitigs; check it):",
                f"- Below 0.64x alone but 1x as a class, one haplotype represented more than once: {Mb(twin)}.",
                f"- Below 0.64x even as a class, truly low depth: {Mb(low)}.",
                f"- About 2x in both primary and class depth, candidates for collapse: {len(two)} unitigs, {Mb(sum(x[0] for x in two))}, of which graph degree at least 2 at both ends: {Mb(sum(x[0] for x in two if deg[(x[1], '+')] >= 2 and deg[(x[1], '-')] >= 2))}."]
        out += [f"  - {u} {Mb(n)} at {p:.0f}x, degree {deg[(u, '+')]},{deg[(u, '-')]}" for n, u, p in two[:8]]
    out.append("Set hap_depth, repeat_depth (about 4.4 x hap_depth) and collapsed. See D01.")
    return out, []
def r_chains(w, c, s):
    out = []; al = []
    for d in ("06_chains", "07_chains2", "08_final_chains"):
        L = fai(f"{w}/{d}/chains.fasta.fai")
        if L: out.append(f"{d}: {len(L)} chains, {Mb(sum(L.values()))}, N50 {Mb(n50(L.values()))}")
    a, b = (sum(fai(f"{w}/{d}/chains.fasta.fai").values()) for d in ("07_chains2", "08_final_chains"))
    if a and b and abs(a - b) > 0.05 * max(a, b): al.append(f"07_chains2 ({Mb(a)}) and 08_final_chains ({Mb(b)}) should hold nearly the same sequence: join verification may have run on incomplete chains.")
    V = list(rows(f"{w}/07_chains2/joins.verify.tsv")); low = [r for r in V if int(r[5]) < 10]; br = [r for r in V if int(r[5]) < 3]
    out.append(f"Joins {len(V)}; under 3 reads {len(br)} (broken by rule); under 10 reads {len(low)} (review in joins.verify.tsv).")
    if os.path.exists(f"{w}/07_chains2/restored_arms.tsv"): out.append(f"Restored arms: {sum(1 for _ in rows(f'{w}/07_chains2/restored_arms.tsv'))}.")
    return out, al
def r_hic(w, c, s):
    p = f"{w}/10_groups/chain_breaks.tsv"; V = list(rows(p)); inv = lambda x: x[:-1] + ("-" if x[-1] == "+" else "+")
    J = {}
    for r in rows(f"{w}/07_chains2/joins.verify.tsv"): J[(r[2], r[3])] = r; J[(inv(r[3]), inv(r[2]))] = r
    def cand(r):
        j = J.get((r[3], r[4]));
        if not j: return None
        return f"- {r[1]} at {int(r[2]):,}: partner_cos {r[7]}, trans {r[8]}/{r[9]}, cross_ratio {r[0]}; reads {j[5]}, at MAPQ>=20 {j[6]}, overlap {int(j[4]):,}"
    both = [x for r in V if min(int(r[8]), int(r[9])) >= 200 and float(r[7]) < 0.8 and J.get((r[3], r[4])) and int(J[(r[3], r[4])][6]) < 10 for x in [cand(r)]]
    hic = [x for r in V if min(int(r[8]), int(r[9])) >= 200 and float(r[7]) < 0.6 for x in [cand(r)] if x]
    out = [f"Candidate positions in {p}: {len(V)}. A break needs two signals: different Hi-C partners each side and weak reads across the join."]
    out += [f"Both signals (trans >= 200 on both sides, partner_cos < 0.8, under 10 reads at MAPQ>=20): {len(both)}."] + both
    out += [f"Different partners alone (partner_cos < 0.6, enough trans), usually overruled by reads: {len(hic)}."] + hic
    out.append("This table only lists. Weak reads with little Hi-C on one side are mostly near-identical sequence, not chimeras.")
    return out, []
def r_evidence(w, c, s):
    L = fai(f"{w}/12_broken/chains.fasta.fai"); K = collections.Counter(r[6] for r in rows(f"{w}/11_ont_chains/chain_class.tsv")); B = collections.Counter()
    for r in rows(f"{w}/11_ont_chains/chain_class.tsv"): B[r[6]] += int(r[1])
    T = sum(1 for r in rows(f"{w}/16_validate/telomeres.chains.tsv") if "T" in r[4])
    return [f"Chains: {len(L)}, {Mb(sum(L.values()))}, N50 {Mb(n50(L.values()))}", "Depth class: " + ", ".join(f"{k} {K[k]} ({Mb(B[k])})" for k in K), f"Chains with a telomere call: {T}."], []
def r_group(w, c, s):
    G = collections.Counter(); B = collections.Counter(); want = c["chromosomes"] * c["ploidy"]
    for r in rows(f"{w}/10_groups/groups.tsv"): G[r[2]] += 1; B[r[2]] += int(r[1])
    gs = sorted(k for k in G if re.match(r"g\d+$", k)); al = []
    if len(gs) != want: al.append(f"{len(gs)} groups, {want} expected. See D05. Do not force the count.")
    sizes = [B[k] for k in gs]; out = [f"Groups: {len(gs)} of {want}; assigned {Mb(sum(sizes))}; unassigned {Mb(B['unassigned'])}; shared {Mb(B['shared'])}.", "Size by group: " + ", ".join(f"{k} {B[k] / 1e6:.1f}" for k in gs)]
    if sizes and max(sizes) > 1.4 * st.median(sizes): al.append("A group is more than 1.4 x the median size: look for satellite chains stacked on it (D06).")
    T = collections.Counter(); tel = {r[0]: r[4] for r in rows(f"{w}/16_validate/telomeres.chains.tsv")}
    for r in rows(f"{w}/10_groups/groups.tsv"):
        if re.match(r"g\d+$", r[2]): T[r[2]] += tel.get(r[0], "--").count("T")
    over = [k for k, v in T.items() if v > 2]
    if over: al.append(f"More than 2 telomeric ends in {over}: a false member or two chromosomes merged (D06, D07).")
    if T: out.append("Telomeric ends by group: " + ", ".join(f"{k} {T[k]}" for k in gs) + ".")
    Un = [r[0] for r in rows(f"{w}/10_groups/groups.tsv") if r[2] == "unassigned"]
    ut = [f"{u} ({tel[u]})" for u in Un if "T" in tel.get(u, "--")]
    if ut and any(T[k] < 2 for k in gs): al.append(f"Groups with under 2 telomeric ends, and unassigned chains with a telomere: {', '.join(ut)}. Test them against those groups (D06).")
    # Long overlaps pull a hub chain's median depth down, so its class can read 1x. The middle third is the fairer figure.
    h = float(s.get("hap_depth", 0)) if isinstance(s, dict) else 0
    if h:
        W = collections.defaultdict(list)
        for r in rows(f"{w}/11_ont_chains/depth.windows.tsv"): W[r[0]].append(float(r[3]))
        hub = [f"{u} ({st.median(W[u][len(W[u]) // 3:2 * len(W[u]) // 3]):.0f}x)" for u in Un if len(W[u]) >= 6 and 1.54 * h <= st.median(W[u][len(W[u]) // 3:2 * len(W[u]) // 3]) <= 2.68 * h]
        if hub: out.append(f"Unassigned with 2x depth in their middle third, candidates for shared (D14): {', '.join(hub)}.")
    if os.path.exists(f"{w}/10_groups/eval_groups.txt"): out.append("Reference alarm table: 10_groups/eval_groups.txt. It says where to look, never what to change.")
    return out, al
def r_bridges(w, c, s):
    B = list(rows(f"{w}/12_shared/graph.bridges.tsv")); y = sum(r[2] == "yes" for r in B); T = list(rows(f"{w}/12_shared/graph.tips.tsv"))
    return [f"Bridges: {len(B)}, within a group {y}. Tips: {len(T)} (12_shared/graph.tips.tsv). A telomeric tip needs a graph link and a free end: add it to edits/group_edits.tsv and edits/order_links.tsv, then `reset 08_bridges`."], []
def r_order(w, c, s):
    out = []; al = []
    for f, lab in ((f"{w}/12_shared/eval_order.txt", "Reference alarm"),):
        if os.path.exists(f): out.append(f"{lab}: " + open(f).read().splitlines()[-1])
    u = f"{w}/12_shared/unhonoured.txt"
    if os.path.exists(u): out += ["Reference-free order tests (12_shared/unhonoured.txt):", "```"] + open(u).read().splitlines()[:25] + ["```"]
    return out, al
def r_arms(w, c, s):
    T = list(rows(f"{w}/12_shared/tails.tsv")); C = list(rows(f"{w}/12_shared/tail_checks.tsv")); L = collections.Counter(r[3] for r in T); placed = sum(int(r[2]) - int(r[1]) + 1 for r in rows(f"{w}/12_shared/pre/{c['sample']}.agp") if r[4] == "W"); out = [f"Pre-build placed {Mb(placed)}."]
    out.append(f"Tail candidates: {len(T)}, {Mb(sum(int(r[5]) for r in T))}; labels {dict(L)}. Tests: 12_shared/tail_checks.tsv; passing {sum(r[-1] == 'pass' for r in C)}.")
    out.append("Nothing is placed without edits/accepted_tails.tsv. Keep the tails.tsv columns; the reason goes in the decision.")
    A = list(rows(f"{w}/12_shared/arm_report.tsv"))
    if A: out.append("Arm labels: " + str(dict(collections.Counter(r[3] for r in A if len(r) > 3))))
    return out, []
def r_build(w, c, s):
    p = f"{w}/13_scaffold/{c['sample']}.graph"; L = fai(p + ".fasta.fai"); J = collections.Counter(r[4] for r in rows(p + ".junctions.tsv") if len(r) > 4)
    return [f"Scaffolds {len(L)}, {Mb(sum(L.values()))}. Junction types: {dict(J)}.", "Build log: 13_scaffold/build.log."], []
def fill_reports(w, c):
    return sorted(glob.glob(f"{w}/14_readfill/h*/fill.h*.report.tsv") + glob.glob(f"{w}/14_readfill/*.graph.report.tsv"))
def r_fill(w, c, s):
    K = collections.Counter(); R = fill_reports(w, c)
    for f in R:
        for r in rows(f): K[(r[4], r[5])] += 1
    out = ["Outcomes: " + ", ".join(f"{k[0]}/{k[1]} {v}" for k, v in K.most_common())]
    st_ = f"{w}/15_final/stats.txt"
    if os.path.exists(st_): out.append(open(st_).read().strip())
    al = []
    lm = []
    for f in R:
        for r in rows(f):
            if r[5] == "no" and int(r[7]) >= 10 and int(r[6]) < 2: lm.append(r[0])
    if lm: al.append(f"Many low-MAPQ spanning reads at {lm}: consider edits/fill_retry.tsv (D11).")
    return out, al
def r_raw(w, c, s):
    a = f"{w}/17_raw_reads/audit.tsv"; V = list(rows(a)); al = []; ok = [r for r in V if r[-1] == "pass"]; out = [f"Joins audited: {len(V)}; suggested pass {len(ok)}. Table: 17_raw_reads/audit.tsv."]
    for r in V:
        if r[-1] != "pass": al.append(f"{r[0]} ({r[1]}, {r[2]}): {r[-1]}")
    C = collections.defaultdict(lambda: [0, 0, 0, 0])
    for r in rows(f"{w}/17_raw_reads/cmp.sites.tsv", "kind"):
        pass
    cmp = f"{w}/17_raw_reads/cmp.sites.tsv"
    if os.path.exists(cmp):
        H = open(cmp).readline().rstrip("\n").split("\t"); S = [dict(zip(H, l.rstrip("\n").split("\t"))) for l in open(cmp).readlines()[1:]]
        for k in ("gap", "ctl"):
            x = [g for g in S if g["kind"] == k]; R = sum(int(g["raw_spanq20"]) for g in x); O = sum(int(g["corr_spanq20"]) for g in x)
            out.append(f"{k}: q20 spanning reads raw {R}, corrected {O}, kept {100 * O / max(R, 1):.0f}%.")
    return out, al
def r_apply(w, c, s):
    f = f"{w}/15_final_raw/stats.txt"; return [open(f).read().strip()] if os.path.exists(f) else ["Final files unchanged."], []
def r_validate(w, c, s):
    V = f"{w}/16_validate"; out = []; al = []
    if os.path.exists(f"{V}/stats.txt"): out.append(open(f"{V}/stats.txt").read().strip())
    T = list(rows(f"{V}/telomeres.final.tsv", "#")); T = [r for r in T if len(r) >= 5 and r[0] != "seq"]
    out.append(f"Telomeric ends {sum(r[4].count('T') for r in T)} of {2 * len(T)}; telomere to telomere {sum(r[4] == 'TT' for r in T)} of {len(T)}.")
    for x in sorted(glob.glob(f"{V}/compleasm/*/summary.txt")):
        name = x.split("/")[-2]; z = [l.strip() for l in open(x) if re.match(r"^[SDFIMN]:", l)]; out.append(f"compleasm {name}: {' '.join(z)}")
    ct = f"{V}/maps/specific/crosstalk.tsv"
    if os.path.exists(ct):
        i = o = 0
        for r in rows(ct, "c"):
            try: i += float(r[2]); o += float(r[3])
            except (ValueError, IndexError): pass
        out.append(f"Specific-pair homolog fraction: {o / max(i + o, 1):.4f}")
    # homologs wholly inverted: a direction flag, not biology
    for f in sorted(glob.glob(f"{V}/sv/inv/*_vs_*.paf")):
        fw = rv = 0
        for l in open(f):
            q = l.split("\t")
            if int(q[10]) >= 100000:
                if q[4] == "+": fw += int(q[10])
                else: rv += int(q[10])
        if rv > fw: al.append(f"{os.path.basename(f)}: more inverted ({Mb(rv)}) than forward ({Mb(fw)}): check direction before calling anything an inversion; edits/flip.txt fixes direction (D12).")
    m = f"{V}/sv/inv/merged.tsv"
    if os.path.exists(m): out.append(f"Inverted homolog blocks to test: {sum(1 for _ in open(m))} (sv/inv/merged.tsv, sv/bp_reads.tsv, sv/inversion_test.txt).")
    out.append("Hand the user the Juicebox files and the SV table (D13).")
    return out, al
FN = {"01_inputs": r_inputs, "03_depth": r_depth, "04_chains": r_chains, "05_hic_chains": r_hic, "06_evidence": r_evidence, "07_group": r_group, "08_bridges": r_bridges, "09_order": r_order, "10_arms_tails": r_arms, "11_build": r_build, "12_fill": r_fill, "13_rawreads": r_raw, "14_apply_joins": r_apply, "15_validate": r_validate}

def write(w, sid, c, s, active=lambda cond: True):
    st_ = BY_ID[sid]; out, al = FN[sid](w, c, s) if sid in FN else ([], [])
    with open(f"{w}/reports/{sid}.md", "w") as o:
        o.write(f"# {sid}: {st_.title}\n\n")
        for l in out: o.write(l + "\n\n" if not l.startswith(("```", "  ")) and not re.match(r"^\s*\d+x ", l) else l + "\n")
        if al:
            o.write("## Alarms\n\n"); [o.write(f"- {x}\n") for x in al]; o.write("\n")
        if st_.gates:
            o.write("## Decisions\n\n"); [o.write(f"- {g}: docs/decisions/{GUIDES[g]}\n") for g, cond in st_.gates if active(cond)]
