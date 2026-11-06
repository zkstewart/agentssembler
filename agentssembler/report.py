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
    W = collections.Counter(); tot = 0
    for r in rows(f"{w}/04_copynumber/class_depth.tsv"):
        d = float(r[2]); b = int(d // 4) * 4; n = int(r[1]) ; W[min(b, 200)] += 10000; tot += 10000
    # own primary depth histogram, weighted by window
    O = collections.Counter()
    for r in rows(f"{w}/04_copynumber/depth.windows.tsv"): O[min(int(float(r[3]) // 4) * 4, 200)] += int(r[2]) - int(r[1])
    t = sum(O.values()) or 1
    top = sorted(O, key=lambda k: -O[k])[:4]
    out = ["Primary depth, percentage of bases per 4x bin:", "```"]
    for k in sorted(O):
        if O[k] / t >= 0.002: out.append(f"{k:4d}x {100 * O[k] / t:5.1f}% {'#' * int(100 * O[k] / t)}")
    out += ["```", f"Largest bins start at: {top}.", "Set hap_depth, repeat_depth (about 4.4 x hap_depth) and collapsed: is there a peak at about twice the haplotype depth in long unitigs?"]
    return out, []
def r_chains(w, c, s):
    out = []; al = []
    for d in ("06_chains", "07_chains2", "08_final_chains"):
        L = fai(f"{w}/{d}/chains.fasta.fai")
        if L: out.append(f"{d}: {len(L)} chains, {Mb(sum(L.values()))}, N50 {Mb(n50(L.values()))}")
    V = list(rows(f"{w}/07_chains2/joins.verify.tsv")); low = [r for r in V if int(r[5]) < 10]; br = [r for r in V if int(r[5]) < 3]
    out.append(f"Joins {len(V)}; under 3 reads {len(br)} (broken by rule); under 10 reads {len(low)} (review in joins.verify.tsv).")
    if os.path.exists(f"{w}/07_chains2/restored_arms.tsv"): out.append(f"Restored arms: {sum(1 for _ in rows(f'{w}/07_chains2/restored_arms.tsv'))}.")
    return out, al
def r_hic(w, c, s):
    p = f"{w}/10_groups/chain_breaks.tsv"; V = list(rows(p))
    return [f"Candidate positions in {p}: {len(V)}. A break needs two signals: different Hi-C partners each side and weak reads across the join."], []
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
