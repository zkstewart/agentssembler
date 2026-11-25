# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""One short report per stage: the numbers to look at, alarms, and what to decide. Reports describe; they never decide."""
import math, collections, glob, gzip, json, os, re, statistics as st
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
    L = fai(f"{w}/00_inputs/unitigs.fasta.fai"); out = [f"Unitigs: {len(L)}, {Mb(sum(L.values()))}, N50 {Mb(n50(L.values()))}", f"Chromosomes {c['chromosomes']}, ploidy {c['ploidy']}."]
    g = f"{w}/01_genes/summary.txt"
    if not os.path.exists(g): out.append("No compleasm on the unitigs (compleasm not set): later stages cannot count the genes each decision gains or loses.")
    else:
        z = [l.strip() for l in open(g) if re.match(r"^[SDFM]:", l)]
        out.append(f"compleasm on the unitigs (01_genes): {' '.join(z)}. This is what the assembly could hold: at ploidy 2 with distinct haplotypes, most genes are D here. "
                   "Stages 05, 07 and 10 lift these hits onto chains.")
    return out, []
def r_depth(w, c, s):
    L = fai(f"{w}/00_inputs/unitigs.fasta.fai")
    C = [(r[0], float(r[2]), float(r[6])) for r in rows(f"{w}/04_copynumber/class_depth.tsv")]
    # one scale for all three histograms, from the data: the top bin (cap, and everything above) sits past 1.5 x the 99th percentile of
    # primary depth, at least 200x, in about 60 bins at most. A fixed 200x cap held 2x of a sample at 84x per haplotype in one bin.
    q = sorted(p for _, p, _c in C)[int(0.99 * (len(C) - 1))] if C else 0
    cap = max(200, 4 * math.ceil(1.5 * q / 4)); bw = 4 * max(1, math.ceil(cap / 240)); cap = bw * math.ceil(cap / bw)
    def hist(title, pairs):
        H = collections.Counter()
        for d, n in pairs: H[min(int(d // bw) * bw, cap)] += n
        t = sum(H.values()) or 1; o = [f"{title} ({Mb(t)}), percentage per {bw}x bin (the last bin holds {cap}x and above):", "```"]
        o += [f"{k:4d}x {100 * H[k] / t:5.1f}% {'#' * int(100 * H[k] / t)}" for k in sorted(H) if H[k] / t >= 0.002]
        return o + ["```"], max(H, key=H.get) if H else 0
    out = []
    o, _ = hist("Primary depth, all unitigs", [(p, 10000) for _, p, _c in C]); out += o
    long = min(1000000, n50(L.values())) if L else 1000000      # a genome whose unitigs never reach 1 Mb still needs a histogram of its long ones
    o, peak = hist(f"Primary depth, unitigs of {'1 Mb' if long == 1000000 else f'{long / 1e3:.0f} kb (their N50)'} or more", [(p, 10000) for u, p, _c in C if L.get(u, 0) >= long]); out += o
    o, _ = hist("Class depth, adding near-identical copies, all unitigs", [(k, 10000) for _, _p, k in C]); out += o
    # Everything below uses the main peak in long unitigs, and only describes. The peak is a candidate for hap_depth, not a decision.
    h = peak + bw // 2 if peak else 0
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
        P = max(int(c["ploidy"]), 2); hi = (P + 0.68) * h; span = "About 2x" if P == 2 else f"From 2x to {P}x"
        two = sorted(((L[u], u, p) for u, (p, k) in M.items() if 1.54 * h <= p <= hi and k <= hi), reverse=True)
        out += [f"Taking {h}x as a candidate haplotype depth (main bin in long unitigs; check it):",
                f"- Below 0.64x alone but 1x as a class, one haplotype represented more than once: {Mb(twin)}.",
                f"- Below 0.64x even as a class, truly low depth: {Mb(low)}.",
                f"- {span} in both primary and class depth, candidates for collapse: {len(two)} unitigs, {Mb(sum(x[0] for x in two))}, of which graph degree at least 2 at both ends: {Mb(sum(x[0] for x in two if deg[(x[1], '+')] >= 2 and deg[(x[1], '-')] >= 2))}."]
        out += [f"  - {u} {Mb(n)} at {p:.0f}x, degree {deg[(u, '+')]},{deg[(u, '-')]}" for n, u, p in two[:8]]
    out.append(f"Set hap_depth, repeat_depth (about {max(4.4, max(int(c['ploidy']), 2) + 1.4):.1f} x hap_depth: above anything {c['ploidy']} haplotypes can share) and collapsed. See D01.")
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
    pin = s.get("length_scale") if isinstance(s, dict) else None; f = f"{w}/state/length_scale.json"
    if pin not in (None, "auto"): out.append(f"Length scale {float(pin):g}, set in config.yaml. Seeds for grouping are chains of {Mb(1e6 * float(pin))} or more.")
    elif os.path.exists(f):
        d = json.load(open(f)); k = d["scale"]
        out.append(f"Length scale {k:g}: mean chromosome {Mb(d['mean_chromosome'])} ({Mb(d['chain_bp'])} of chains over {c['chromosomes']} x {c['ploidy']}) against 36 Mb for citrus, chain N50 {Mb(d['chain_n50'])}; set by {d['limited_by']}. "
                   f"It multiplies the chromosome-scale windows: seeds for grouping are chains of {1000 * k:.0f} kb or more, chimera tests look {500 * k:.0f} kb each side of a junction, contact maps use {250 * k:.0f} kb bins (SETTINGS, `length_scale`).")
        if d["limited_by"] == "chain N50": al.append(f"The length scale is held down by a chain N50 of {Mb(d['chain_n50'])}: chains are short for chromosomes of {Mb(d['mean_chromosome'])}, so expect many chains per group and a slow ordering.")
        if k > 1.5 or k < 0.67: al.append(f"Length scale {k:g} is outside anything run so far (citrus, 0.97 to 1.03). The scaled windows are untested: read the chimera table and the seed count with that in mind.")
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
    al = []; H = list(rows(f"{w}/10_groups/hub_cores.d04.tsv")); yes = [r for r in H if r[-1] == "yes"]
    if yes: al.append(f"{len(yes)} breakable collapsed hubs (10_groups/hub_cores.d04.tsv): " + ", ".join(f"{r[0]} at {r[10]}/{r[11]}" for r in yes[:12]) + (" ..." if len(yes) > 12 else "")
                      + ". A 2x core of unique reads inside a chain whose arms read 1x. Breaking now, with the chimera breaks, saves a `reset 06_evidence --remap` at D06. "
                      "Positions are on 08_final_chains, the chains D04 breaks; 07_chains2/depth.windows.tsv holds the depth on stage 04's chains, before D03's breaks, so "
                      "check each edge there by the unitig at the break (08_final_chains/chains.agp), then mark the core shared at D14.")
    if len(H) > len(yes): out.append(f"Hub cores not breakable here: {len(H) - len(yes)} (inside one unitig, or the freed piece would carry 1x arms); stage 06 lists them again for D06 and D14.")
    rs = f"{w}/09_hic_map/rescue_signal.txt"
    if os.path.exists(rs):
        L_ = open(rs).read().splitlines(); out += ["Hi-C allele rescue signal (09_hic_map/rescue_signal.txt, D04):"] + [x for x in L_ if not x.startswith("ALARM")]
        on = isinstance(s, dict) and s.get("hic_allele_rescue")
        al += [x[len("ALARM: "):] + f" The setting is now {'on' if on else 'off'}." for x in L_ if x.startswith("ALARM")]
    return out, al
def r_evidence(w, c, s):
    L = fai(f"{w}/12_broken/chains.fasta.fai"); K = collections.Counter(r[6] for r in rows(f"{w}/11_ont_chains/chain_class.tsv")); B = collections.Counter()
    for r in rows(f"{w}/11_ont_chains/chain_class.tsv"): B[r[6]] += int(r[1])
    T = sum(1 for r in rows(f"{w}/16_validate/telomeres.chains.tsv") if "T" in r[4])
    out = [f"Chains: {len(L)}, {Mb(sum(L.values()))}, N50 {Mb(n50(L.values()))}", "Depth class: " + ", ".join(f"{k} {K[k]} ({Mb(B[k])})" for k in K), f"Chains with a telomere call: {T}."]; al = []
    pin = s.get("min_within") if isinstance(s, dict) else None; f = f"{w}/state/hic_density.json"; mw = None
    if pin not in (None, "auto"): mw = float(pin); out.append(f"min_within {mw:g} pairs/Mb, set in config.yaml or by a decision.")
    elif os.path.exists(f):
        d = json.load(open(f)); mw = d["min_within"]
        if "hic_density" in d:
            out.append(f"Hi-C density H = {d['hic_density']:.0f} pairs/Mb within a chain (length-weighted median of the {d['chains_used']} chains of {d['seed_len'] / 1e3:.0f} kb or more). min_within = {d['share']} x H = {mw:g}: a chain below it carries no placement information (SETTINGS, `min_within`).")
            if d["hic_density"] < 2350: al.append(f"Hi-C density {d['hic_density']:.0f} pairs/Mb is under half the lowest run so far (4,700). Few pairs place anything: expect weak grouping and much left unplaced, and read merge scores with that in mind (D05).")
        else: out.append(f"min_within {mw:g} pairs/Mb: {d['source']}.")
    if mw is not None:
        lo = [(int(r[1]), r[0]) for r in rows(f"{w}/10_groups/chain_stats.tsv") if float(r[2]) < mw]
        out.append(f"Chains below min_within: {len(lo)}, {Mb(sum(n for n, _ in lo))}; of 1 Mb or more: {sum(n >= 1e6 for n, _ in lo)}. Mostly satellite arrays: they do not seed grouping, and are ordered after the rest.")
    return out, al
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
    conn = {r[0] for r in rows(f"{w}/12_broken/connector_chains.tsv")}
    Un = [r[0] for r in rows(f"{w}/10_groups/groups.tsv") if r[2] == "unassigned" and r[0] not in conn]
    if conn: out.append(f"Connector chains (restored as possible bridge paths, not in any group): {len(conn)}.")
    ut = [f"{u} ({tel[u]})" for u in Un if "T" in tel.get(u, "--")]
    if ut and any(T[k] < 2 for k in gs): al.append(f"Groups with under 2 telomeric ends, and unassigned chains with a telomere: {', '.join(ut)}. Test them against those groups (D06).")
    # Long overlaps pull a hub chain's median depth down, so its class can read 1x. The middle third is the fairer figure.
    h = float(s.get("hap_depth", 0)) if isinstance(s, dict) else 0
    if h:
        W = collections.defaultdict(list)
        for r in rows(f"{w}/11_ont_chains/depth.windows.tsv"): W[r[0]].append(float(r[3]))
        hub = [f"{u} ({st.median(W[u][len(W[u]) // 3:2 * len(W[u]) // 3]):.0f}x)" for u in Un if len(W[u]) >= 6 and 1.54 * h <= st.median(W[u][len(W[u]) // 3:2 * len(W[u]) // 3]) <= (max(int(c["ploidy"]), 2) + 0.68) * h]
        if hub: out.append(f"Unassigned with {'2x' if int(c['ploidy']) <= 2 else '2x or more'} depth in their middle third, candidates only; test each with review/hub_depth.py (D14): {', '.join(hub)}.")
    dos = list(rows(f"{w}/10_groups/groups.dosage.tsv"))
    if dos:
        K = collections.Counter(r[8] for r in dos); arb = [f"{r[0]} ({r[4]} of {r[4] if r[6] == '.' else int(r[4]) + len(r[6].split(','))}, {r[3]}x)" for r in dos if r[8] == "arbitrary"]
        out.append(f"Shared chains in fewer haplotypes than their chromosome has (10_groups/groups.dosage.tsv): {len(dos)}, {Mb(sum(int(r[1]) for r in dos))}; settled by " + ", ".join(f"{k} {v}" for k, v in sorted(K.items())) + ".")
        if arb: al.append(f"{len(arb)} shared chains whose haplotypes nothing separates: {', '.join(arb[:12])}{' ...' if len(arb) > 12 else ''}. Each is a guess. Test with graph position and interior depth, then name the haplotypes with a shared_in edit or unassign the chain (D14).")
    forks = f"{w}/10_groups/chain_forks.txt"
    if os.path.exists(forks):
        n = sum(1 for l in open(forks) if l.startswith("HUB")); b = sum(1 for l in open(forks) if l.startswith("BUBBLE"))
        if b: out.append(f"{b} BUBBLE rows in 10_groups/chain_forks.txt: 2x-class hubs that read under 2x on the chain itself (hub_total), so allele bubbles, not collapses. No action.")
        if n: al.append(f"{n} HUB rows in 10_groups/chain_forks.txt: 2x unitigs in a chain whose other exit reaches a homolog's chain, reading 2x on the chain itself (or with no depth to say). "
                        "Read them before going on: a break costs a read mapping (D04, D06).")
    G = {r[0]: r[2] for r in rows(f"{w}/10_groups/groups.tsv")}
    HC = collections.defaultdict(list)
    for r in rows(f"{w}/10_groups/hub_cores.tsv"): HC[r[0]].append(r)
    if HC:
        fmt = lambda kv: f"{kv[0]} ({sum(int(r[6]) for r in kv[1])} kb, {G.get(kv[0], '?')})"
        hc = sorted(HC.items(), key=lambda kv: -sum(int(r[6]) for r in kv[1])); yes = [kv for kv in hc if any(r[-1] == "yes" for r in kv[1])]; no = [kv for kv in hc if kv not in yes]
        if yes: al.append(f"{len(yes)} chains classed single-copy hold a breakable core of 2x unique reads (10_groups/hub_cores.tsv): " + ", ".join(map(fmt, yes[:15])) + (" ..." if len(yes) > 15 else "")
                          + ". Each is a collapsed hub whose arms read 1x on MAPQ 0 reads. Unless a chain is already shared, break it at break_left/break_right, mark the core shared, and `reset 06_evidence --remap` in one round (D04, D06, D14).")
        why = lambda kv: "inside one unitig" if any("inside" in r[12] for r in kv[1]) else "freed piece would carry 1x arms"
        if no: al.append(f"{len(no)} chains classed single-copy hold a 2x unique-read core that no D04 break can separate cleanly (10_groups/hub_cores.tsv): "
                         + ", ".join(f"{fmt(kv)[:-1]}, {why(kv)})" for kv in no[:15]) + (" ..." if len(no) > 15 else "")
                         + ". Inside one unitig: break_chain splits only between unitigs. Freed piece with 1x arms: a break would put those arms on both haplotypes (core_2x_frac below 0.75). "
                         "Each chain's core is missing from one haplotype as placed. The choice is D14's: leave it on one haplotype, or mark the chain shared and duplicate its 1x arms; weigh the genes on each.")
    RS = [r for r in rows(f"{w}/10_groups/groups.rescued.tsv") if r[4] in ("collapsed", "redundant")]
    if RS:
        k = collections.defaultdict(list)
        for r in RS: k[r[4]].append(r[0])
        out.append("Rescue left unassigned (10_groups/groups.rescued.tsv): " + "; ".join(f"{x} {len(v)} ({', '.join(v[:10])}{' ...' if len(v) > 10 else ''})" for x, v in sorted(k.items()))
                   + ". collapsed: a breakable 2x unique-read core, both haplotypes, for D06 to break (10_groups/hub_cores.tsv). redundant: a near-copy of sequence every group already holds.")
    for f, lab in ((f"{w}/10_groups/rescue_preview.txt", "Rescue preview (hic_allele_rescue off; what it would place, in 10_groups/rescue_preview/)"), (f"{w}/10_groups/genes_placed.txt", "Genes in the groups")):
        if os.path.exists(f): out.append(f"{lab}: " + " ".join(open(f).read().split()))
    gp = f"{w}/10_groups/genes_placed.txt"
    if os.path.exists(gp) and not (s.get("hic_allele_rescue") if isinstance(s, dict) else False):
        m = re.search(r"short genes made whole: (\d+) of (\d+)", open(gp).read())
        if m and int(m.group(1)) >= max(10, 0.1 * int(m.group(2))):
            al.append(f"The Hi-C allele rescue would make {m.group(1)} of {m.group(2)} short genes whole (10_groups/genes_placed.txt). If D04 did not settle it, decide now: "
                      "record D05 with `--set hic_allele_rescue=true` and `reset 07_group --keepDecision`, which keeps that setting (stages 05-06 already made its pairs). Later, it costs everything from here.")
    if os.path.exists(f"{w}/10_groups/eval_groups.txt"): out.append("Reference alarm table: 10_groups/eval_groups.txt. It says where to look, never what to change.")
    return out, al
def r_bridges(w, c, s):
    B = list(rows(f"{w}/12_shared/graph.bridges.tsv")); y = sum(r[2] == "yes" for r in B); T = list(rows(f"{w}/12_shared/graph.tips.tsv"))
    al = []; cl = f"{w}/12_shared/copies.log"; wn = [l for l in open(cl) if l.startswith("# WARNING")] if os.path.exists(cl) else []
    if wn: al.append(f"{len(wn)} warnings in 12_shared/copies.log about the haplotypes given to shared chains, such as a graph link into a haplotype that was left out. Read them (D14).")
    return [f"Bridges: {len(B)}, within a group {y}. Tips: {len(T)} (12_shared/graph.tips.tsv). A telomeric tip needs a graph link and a free end: add it to edits/group_edits.tsv and edits/order_links.tsv, then `reset 08_bridges`."], al
def r_order(w, c, s):
    out = []; al = []
    for f, lab in ((f"{w}/12_shared/eval_order.txt", "Reference alarm"),):
        if os.path.exists(f): out.append(f"{lab}: " + open(f).read().splitlines()[-1])
    u = f"{w}/12_shared/unhonoured.txt"
    if os.path.exists(u): out += ["Reference-free order tests (12_shared/unhonoured.txt):", "```"] + open(u).read().splitlines()[:25] + ["```"]
    # random restarts of the order search, per group (hic_order.py): improved = the first search was a local optimum
    R = []; N = {}
    for f in sorted(glob.glob(f"{w}/12_shared/by_chrom/c*.log")):
        for l in open(f):
            m = re.match(r"^# (\S+): (\d+) ordered", l)
            if m: N[m.group(1)] = int(m.group(2))
            m = re.match(r"^# (\S+): restarts (\d+) of at most \d+ \(patience \d+\), improved (\d+), reached the best order (\d+)", l)
            if m: R.append((m.group(1), int(m.group(2)), int(m.group(3)), int(m.group(4))))
    if R:
        imp = [f"{g} ({i} of {n})" for g, n, i, r in R if i]; none = [f"{g} ({N.get(g, '?')} chains)" for g, n, i, r in R if not r]
        out.append(f"Order restarts (12_shared/by_chrom/*.log): {len(R)} groups, {sum(n for _, n, _, _ in R)} restarts; improved the first search in {len(imp)}"
                   + (f": {', '.join(imp)}" if imp else "") + f". Groups no restart reached the best order: {len(none)}" + (f" ({', '.join(none[:12])})" if none else "") + ".")
        if imp: al.append(f"Random restarts found better orders than the first search in {', '.join(imp)}: the search has local optima there, so the order may still be one. "
                          "Read those groups' orientation_llr and position_llr (by_chrom/cN.order.tsv) before trusting junctions with small values (D09).")
    return out, al
def r_arms(w, c, s):
    T = list(rows(f"{w}/12_shared/tails.tsv")); C = list(rows(f"{w}/12_shared/tail_checks.tsv")); L = collections.Counter(r[3] for r in T); placed = sum(int(r[2]) - int(r[1]) + 1 for r in rows(f"{w}/12_shared/pre/{c['sample']}.agp") if r[4] == "W"); out = [f"Pre-build placed {Mb(placed)}."]
    out.append(f"Tail candidates: {len(T)}, {Mb(sum(int(r[5]) for r in T))}; labels {dict(L)}. Tests: 12_shared/tail_checks.tsv; passing {sum(r[-1] == 'pass' for r in C)}.")
    out.append("Nothing is placed without edits/accepted_tails.tsv. Keep the tails.tsv columns; the reason goes in the decision.")
    tg = list(rows(f"{w}/12_shared/tail_genes.tsv"))
    if tg: out.append("Genes on tails (12_shared/tail_genes.tsv; short = no full set of copies placed elsewhere): " + ", ".join(f"{r[0]} {r[4]} ({r[5]} short)" for r in tg if int(r[4])) + ". "
                      "A tail that fails the tests keeps its genes as unlocalized sequence (stage 15).")
    UG = [r for r in rows(f"{w}/12_shared/unplaced_genes.tsv") if int(r[4])]
    if UG: out.append(f"Unplaced chains, off every tail, carrying genes short of a full set placed (12_shared/unplaced_genes.tsv): {len(UG)} chains, {sum(int(r[4]) for r in UG)} genes; most: "
                      + ", ".join(f"{r[0]} {r[4]}{' via LOST from ' + r[5] if r[5] != '.' else ''}" for r in UG[:8]) + ". The cost of leaving them out, for the D16 tangle check; most stay unlocalized sequence (stage 15).")
    PH = list(rows(f"{w}/12_shared/phase_reads.pre.tsv", "scaffold"))
    if PH:
        K = collections.Counter(r[-1] for r in PH); out.append("Arbitrary arms in the pre-build, by raw reads (12_shared/phase_reads.pre.tsv, settings.raw_phase_early): " + ", ".join(f"{k} {v}" for k, v in K.most_common()) + ".")
    A = list(rows(f"{w}/12_shared/arm_report.tsv"))
    if A: out.append("Arm labels: " + str(dict(collections.Counter(r[3] for r in A if len(r) > 3))))
    al = direction_alarms(w, "12_shared/homolog_direction.tsv")
    al += swap_alarms(f"{w}/12_shared/rescue_swaps.tsv", pre=True)
    al += [f"{r[0]}:{r[1]}: arm phase conflict by raw reads, reads to specific sequence here {int(r[8]) + int(r[9])}, on the homolog {r[11]} (12_shared/phase_reads.pre.tsv). "
           "Look before D15: the arm's choice may be the other haplotype's." for r in PH if r[-1] == "conflict"]
    FB = sum(1 for r in rows(f"{w}/12_shared/chain_forks.tsv") if r[0] == "BUBBLE")
    if FB: out.append(f"{FB} BUBBLE rows in 12_shared/chain_forks.tsv: forks that would be HUB or LOST, but the chain reads under 2x over the hub (allele bubbles). No action.")
    F = [r for r in rows(f"{w}/12_shared/chain_forks.tsv") if r[0] == "LOST"]
    if F:
        pairs = sorted({(r[1], r[9], int(r[13])) for r in F}, key=lambda x: -x[2])
        al.append(f"{len(pairs)} lost-route forks in 12_shared/chain_forks.tsv (LOST): a placed or tail chain runs through a 2x hub whose other exit reaches an unplaced 1x chain: "
                  + ", ".join(f"{c}->{t} ({n // 1000} kb)" for c, t, n in pairs[:12]) + ". Do the tangle check in D16 (docs/TANGLES.md) before choosing tails: a break at break_at needs `reset 06_evidence --remap`.")
    return out, al
def r_build(w, c, s):
    p = f"{w}/13_scaffold/{c['sample']}.graph"; L = fai(p + ".fasta.fai"); J = collections.Counter(r[4] for r in rows(p + ".junctions.tsv") if len(r) > 4)
    al = direction_alarms(w, "13_scaffold/homolog_direction.tsv")
    return [f"Scaffolds {len(L)}, {Mb(sum(L.values()))}. Junction types: {dict(J)}.", "Build log: 13_scaffold/build.log."], al + swap_alarms(f"{w}/13_scaffold/rescue_swaps.tsv")
def direction_alarms(w, p):
    """Homologs written in opposite directions (review/homolog_direction.py), less those edits/flip.txt already lists."""
    fl = {r[0] for r in rows(f"{w}/edits/flip.txt")}
    inv = [r for r in rows(f"{w}/{p}") if r[4] == "inverted" and not {r[0], r[1]} & fl]
    return [f"{r[0]} is written in reverse against {r[1]} ({Mb(int(r[3]))} inverted, {Mb(int(r[2]))} forward; {p}). A direction flag, not biology: "
            "list one of the two in edits/flip.txt now, before 12_fill, which applies it. Found after stage 12, it costs `reset 12_fill` (D12)." for r in inv]
def swap_alarms(p, pre=False):
    """Hi-C gaps that swapping two chains between homologs would close (review/rescue_swaps.py). pre: the stage 10 pre-build, before any
    scaffold is built, when a swap costs only stages 08-10; ask the user then rather than at D13."""
    S = list(rows(p))
    if not S: return []
    when = ("Ask the user now (D16): before the build a swap is a group edit and `reset 08_bridges` (stages 08-10); after it, the same swap also reruns stages 11-15." if pre
            else "Hand the choice to D13, unless the user already answered at D16 (a swap is a group edit and `reset 08_bridges`).")
    return [f"{len(S)} Hi-C gaps would close if two chains swapped homologs ({'/'.join(p.split('/')[-2:])}): " + ", ".join(f"{r[0]}:{r[1]} {r[4]} ({r[5]}) <-> {r[6]} ({r[8]})" for r in S[:10])
            + (" ..." if len(S) > 10 else "") + ". The graph links both pairs the other way. Neither the graph link nor a Hi-C allele call is proof where haplotypes are near-identical: "
            "read the junction reads. " + when]
def fill_reports(w, c):
    return sorted(glob.glob(f"{w}/14_readfill/h*/fill.h*.report.tsv") + glob.glob(f"{w}/14_readfill/*.graph.report.tsv"))
def homolog_lengths(fa_fai, agp, minBp=1000000, minFrac=0.05):
    """Scaffolds much shorter than their longest homolog. Sequence may lie behind a gap on the short side (D11)."""
    H, gaps = collections.defaultdict(dict), collections.Counter()
    if not os.path.exists(fa_fai): return []
    for n, L in fai(fa_fai).items():
        k, _, h = n.rpartition("_h")
        if k: H[k][n] = L
    if os.path.exists(agp):
        for r in rows(agp):
            if len(r) > 4 and r[4] in ("N", "U"): gaps[r[0]] += 1
    out = []
    for k, v in sorted(H.items()):
        top = max(v.values())
        for n, L in sorted(v.items()):
            if top - L >= minBp and (top - L) / top >= minFrac: out.append(f"{n} is {Mb(top - L)} ({100 * (top - L) / top:.0f}%) shorter than its longest homolog, " + (f"with {gaps[n]} gap{'s' if gaps[n] > 1 else ''}: look behind, or say why not (D11)." if gaps[n] else "with no gaps: check its telomeres before calling the difference real (D11)."))
    return out
def tail_alarms(p):
    """Gap edges where at least 2 reads, and at least half of them, continue into one unplaced sequence. Tails entering a scaffold are mostly homologs at divergent ends, or repeats: counted, not listed."""
    out = []; other = set(); dup = []
    if os.path.exists(p):
        for r in rows(p):
            if r[4] == "." or r[4] == r[0] or int(r[5]) < 2 or 2 * int(r[5]) < int(r[3]): continue
            if r[4].startswith("unplaced") and len(r) > 11 and 2 * int(r[11]) >= int(r[5]): dup.append(f"{r[0]}:{r[1]} {r[2]} -> {r[4][9:]}")
            elif r[4].startswith("unplaced"): out.append(f"{r[0]}:{r[1]} {r[2]} edge: {r[5]} of {r[3]} q20 read tails enter {r[4]} ({p.split('/')[-2]}/gap_tails.tsv). A candidate only: test the graph path, homolog and Hi-C (D11).")
            else: other.add((r[0], r[1], r[2]))
    if dup: out.append(f"{len(dup)} gap edges have read tails entering unplaced sequence that duplicates placed sequence (also_placed, {p.split('/')[-2]}/gap_tails.tsv): {', '.join(dup[:8])}{' ...' if len(dup) > 8 else ''}. Not candidates for the gap.")
    if other: out.append(f"{len(other)} gap edges have read tails mostly entering another scaffold ({p.split('/')[-2]}/gap_tails.tsv). Usually homologs at divergent ends, or repeats.")
    return out
def component_alarms(p, minMb=1.0):
    """Large unplaced tangles with graph contact to a scaffold."""
    out = []
    if os.path.exists(p):
        for l in open(p):
            m = re.match(r"(\d+) chains, ([\d.]+) Mb; graph contacts: (\{.*\})", l)
            if m and float(m[2]) >= minMb: out.append(f"Unplaced tangle of {m[1]} chains, {m[2]} Mb, graph contacts {m[3]} (unplaced_components.txt). Report it; do not place it on this alone.")
    return out
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
    sm = c["sample"]; al += homolog_lengths(f"{w}/15_final/{sm}.final.fasta.fai", f"{w}/15_final/{sm}.final.agp") + tail_alarms(f"{w}/14_readfill/gap_tails.tsv")
    return out, al
def r_raw(w, c, s):
    a = f"{w}/17_raw_reads/audit.tsv"; V = list(rows(a)); al = []; ok = [r for r in V if r[-1] == "pass"]; out = [f"Joins audited: {len(V)}; suggested pass {len(ok)}. Table: 17_raw_reads/audit.tsv."]
    for r in V:
        if r[-1] != "pass": al.append(f"{r[0]} ({r[1]}, {r[2]}): {r[-1]}")
    cmp = f"{w}/17_raw_reads/cmp.sites.tsv"
    if os.path.exists(cmp):
        H = open(cmp).readline().rstrip("\n").split("\t"); S = [dict(zip(H, l.rstrip("\n").split("\t"))) for l in open(cmp).readlines()[1:]]
        for k in ("gap", "ctl"):   # per site, not pooled: one repeat array once held 37 of 40 spanners (T13)
            x = [g for g in S if g["kind"] == k]
            if not x: continue
            R = [int(g["raw_spanq20"]) for g in x]; O = [int(g["corr_spanq20"]) for g in x]
            out.append(f"{k}: {len(x)} sites; with a q20 spanning read: raw {sum(r > 0 for r in R)}, corrected {sum(o > 0 for o in O)}; "
                       f"median q20 spanners per site raw {st.median(R):g}, corrected {st.median(O):g}.")
    al += tail_alarms(f"{w}/17_raw_reads/gap_tails.tsv")
    ph = list(rows(f"{w}/17_raw_reads/phase_reads.tsv", "scaffold"))
    if ph:
        K = collections.Counter(r[-1] for r in ph); out.append(f"Arbitrary arms, by raw reads (phase_reads.tsv): " + ", ".join(f"{k} {v}" for k, v in K.most_common()) + ".")
        al += [f"{r[0]}:{r[1]}: {r[-1]}, reads to specific sequence here {int(r[8]) + int(r[9])}, on the homolog {r[11]} (D15)." for r in ph if r[-1] == "conflict"]
    return out, al
def r_apply(w, c, s):
    f = f"{w}/15_final_raw/stats.txt"; return [open(f).read().strip()] if os.path.exists(f) else ["Final files unchanged."], []
def completeness_alarms(V, c, minComplete=95.0, maxRecoverable=1.0):
    """Low gene completeness on the placed scaffolds, and gene copies that unplaced sequence could add: a placement problem, not a sequence one."""
    al = []; bus = c.get("compleasm") or {}; lo = float(bus.get("min_complete", minComplete))
    for x in sorted(glob.glob(f"{V}/compleasm/h*/summary.txt")):
        name = x.split("/")[-2]
        if "_" in name: continue
        z = {m.group(1): float(m.group(2)) for m in (re.match(r"^([SD]):([\d.]+)%", l) for l in open(x)) if m}
        u = f"{V}/compleasm/{name}_unloc/summary.txt"; zu = {}
        if os.path.exists(u): zu = {m.group(1): float(m.group(2)) for m in (re.match(r"^([SD]):([\d.]+)%", l) for l in open(u)) if m}
        cu = zu.get("S", 0) + zu.get("D", 0) if zu else None
        if z and (cu if cu is not None else z.get("S", 0) + z.get("D", 0)) < lo:
            al.append(f"compleasm {name}: complete {z.get('S', 0) + z.get('D', 0):.2f}%" + (f", {cu:.2f}% with its unlocalized scaffolds" if cu is not None else "") + f", under {lo:g}% (compleasm.min_complete). Read complement.txt: genes complete in unplaced "
                      "sequence point to placement, not to missing sequence (D09, D16; unloc/ for chains with a known chromosome).")
    cg = f"{V}/compleasm/complement.txt"
    if os.path.exists(cg):
        m = re.search(r"short of \d+ per gene, on any scaffold: (\d+) of (\d+) \(([\d.]+)%\), in (\d+) genes; (\d+) of them complete in unplaced", open(cg).read())
        if m and float(m.group(3)) > maxRecoverable and int(m.group(5)) > 0:
            al.append(f"{m.group(1)} gene copies ({m.group(3)}%) short on the scaffolds, in {m.group(4)} genes; {m.group(5)} of those genes are complete in unplaced sequence. "
                      "Find which unplaced chains carry them (compleasm/unplaced) and why they were not placed: Hi-C rescue, tails (D16), tangles, or unloc/.")
    return al
def r_validate(w, c, s):
    V = f"{w}/16_validate"; out = []; al = []
    if os.path.exists(f"{V}/stats.txt"): out.append(open(f"{V}/stats.txt").read().strip())
    T = list(rows(f"{V}/telomeres.final.tsv", "#")); T = [r for r in T if len(r) >= 5 and r[0] != "seq"]
    out.append(f"Telomeric ends {sum(r[4].count('T') for r in T)} of {2 * len(T)}; telomere to telomere {sum(r[4] == 'TT' for r in T)} of {len(T)}.")
    for x in sorted(glob.glob(f"{V}/compleasm/*/summary.txt")):
        name = x.split("/")[-2]; z = [l.strip() for l in open(x) if re.match(r"^[SDFIMN]:", l)]; out.append(f"compleasm {name}: {' '.join(z)}")
    cg = f"{V}/compleasm/complement.txt"
    if os.path.exists(cg): out.append(open(cg).read().strip().replace("\n", "\n\n"))
    al += completeness_alarms(V, c)
    sc = f"{V}/juicebox/out_JBAT.scale"
    if os.path.exists(sc) and open(sc).read().strip() not in ("", "1"):
        al.append(f"The Juicebox map is scaled by {open(sc).read().strip()} (assembly over 2.1 Gb; 16_validate/juicebox/out_JBAT.scale). Say so in the hand-over, and check that scaffold edges fall where the .assembly puts them: a scaled map has not been tried here.")
    ut = f"{V}/unloc/unloc.txt"
    if os.path.exists(ut): out.append("Unlocalized scaffolds, part of the assembly handed over (leaving them out is the user's call, D13): " + open(ut).read().strip().replace("\n", "\n\n"))
    ct = f"{V}/maps/specific/crosstalk.tsv"
    if os.path.exists(ct):
        i = o = 0
        for r in rows(ct):
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
    K = collections.defaultdict(list)
    for r in rows(f"{V}/depth.2x_chains.tsv"): K[(r[4], r[0])].append(float(r[3]))
    HC = {r[0] for r in rows(f"{w}/10_groups/hub_cores.tsv")}
    k2 = sorted(((c, sc, v) for (c, sc), v in K.items() if len(v) >= 2), key=lambda t: -len(t[2]))
    if k2: al.append(f"{len(k2)} chains placed as haplotype-specific read about 2x over 100 kb or more (16_validate/depth.2x_chains.tsv): "
                     + ", ".join(f"{c} on {sc} ({50 * len(v)} kb{', hub core' if c in HC else ''})" for c, sc, v in k2[:15]) + (" ..." if len(k2) > 15 else "")
                     + ". Each holds both haplotypes collapsed: a hub to break and share (10_groups/hub_cores.tsv, D06, D14), or a shared chain placed once (D14). "
                     f"Single windows ({sum(1 for v in K.values() if len(v) == 1)}) are counted, not listed.")
    m = f"{V}/sv/inv/merged.tsv"
    if os.path.exists(m): out.append(f"Inverted homolog blocks to test: {sum(1 for _ in open(m))} (sv/inv/merged.tsv, sv/bp_reads.tsv, sv/inversion_test.txt).")
    al += swap_alarms(f"{w}/13_scaffold/rescue_swaps.tsv")
    F = "15_final_raw" if os.path.exists(f"{w}/15_final_raw/{c['sample']}.final.agp") else "15_final"
    al += homolog_lengths(f"{w}/{F}/{c['sample']}.final.fasta.fai", f"{w}/{F}/{c['sample']}.final.agp") + component_alarms(f"{V}/unplaced_components.txt")
    # Hand tables set aside by a reset. Coordinates and haplotype numbers can change in a re-run, so these are lists to re-check, never rows to copy back.
    old = [os.path.relpath(f, w) for pat in ("16_validate/sv/sv_calls.tsv", "16_validate/curation_targets.tsv") for f in sorted(glob.glob(f"{w}/archive/*/{pat}")) if any(True for _ in rows(f))]
    if old: al.append("Hand tables from before a reset: " + ", ".join(old) + ". Re-test each row against this build and write the result in the new tables. Do not copy rows across: scaffold names and coordinates may have changed (D12, D13).")
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
