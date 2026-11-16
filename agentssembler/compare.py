# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Side-by-side metrics of two finished work directories, for choosing between builds at D13 (for example a TTT build,
docs/TANGLES.md): validation report lines, unplaced sequence, reused arms, depth by class, and length, gaps and telomeres per
scaffold. Scaffold names must match for the per-scaffold rows to line up."""
import os, re, collections

def final_dir(w, sample):
    return f"{w}/15_final_raw" if os.path.exists(f"{w}/15_final_raw/{sample}.final.agp") else f"{w}/15_final"

def grab(path, pats):
    out = {}
    if not os.path.exists(path): return out
    txt = open(path).read()
    for k, pat in pats.items():
        m = re.search(pat, txt, re.M); out[k] = m.group(1).strip() if m else "-"
    return out

def metrics(w, sample):
    M = collections.OrderedDict(); F = final_dir(w, sample)
    M.update(grab(f"{w}/reports/15_validate.md", {
        "assembly": r"final\.agp (scaffolds=.*)$", "telomeres": r"^(Telomeric ends.*)$",
        "compleasm h1": r"^compleasm h1: (.*)$", "compleasm h2": r"^compleasm h2: (.*)$",
        "compleasm unplaced": r"^compleasm unplaced: (.*)$", "complete per haplotype": r"^Complete per haplotype number: (.*)$",
        "gene copies short": r"^Gene copies short of 2 per gene, on any scaffold: (.*)$",
        "specific-pair homolog fraction": r"^Specific-pair homolog fraction: (.*)$",
        "inverted blocks": r"^Inverted homolog blocks to test: (\d+)"}))
    uf = f"{F}/{sample}.final.unplaced.fasta.fai"
    if os.path.exists(uf):
        L = [int(l.split("\t")[1]) for l in open(uf)]; M["unplaced"] = f"{len(L)} seqs, {sum(L)/1e6:.2f} Mb"
    au = f"{w}/12_shared/arm_usage.tsv"
    if os.path.exists(au):
        ov = [l.split("\t") for l in open(au) if l.rstrip().endswith("over")]
        M["reused arms (over)"] = f"{len(ov)} chains, {sum(int(x[1]) for x in ov)/1e6:.2f} Mb: " + ",".join(x[0][5:] for x in ov)
    dl = f"{w}/16_validate/depth_classes.log"
    if os.path.exists(dl):
        for l in open(dl):
            m = re.match(r"\s*(specific|shared|mixed)\b.*?(\d+) windows, median ([\d.]+).*?>1.5x median ([\d.]+), CV ([\d.]+)", l)
            if m: M[f"depth {m.group(1)}"] = f"{m.group(2)} windows, median {m.group(3)}, >1.5x {m.group(4)}, CV {m.group(5)}"
            m = re.match(r"\s+(\S+)\t([\d.]+) Mb\t([\d.]+) Mb\t([\d.]+)", l)
            if m: M[f"in both haplotypes {m.group(1)}"] = f"{m.group(3)} Mb ({m.group(4)})"
    agp = f"{F}/{sample}.final.agp"
    if os.path.exists(agp):
        ln = collections.Counter(); gp = collections.Counter()
        for l in open(agp):
            s = l.split("\t")
            if l[0] == "#": continue
            ln[s[0]] = max(ln[s[0]], int(s[2]))
            if s[4] in "NU": gp[s[0]] += 1
        for s in sorted(ln): M[f"len {s}"] = f"{ln[s]/1e6:.2f} Mb, {gp[s]} gaps"
    tf = f"{w}/16_validate/telomeres.final.tsv"
    if os.path.exists(tf):
        for l in open(tf):
            x = l.split()
            if l[0] != "#" and len(x) >= 5 and f"len {x[0]}" in M: M[f"len {x[0]}"] += f", tel {x[-1]}"
    return M

def table(wa, wb, sa, sb):
    """Markdown table, one row per metric, marking rows that differ."""
    A, B = metrics(wa, sa), metrics(wb, sb)
    keys = list(A) + [k for k in B if k not in A]
    out = [f"| Metric | {os.path.basename(wa.rstrip('/'))} | {os.path.basename(wb.rstrip('/'))} |", "|---|---|---|"]
    for k in keys:
        va, vb = A.get(k, "-"), B.get(k, "-")
        out.append(f"| {k} | {va} | {vb}{' **(diff)**' if va != vb else ''} |")
    return "\n".join(out)
