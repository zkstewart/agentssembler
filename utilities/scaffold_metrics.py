#!/usr/bin/env python3
'''
scaffold_metrics.py -- contiguity metrics to rank scaffolding methods fairly.

Subcommands

  whole   Goal 1: contiguity of each tool's complete output, with all haplotypes and
          unplaced sequence pooled together.
          Reads annotarium 'fasta stats' TSVs named <genome>.[<part>.]<tool>.stats.tsv
          from one directory. Every file of a genome+tool is pooled, e.g.
          Busk.hap1.agentssembler + Busk.hap2.agentssembler + Busk.unplaced.agentssembler.
          <genome>.p_utg.stats.tsv is the shared input. Its total length is the
          denominator G for the NG and auNG values, unless --genome_size overrides it.
          Input contig sets (p_utg, p_ctg) are listed as baselines; --input_map says which
          one each tool was given (default p_utg), for pct_of_input.
          n_gaps = N bases / --gap_size (assumes one run of exactly gap_size Ns per gap).

          auN  = sum(L^2) / sum(L)  (area under the Nx curve; no single threshold)
          auNG = sum(L^2) / G       (the same, with one fixed denominator per genome)
          Within a genome G is the same for every tool, so the auNG ranking does not depend
          on G; only the absolute values do. NG50 rankings can depend on G.

  placed  Goal 2: contiguity of the chromosome-scale sequence, from the labels.v2.tsv of
          'hic_crosstalk.py prepare' (xt/<tool>/prepare/labels.v2.tsv). There every sequence
          has a chromosome (by alignment, the same way for every tool) and a class: backbone
          (chromosome-scale, >= its --scaffold_len), fragment, or unplaced. No haplotype
          labels are used. Per chromosome, the --ploidy longest backbones are its copies
          ("anchored"); further backbones and fragments are pieces left unjoined.
            anchored_bp      total length of the copies
            anchored_dom_bp  the same, times each copy's share of aligned bases on its
                             chromosome (discounts sequence joined from other chromosomes)
            each as % of the tool's output, of its input (--input_map) and of G (p_utg).
            extra_backbones, fragments, pieces_per_copy (1.0 = one piece per copy),
            possible_chimeras (copies with dominance < --min_dominance), and with --stats,
            anchored_gaps and contigs.
          With --compleasm (ideally one run on each tool's whole output), BUSCO genes are
          counted within the anchored sequences only: busco_any (>= 1 copy),
          busco_all_copies (>= ploidy copies), busco_excess (> ploidy: duplication),
          busco_copy_ratio (mean of min(copies, ploidy) / ploidy), busco_frag_only, busco_missing.

Python 3 standard library only.
'''

import argparse, collections, os, re, sys

STATS_RE = re.compile(r"^(?P<genome>[^.]+)\.(?:(?P<part>.+)\.)?(?P<tool>[^.]+)\.stats\.tsv$")
COMPLETE = {"Single", "Duplicated"}
PARTIAL = {"Fragmented", "Interspersed"}

# ------------------------------------------------------------------ utilities ----
def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)

def read_stats(path):
    '''annotarium fasta stats TSV -> list of (id, length, N bases)'''
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("\t") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            n = int(float(cols[3])) if len(cols) > 3 and cols[3] else 0
            rows.append((cols[0], int(float(cols[1])), n))
    return rows

def read_stats_lengths(path):
    '''annotarium fasta stats TSV -> (list of lengths, total N bases)'''
    rows = read_stats(path)
    return [r[1] for r in rows], sum(r[2] for r in rows)

GAP_WARNED = set()
def count_gaps(nBases, gapSize, where):
    '''Gap count inferred as N bases / gapSize. Exact only if every gap is one run of exactly
    gapSize Ns; warns (once per source) when the N count is not a multiple of gapSize.'''
    if nBases % gapSize != 0 and where not in GAP_WARNED:
        GAP_WARNED.add(where)
        print(f"# WARNING: {where} has {nBases} N bases, not a multiple of --gap_size {gapSize}; "
              f"its gap count is not exact", file=sys.stderr)
    return nBases // gapSize

def gap_assumption_note(gapSize):
    print(f"# NOTE: gap counts = N bases / {gapSize}. This assumes every gap is a single run of exactly "
          f"{gapSize} Ns (no longer or merged gaps, no Ns inside contigs). Verify on the FASTA, e.g.\n"
          f"#   grep -v '>' X.fa | tr -d '\\n' | grep -oE '[Nn]+' | awk '{{print length($0)}}' | sort -n | uniq -c",
          file=sys.stderr)

def nx(lengths, threshold):
    '''(Nx length, Lx count): smallest length at which cumulative sum >= threshold'''
    cum = 0
    for i, length in enumerate(sorted(lengths, reverse=True), 1):
        cum += length
        if cum >= threshold:
            return length, i
    return 0, None  # threshold never reached (only possible for NGx)

def contiguity(lengths, G):
    total = sum(lengths)
    sq = sum(length * length for length in lengths)
    n50, l50 = nx(lengths, total * 0.5)
    n90, _ = nx(lengths, total * 0.9)
    ng50, lg50 = nx(lengths, G * 0.5)
    ng90, _ = nx(lengths, G * 0.9)
    return {"n_seqs": len(lengths), "total_bp": total, "pct_of_G": 100 * total / G,
            "largest": max(lengths), "N50": n50, "L50": l50, "N90": n90,
            "NG50": ng50, "LG50": lg50 if lg50 is not None else ".", "NG90": ng90,
            "auN": sq / total, "auNG": sq / G}

def fmt(v):
    if isinstance(v, float):
        return f"{v:.2f}" if abs(v) < 1000 else f"{v:.0f}"
    return str(v)

def write_table(rows, cols, path):
    with (open(path, "w") if path else sys.stdout) as out:
        out.write("\t".join(cols) + "\n")
        for r in rows:
            out.write("\t".join(fmt(r.get(c, ".")) for c in cols) + "\n")

def print_pretty(rows, cols, title):
    print(f"\n## {title}", file=sys.stderr)
    cells = [cols] + [[fmt(r.get(c, ".")) for c in cols] for r in rows]
    widths = [max(len(row[i]) for row in cells) for i in range(len(cols))]
    for row in cells:
        print("  ".join(v.rjust(w) for v, w in zip(row, widths)), file=sys.stderr)

def discover_stats(directory):
    '''-> {genome: {tool: [paths]}}'''
    found = collections.defaultdict(lambda: collections.defaultdict(list))
    for name in sorted(os.listdir(directory)):
        m = STATS_RE.match(name)
        if m:
            found[m.group("genome")][m.group("tool")].append(os.path.join(directory, name))
    return found

def parse_sizes(specs):
    sizes = {}
    for spec in specs or []:
        genome, value = spec.split("=", 1)
        sizes[genome] = int(float(value))
    return sizes

# ---------------------------------------------------------------------- whole ----
def parse_input_map(specs):
    '''TOOL=INPUT, e.g. yahs=p_ctg -> {tool: input}'''
    out = {}
    for spec in specs or []:
        tool, inp = spec.split("=", 1)
        out[tool] = inp
    return out

def pooled_lengths(paths):
    lengths, nBases = [], 0
    for path in paths:
        L, N = read_stats_lengths(path)
        lengths += L
        nBases += N
    return lengths, nBases

def cmd_whole(a):
    found = discover_stats(a.dir)
    overrides = parse_sizes(a.genome_size)
    inputMap = parse_input_map(a.input_map)
    inputs = set(a.inputs)
    gap_assumption_note(a.gap_size)
    rows = []
    for genome in sorted(found):
        tools = found[genome]
        if "p_utg" not in tools and genome not in overrides:
            die(f"no {genome}.p_utg.stats.tsv in {a.dir} and no --genome_size {genome}=...")
        G = overrides.get(genome) or sum(pooled_lengths(tools["p_utg"])[0])
        for tool in sorted(tools, key=lambda t: (t not in inputs, t)):
            lengths, nBases = pooled_lengths(tools[tool])
            for path in tools[tool]:
                for sid, _, n in read_stats(path):
                    count_gaps(n, a.gap_size, f"{os.path.basename(path)}:{sid}")
            inp = None if tool in inputs else inputMap.get(tool, "p_utg")
            if inp is not None and inp not in tools:
                die(f"{genome} {tool}: input set '{inp}' has no {genome}.[*.]{inp}.stats.tsv in {a.dir}")
            inBp = sum(pooled_lengths(tools[inp])[0]) if inp else None
            r = {"genome": genome, "tool": tool, "G": G, "files": len(tools[tool]),
                 "n_bases": nBases, "n_gaps": nBases // a.gap_size,
                 "input": inp or "(input)", "input_bp": inBp if inBp else ".",
                 "pct_of_input": 100 * sum(lengths) / inBp if inBp else "."}
            r.update(contiguity(lengths, G))
            rows.append(r)
            print(f"# {genome} {tool}: pooled {', '.join(os.path.basename(p) for p in tools[tool])}"
                  + (f"; input = {inp}" if inp else ""), file=sys.stderr)
    cols = ["genome", "tool", "files", "n_seqs", "total_bp", "input", "input_bp", "pct_of_input",
            "G", "pct_of_G", "n_bases", "n_gaps", "largest",
            "N50", "L50", "N90", "NG50", "LG50", "NG90", "auN", "auNG"]
    write_table(rows, cols, a.out)
    pretty = [{**r, **{k: r[k] / 1e6 for k in ("total_bp", "largest", "N50", "NG50", "NG90", "auN", "auNG")}}
              for r in rows]
    print_pretty(pretty, ["genome", "tool", "n_seqs", "n_gaps", "total_bp", "input", "pct_of_input",
                          "pct_of_G", "N50", "NG50", "LG50", "auN", "auNG"],
                 "Whole-output contiguity (lengths in Mb; G = p_utg total unless overridden)")

# --------------------------------------------------------------------- placed ----
def read_labels(path):
    '''labels.v2.tsv from hic_crosstalk.py prepare'''
    rows = []
    with open(path) as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        if "class" not in header:
            die(f"{path} is not a labels.v2.tsv (no 'class' column); rerun the updated hic_crosstalk.py prepare")
        for line in fh:
            if line.strip():
                r = dict(zip(header, line.rstrip("\n").split("\t")))
                r["length"] = int(r["length"])
                r["dominance"] = float(r["dominance"])
                rows.append(r)
    return rows

def parse_assembly_specs(specs, what):
    '''NAME=PATH[,PATH...] -> {name: [paths]}'''
    out = collections.OrderedDict()
    for spec in specs or []:
        if "=" not in spec:
            die(f"{what} must be NAME=PATH, not '{spec}'")
        name, paths = spec.split("=", 1)
        out.setdefault(name, []).extend(paths.split(","))
    return out

def build_id_lookup(labels):
    '''Sequence ID -> label row, by the tool's own source_id first, then by labels name/target_id.'''
    byName = {}
    for r in labels:
        byName[r["name"]] = r
        byName[r["target_id"]] = r
    bySource = collections.defaultdict(list)
    for r in labels:
        bySource[r["source_id"]].append(r)
    return byName, bySource

def resolve(seqID, byName, bySource, srcFilter):
    hits = bySource.get(seqID, [])
    if srcFilter:
        hits = [r for r in hits if srcFilter in r["source_file"]]
    if len(hits) == 1:
        return hits[0]
    if len(hits) > 1:
        raise LookupError(f"'{seqID}' matches {len(hits)} source sequences; "
                          f"give the file as FILE@<source file substring>")
    return byName.get(seqID)

def score_compleasm(tables, labels, anchoredIDs, ploidy):
    '''Copies of each BUSCO gene within the anchored sequences -> percentages'''
    byName, bySource = build_id_lookup(labels)
    complete = collections.Counter()
    partial = set()
    genes, unmatched = set(), collections.Counter()
    for spec in tables:
        path, _, srcFilter = spec.partition("@")
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            iG, iS, iQ = header.index("Gene"), header.index("Status"), header.index("Sequence")
            for line in fh:
                c = line.rstrip("\n").split("\t")
                genes.add(c[iG])
                if c[iS] not in COMPLETE and c[iS] not in PARTIAL:
                    continue
                try:
                    r = resolve(c[iQ], byName, bySource, srcFilter)
                except LookupError as e:
                    die(f"{path}: {e}")
                if r is None:
                    unmatched[c[iQ]] += 1
                    continue
                if r["target_id"] not in anchoredIDs:
                    continue
                if c[iS] in COMPLETE:
                    complete[c[iG]] += 1
                else:
                    partial.add(c[iG])
    if unmatched:
        ex = ", ".join(list(unmatched)[:3])
        print(f"# WARNING: {len(unmatched)} compleasm sequence ID(s) not found in labels (e.g. {ex}); "
              f"their hits are ignored", file=sys.stderr)
    N = len(genes)
    if not N:
        return {}
    k = {g: complete.get(g, 0) for g in genes}
    return {"busco_any": 100 * sum(1 for g in genes if k[g] >= 1) / N,
            "busco_all_copies": 100 * sum(1 for g in genes if k[g] >= ploidy) / N,
            "busco_excess": 100 * sum(1 for g in genes if k[g] > ploidy) / N,
            "busco_copy_ratio": sum(min(k[g], ploidy) for g in genes) / (N * ploidy),
            "busco_frag_only": 100 * sum(1 for g in genes if k[g] == 0 and g in partial) / N,
            "busco_missing": 100 * sum(1 for g in genes if k[g] == 0 and g not in partial) / N}

def natural(c):
    return [int(x) if x.isdigit() else x for x in re.split(r"(\d+)", c)]

def cmd_placed(a):
    labelSpecs = parse_assembly_specs(a.labels, "--labels")
    compSpecs = parse_assembly_specs(a.compleasm, "--compleasm")
    statSpecs = parse_assembly_specs(a.stats, "--stats")
    inputMap = parse_input_map(a.input_map)
    found = discover_stats(a.stats_dir) if a.stats_dir else {}
    if statSpecs:
        gap_assumption_note(a.gap_size)
    P = a.ploidy

    chromRows, asmRows = [], []
    for name, paths in labelSpecs.items():
        if len(paths) != 1:
            die(f"--labels {name} needs exactly one labels.v2.tsv")
        genome, _, tool = name.partition(":")
        if not tool:
            die(f"--labels name must be GENOME:TOOL, not '{name}'")
        labels = read_labels(paths[0])
        allBp = sum(r["length"] for r in labels)

        hasGaps = name in statSpecs
        if hasGaps:
            byName, bySource = build_id_lookup(labels)
            for spec in statSpecs[name]:
                path, _, srcFilter = spec.partition("@")
                for sid, _, n in read_stats(path):
                    try:
                        r = resolve(sid, byName, bySource, srcFilter)
                    except LookupError as e:
                        die(f"{path}: {e}")
                    if r is not None:
                        r["gaps"] = count_gaps(n, a.gap_size, f"{os.path.basename(path)}:{sid}")
            missing = [r["source_id"] for r in labels if r["class"] != "unplaced" and "gaps" not in r]
            if missing:
                die(f"{name}: {len(missing)} placed sequence(s) not found in its --stats file(s), "
                    f"e.g. '{missing[0]}'")

        byChrom = collections.defaultdict(list)
        for r in labels:
            if r["class"] != "unplaced":
                byChrom[r["chrom"]].append(r)
        anchoredIDs, tot = set(), collections.Counter()
        for chrom in sorted(byChrom, key=natural):
            rows = byChrom[chrom]
            bbs = sorted((r for r in rows if r["class"] == "backbone"), key=lambda r: -r["length"])
            frs = [r for r in rows if r["class"] == "fragment"]
            anchored = bbs[:P]  # the ploidy longest chromosome-scale sequences = the chromosome copies
            anchoredIDs.update(r["target_id"] for r in anchored)
            chromBp = sum(r["length"] for r in rows)
            aBp = sum(r["length"] for r in anchored)
            aDom = sum(r["length"] * r["dominance"] for r in anchored)
            chim = sum(1 for r in anchored if r["dominance"] < a.min_dominance)
            row = {"genome": genome, "tool": tool, "chrom": chrom, "copies": len(anchored),
                   "extra_backbones": len(bbs) - len(anchored), "fragments": len(frs),
                   "pieces": len(rows), "chrom_bp": chromBp, "anchored_bp": aBp,
                   "anchored_frac": aBp / chromBp if chromBp else 0.0,
                   "copy_lengths": ",".join(str(r["length"]) for r in anchored) or ".",
                   "copy_sources": ",".join(r["source_id"] for r in anchored) or ".",
                   "possible_chimeras": chim}
            if hasGaps:
                row["anchored_gaps"] = sum(r["gaps"] for r in anchored)
                row["contigs"] = sum(r["gaps"] + 1 for r in rows)
                tot["anchored_gaps"] += row["anchored_gaps"]
                tot["contigs"] += row["contigs"]
            chromRows.append(row)
            for k in ("copies", "extra_backbones", "fragments", "pieces", "chrom_bp", "possible_chimeras"):
                tot[k] += row[k]
            tot["anchored_bp"] += aBp
            tot["anchored_dom_bp"] += aDom
            tot["chroms"] += 1

        r = {"genome": genome, "tool": tool, "chroms": tot["chroms"], "expected_copies": tot["chroms"] * P,
             "copies": tot["copies"], "extra_backbones": tot["extra_backbones"], "fragments": tot["fragments"],
             "pieces_per_copy": tot["pieces"] / (tot["chroms"] * P) if tot["chroms"] else 0.0,
             "possible_chimeras": tot["possible_chimeras"], "all_bp": allBp, "placed_bp": tot["chrom_bp"],
             "anchored_bp": tot["anchored_bp"], "anchored_dom_bp": round(tot["anchored_dom_bp"]),
             "placed_pct_output": 100 * tot["chrom_bp"] / allBp,
             "anchored_pct_output": 100 * tot["anchored_bp"] / allBp,
             "anchored_dom_pct_output": 100 * tot["anchored_dom_bp"] / allBp}
        if hasGaps:
            r.update(anchored_gaps=tot["anchored_gaps"], contigs=tot["contigs"])
        tools = found.get(genome, {})
        if "p_utg" in tools:
            G = sum(pooled_lengths(tools["p_utg"])[0])
            r.update(G=G, anchored_pct_G=100 * tot["anchored_bp"] / G,
                     anchored_dom_pct_G=100 * tot["anchored_dom_bp"] / G)
        inp = inputMap.get(tool, "p_utg")
        if inp in tools:
            inBp = sum(pooled_lengths(tools[inp])[0])
            r.update(input=inp, input_bp=inBp, anchored_pct_input=100 * tot["anchored_bp"] / inBp)
        if name in compSpecs:
            r.update(score_compleasm(compSpecs[name], labels, anchoredIDs, P))
        asmRows.append(r)

    os.makedirs(a.outdir, exist_ok=True)
    cCols = ["genome", "tool", "chrom", "copies", "extra_backbones", "fragments", "pieces", "contigs",
             "anchored_gaps", "chrom_bp", "anchored_bp", "anchored_frac", "possible_chimeras",
             "copy_lengths", "copy_sources"]
    aCols = ["genome", "tool", "chroms", "expected_copies", "copies", "extra_backbones", "fragments",
             "pieces_per_copy", "contigs", "anchored_gaps", "possible_chimeras", "all_bp", "input", "input_bp",
             "G", "placed_bp", "anchored_bp", "anchored_dom_bp", "placed_pct_output", "anchored_pct_output",
             "anchored_dom_pct_output", "anchored_pct_input", "anchored_pct_G", "anchored_dom_pct_G",
             "busco_any", "busco_all_copies", "busco_copy_ratio", "busco_excess", "busco_frag_only",
             "busco_missing"]
    write_table(chromRows, cCols, os.path.join(a.outdir, "placed.per_chromosome.tsv"))
    write_table(asmRows, aCols, os.path.join(a.outdir, "placed.summary.tsv"))
    print_pretty(asmRows, ["genome", "tool", "copies", "expected_copies", "extra_backbones", "fragments",
                           "anchored_gaps", "possible_chimeras", "anchored_pct_output", "anchored_pct_input",
                           "anchored_pct_G", "busco_all_copies", "busco_excess"],
                 f"Anchored-sequence summary (tables in {a.outdir})")

# ----------------------------------------------------------------------- main ----
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    w = sub.add_parser("whole", help="goal 1: pooled whole-output contiguity")
    w.add_argument("-d", "--dir", required=True, help="directory of *.stats.tsv files")
    w.add_argument("-o", "--out", help="output TSV (default stdout)")
    w.add_argument("--genome_size", action="append", metavar="GENOME=BP",
                   help="override G for a genome (default: total length of <genome>.p_utg)")
    w.add_argument("--inputs", nargs="+", default=["p_utg", "p_ctg"],
                   help="'tool' names that are input contig sets, listed as baselines [%(default)s]")
    w.add_argument("--input_map", action="append", metavar="TOOL=INPUT",
                   help="input contig set a tool was given, for pct_of_input (default p_utg), e.g. yahs=p_ctg")
    w.add_argument("--gap_size", type=int, default=100, help="length of one scaffolding gap [%(default)s]")
    w.set_defaults(func=cmd_whole)

    p = sub.add_parser("placed", help="goal 2: chromosome-scale contiguity, from hic_crosstalk labels.v2.tsv")
    p.add_argument("--labels", action="append", required=True, metavar="GENOME:TOOL=labels.v2.tsv",
                   help="repeat for each assembly of one genome, e.g. Busk:yahs=xt/yahs/prepare/labels.v2.tsv")
    p.add_argument("--compleasm", action="append", metavar="GENOME:TOOL=TABLE[@SRC][,TABLE[@SRC]...]",
                   help="""compleasm full_table.tsv file(s) for that assembly. Sequence IDs are matched
                   to the tool's own IDs (labels source_id), else to labels name/target_id; when a
                   source_id occurs in several input files (YaHS hap1 and hap2 both have
                   scaffold_1), add @SRC, a substring of that table's source file path""")
    p.add_argument("--ploidy", type=int, required=True)
    p.add_argument("--min_dominance", type=float, default=0.9,
                   help="a copy counts as a possible chimera below this dominance [%(default)s]")
    p.add_argument("--stats_dir", help="directory of *.stats.tsv: <genome>.p_utg gives G, and the "
                   "--input_map set gives each tool's input, for anchored %% of G / of input")
    p.add_argument("--input_map", action="append", metavar="TOOL=INPUT",
                   help="input contig set a tool was given (default p_utg), e.g. yahs=p_ctg")
    p.add_argument("--stats", action="append", metavar="GENOME:TOOL=TSV[@SRC][,TSV[@SRC]...]",
                   help="""annotarium fasta stats TSV(s) of that tool's output, to count gaps per
                   chromosome copy; IDs are matched like --compleasm tables""")
    p.add_argument("--gap_size", type=int, default=100, help="length of one scaffolding gap [%(default)s]")
    p.add_argument("-o", "--outdir", required=True)
    p.set_defaults(func=cmd_placed)

    a = ap.parse_args()
    a.func(a)

if __name__ == "__main__":
    main()
