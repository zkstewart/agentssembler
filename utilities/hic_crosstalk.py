#!/usr/bin/env python3
'''
hic_crosstalk.py -- compare Hi-C phasing quality (homolog crosstalk, switch errors) and
inter-chromosome crosstalk between phased assemblies, using one method for all of them,
and draw per-chromosome contact maps that also show fragmentation.

Design
  * The reference (default: the agentssembler sequences named <chr>_h1) is used ONLY to give
    every sequence a common chromosome name and an axis position. It is never used to decide
    haplotypes, so no assembler's haplotype labels are imposed on another.
  * Every assembly is labelled the same way, by alignment. Names in the inputs (hapN=,
    backbone= <chr>_h<N>) are recorded as claims and compared, but are not used.
  * Crosstalk measures phasing, not fragmentation. It is scored on each assembly's own
    chromosome-scale sequences ("backbones"); fragments are drawn in the maps and never enter
    the contact rates. Only hic_score counts their bp, as off-backbone, so that leaving hard
    sequence unscaffolded cannot raise the score.

Subcommands (run in this order, or use 'all'; each step is skipped when already done)

  prepare  1. The reference haplotype (--ref_regex) is extracted; its sequences name the
              chromosomes.
           2. annotarium 'assembly phase' aligns the assembly to it (minimap2). Only the
              per-sequence alignment summary is used; its haplotype partitions are ignored.
           3. Labelling, by alignment for every assembly:
                backbone   >= --scaffold_len; chromosome = the reference chromosome with most
                           aligned bases ('chimeric?' when that is < --min_dominance of them:
                           a misjoin, or a real translocation relative to the reference; the
                           label only names and orders the backbone, it is never penalised)
                fragment   shorter, >= --min_part_len, >= --min_dominance to one chromosome
                           and >= --min_aligned_frac of its length aligned
                unplaced   everything else (still in the mapping target)
              Names: <chr>_s<N> for backbones (N by length) and <chr>_f<NNN> for fragments
              (NNN by reference position).
           4. Mapping target = all sequences + decoys (named decoy_*). When an earlier run's
              prepare/mapping_target.fa and labels.tsv exist they are REUSED unchanged (same
              sequences, so the Hi-C mapping stays valid), and only the labels change.
           5. Homology, from the assembly's own sequences (no reference, no chromosome labels):
              each backbone is the minimap2 target for every other chromosome-assigned
              sequence, so a homolog moved to another chromosome by a translocation is found.
                homologous regions  alignment blocks >= --mask_min_aln chained along both
                           sequences (gaps <= --homolog_gap), kept when the chain aligns
                           >= min(--min_homolog_region, --partner_frac x the other sequence's
                           length). Region-level, not sequence-level: a backbone carrying a
                           translocated block is homologous to the other chromosome only in
                           that block (homology.v4.tsv). Pieces of one copy do not overlap, so
                           contacts between them are never penalised.
                testable   backbone sequence covered by >= ploidy-1 other sequences (so not
                           collapsed or hemizygous) AND by >= 1 other backbone (so a phasing
                           error there could be seen). Only testable loci are scored.
                excess     backbone sequence covered by >= ploidy other backbones: more copies
                           than the ploidy allows (duplicated or unpurged), written to excess.v4.bed.
  map      bwa-mem2 mem -5SP -T0 | samtools fixmate -m | sort | markdup -r |
           view -F 0x904 | sort -n, then pairs.tsv.gz: ctgA posA mapqA ctgB posB mapqB
           (read 1 first; both mates MAPQ >= 1 and NM < --max_nm).
  score    Phase-testable = testable minus excess (one mask). Pairs with both MAPQ >= --mapq,
           no decoy and no unplaced end are drawn in the maps; of those, backbone pairs are:
             homologous     two backbones, each end in its backbone's region homologous to
                            the other, whatever the chromosome labels (so a translocation in
                            one haplotype is scored as phasing, not structure)
             interchrom     different chromosome labels and not homologous, NOT masked (a
                            misjoined segment has no homolog on its scaffold's chromosome,
                            so the allelic mask would hide it)           -> structure
             unresolved     same label and not homologous, e.g. two pieces of one copy
                                                                            [not penalised]
             cis            same backbone, both ends phase-testable
             homolog        homologous, both ends phase-testable         -> phasing error
             dropped_masked same backbone or homologous, an end not phase-testable
           Pairs with a fragment end are counted (fragment_pairs) but never scored.
           homolog_fraction = homolog / (cis + homolog)
           interchrom_fraction = interchrom / backbone_pairs
           Bin scans (--bin), one rule for both: per bin, own = same-backbone contacts >= --far
           away, other = contacts with the single most-contacted other sequence. Baseline =
           the assembly's median other/own over bins with own + other >= --min_hom. Flagged
           when other >= --min_hom and other/own > --fold * baseline (a self-calibrated
           threshold: library, ploidy and backbone length shift the baseline, not the flags).
             switch   other = one partner backbone; bins >= --min_testable phase-testable;
                      own and other from phase-testable pairs only   -> switch_runs.tsv
             misjoin  other = one other chromosome (all its backbones); every backbone bin;
                      no mask                                        -> misjoin_runs.tsv
           Adjacent flagged bins form runs (events). A switch flags both haplotypes, and
           up to about --far either side of its breakpoint.
           hic_score: every bp of the assembly (or --genome_size) in one category, BUSCO-like:
             C consistent   backbone, phase-testable, not flagged          <- hic_score = C%
             X conflict     backbone, in a flagged misjoin or switch bin
             D duplicated   backbone, excess copy (more copies than --ploidy), not flagged
             U untestable   other backbone bp: collapsed, hemizygous, or no partner backbone
             O off-backbone fragments, unplaced (and absent, when --genome_size is larger)
           Maps (--contact_map): one per chromosome with its backbones (longest first) then
           its fragments (by reference position), from all MAPQ-filtered pairs without
           masking, so fragmentation and misplacement are visible; plus one map of all
           backbones (--genome_map) for inter-chromosome structure.
  compare  One table of the summary values of several assembly directories.

Requirements: python3 (stdlib), annotarium, minimap2, bwa-mem2 (or bwa), samtools.

Example (one assembly per call/job)
  hic_crosstalk.py all -o xt/agentssembler --annotarium ~/annotarium/annotarium.py \\
      --reference Busk.final.fasta --ploidy 2 \\
      --input Busk.final.fasta --input Busk.final.unplaced.fasta \\
      --decoys decoys.fasta -1 HiC_R1.fq.gz -2 HiC_R2.fq.gz -t 32
  hic_crosstalk.py all -o xt/yahs ... --input hap1_scaffolds_final.fa --input hap2_scaffolds_final.fa
  hic_crosstalk.py compare xt/agentssembler xt/yahs --mapq 20 --genome_size 1300000000
'''

import argparse, bisect, collections, gzip, os, re, shlex, shutil, subprocess, sys, time

COMP = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv", "TGCANtgcanYRMKSWVHDByrmkswvhdb")

# ------------------------------------------------------------------ utilities ----
def log(msg):
    print(f"[{time.strftime('%F %T')}] {msg}", file=sys.stderr, flush=True)

def die(msg):
    log(f"ERROR: {msg}")
    sys.exit(1)

def opener(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def read_fasta(path):
    name, chunks = None, []
    with opener(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)

def write_fasta(fh, name, seq, width=80):
    fh.write(f">{name}\n")
    for i in range(0, len(seq), width):
        fh.write(seq[i:i + width] + "\n")

def revcomp(seq):
    return seq.translate(COMP)[::-1]

def run(cmd, log_file=None):
    '''Run a command list, or a bash pipeline string with pipefail.'''
    shown = cmd if isinstance(cmd, str) else " ".join(shlex.quote(c) for c in cmd)
    log(f"$ {shown}")
    err = open(log_file, "a") if log_file else None
    try:
        if isinstance(cmd, str):
            subprocess.run(["bash", "-c", "set -euo pipefail; " + cmd], check=True, stderr=err)
        else:
            subprocess.run(cmd, check=True, stderr=err)
    except subprocess.CalledProcessError as e:
        die(f"command failed (exit {e.returncode}){'; see ' + log_file if log_file else ''}: {shown}")
    finally:
        if err:
            err.close()

def need(exe):
    path = shutil.which(exe)
    if path is None:
        die(f"'{exe}' not found in PATH")
    return path

def merge(intervals):
    out = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [tuple(x) for x in out]

def intersect(a, b):
    out, i, j = [], 0, 0
    while i < len(a) and j < len(b):
        s, e = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        if s < e:
            out.append((s, e))
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return out

def done(path):
    return os.path.exists(path)

def touch(path):
    open(path, "w").close()


def fasta_names(path):
    with opener(path) as fh:
        return [line[1:].split()[0] for line in fh if line.startswith(">")]

def depth_at_least(intervalLists, k, length):
    '''Regions covered by >= k of the (already merged) interval lists.'''
    if k <= 0:
        return [(0, length)]
    events = []
    for iv in intervalLists:
        for s, e in iv:
            events.append((s, 1))
            events.append((e, -1))
    events.sort(key=lambda x: (x[0], x[1]))  # at one position, ends before starts
    out, d, start = [], 0, None
    for pos, delta in events:
        prev, d = d, d + delta
        if prev < k <= d:
            start = pos
        elif prev >= k > d and pos > start:
            out.append((start, pos))
    return merge(out)

def total(iv):
    return sum(e - s for s, e in iv)

LABEL_COLS = ["name", "target_id", "length", "class", "chrom", "ref_start", "strand", "dominance",
              "aligned_frac", "partners", "testable_bp", "claimed_chrom", "claimed_hap", "source_file",
              "source_id", "note"]

def read_labels(outdir):
    path = os.path.join(outdir, "prepare", "labels.v2.tsv")
    if not os.path.exists(path):
        die(f"{path} not found; run 'prepare' first")
    rows = []
    with open(path) as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            r["length"] = int(r["length"])
            rows.append(r)
    return rows

# -------------------------------------------------------------------- prepare ----
def parse_inputs(specs):
    '''--input [hapN=|backbone=]FILE -> list of (tag, kind, hap, path). Labels are kept only so
    that earlier runs' sequence IDs (and so their mapping targets) are reproduced; the claims they
    carry are recorded, not used.'''
    inputs = []
    for i, spec in enumerate(specs, 1):
        label, path = spec.split("=", 1) if "=" in spec else ("", spec)
        if not os.path.isfile(path):
            die(f"input not found: {path}")
        m = re.fullmatch(r"hap(\d+)", label)
        if m:
            inputs.append((f"{label}in{i}", "hap", int(m.group(1)), path))
        elif label == "backbone":
            inputs.append((f"bb{i}", "backbone", None, path))
        elif label == "":
            inputs.append((f"in{i}", "auto", None, path))
        else:
            die(f"--input label must be hapN or backbone, not '{label}'")
    return inputs

def parse_phasing(path):
    rows = {}
    with open(path) as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            regions = []
            if r["reference_regions"] != ".":
                for reg, bp, st in zip(r["reference_regions"].split(","), r["reference_bp"].split(","),
                                       r["reference_strand"].split(",")):
                    ref, coords = reg.rsplit(":", 1)
                    regions.append((ref, int(bp), st, int(coords.split("-")[0])))
            r["regions"] = regions
            rows[r["contig"]] = r
    return rows

def read_old_labels(path):
    '''labels.v2.tsv of an earlier run of this version, or labels.tsv of the first version:
    (source_file, source_id) -> (target_id, strand in the mapping target)'''
    out, byID = {}, collections.defaultdict(list)
    v2 = path.endswith("labels.v2.tsv")
    with open(path) as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            tid, st = (r["target_id"], r["strand"]) if v2 else (r["new_id"], r["strand"])
            val = (tid, st if st in ("+", "-") else "+")
            out[(r["source_file"], r["source_id"])] = val
            byID[r["source_id"]].append(val)
    return out, byID

def prepare(a):
    P = os.path.join(a.outdir, "prepare")
    os.makedirs(P, exist_ok=True)
    if done(os.path.join(P, "done.v4.ok")):
        log("prepare: already done")
        return
    if not a.reference or not a.input or not a.annotarium:
        die("prepare needs --reference, --input and --annotarium")
    minimap2 = a.minimap2 or need("minimap2")

    # 1. reference haplotype: chromosome names and order only
    refFa = os.path.join(P, "reference_hap.fa")
    if os.path.exists(refFa):
        refOrder = fasta_names(refFa)
        log(f"reusing {refFa}")
    else:
        refRe = re.compile(a.ref_regex)
        refOrder = []
        with open(refFa + ".part", "w") as out:
            for name, seq in read_fasta(a.reference):
                m = refRe.fullmatch(name)
                if m:
                    refOrder.append(m.group(1))
                    write_fasta(out, m.group(1), seq)
        os.replace(refFa + ".part", refFa)
    if not refOrder:
        die(f"no reference sequence matches --ref_regex '{a.ref_regex}'")
    log(f"reference: {len(refOrder)} chromosomes ({', '.join(refOrder[:5])}{', ...' if len(refOrder) > 5 else ''})")

    # 2. inputs with unique, traceable IDs (identical to earlier versions, so old outputs stay valid)
    inputs = parse_inputs(a.input)
    bbRe = re.compile(a.backbone_regex)
    queryFa = os.path.join(P, "query.fa")
    writeQuery = not os.path.exists(queryFa)
    meta = collections.OrderedDict()
    out = open(queryFa + ".part", "w") if writeQuery else None
    for tag, kind, hap, path in inputs:
        for name, seq in read_fasta(path):
            qid = f"{tag}__{name}"
            if qid in meta:
                die(f"duplicate sequence ID {name} in {path}")
            claimedChrom, claimedHap = None, hap
            if kind == "backbone":
                m = bbRe.fullmatch(name)
                if m:
                    claimedChrom, claimedHap = m.group(1), int(m.group(2))
            meta[qid] = {"file": path, "orig": name, "length": len(seq),
                         "claimed_chrom": claimedChrom, "claimed_hap": claimedHap}
            if out:
                write_fasta(out, qid, seq)
    if out:
        out.close()
        os.replace(queryFa + ".part", queryFa)
    log(f"inputs: {len(meta)} sequences from {len(inputs)} file(s)")

    # 3. annotarium phase, for its alignment summary only
    phaseDir = os.path.join(P, "phase")
    phasing = os.path.join(phaseDir, "phasing.tsv")
    if not os.path.exists(phasing):
        cmd = [sys.executable if a.annotarium.endswith(".py") else a.annotarium]
        if a.annotarium.endswith(".py"):
            cmd.append(a.annotarium)
        cmd += ["assembly", "phase", "-i", queryFa, "-r", refFa, "-o", phaseDir,
                "--ploidy", str(a.ploidy), "--preset", a.preset, "--threads", str(a.threads),
                "--minAlignLen", str(a.phase_min_aln), "--minimap2", minimap2]
        run(cmd, os.path.join(P, "annotarium_phase.log"))
    else:
        log(f"reusing {phasing}")
    ph = parse_phasing(phasing)
    missing = [q for q in meta if q not in ph]
    if missing:
        die(f"{len(missing)} input sequence(s) are not in {phasing} (e.g. {missing[0]}); "
            f"were the same --input files given in the same order?")

    # 4. labelling, by alignment, identically for every assembly
    lab = collections.OrderedDict()
    for qid, m in meta.items():
        bpByRef, strandOf, startOf = collections.Counter(), {}, {}
        for ref, bp, st, start in ph[qid]["regions"]:
            bpByRef[ref] += bp
            strandOf.setdefault(ref, st)
            startOf[ref] = min(startOf.get(ref, start), start)
        tot = sum(bpByRef.values())
        best, bestBp = (bpByRef.most_common(1)[0] if bpByRef else (None, 0))
        dom = bestBp / tot if tot else 0.0
        alnFrac = tot / m["length"] if m["length"] else 0.0
        note, cls, chrom = [], "unplaced", None
        if best is None:
            note.append("unaligned")
        elif m["length"] >= a.scaffold_len:
            cls, chrom = "backbone", best
            if dom < a.min_dominance:
                note.append("chimeric?")
        elif dom < a.min_dominance:
            note.append("ambiguous_chromosome")
        elif alnFrac < a.min_aligned_frac:
            note.append("low_aligned_fraction")
        elif m["length"] < a.min_part_len:
            note.append("too_short")
        else:
            cls, chrom = "fragment", best
        if m["claimed_chrom"] and chrom and m["claimed_chrom"] != chrom:
            note.append(f"claimed_{m['claimed_chrom']}")
        lab[qid] = {"class": cls, "chrom": chrom, "dom": dom, "alnFrac": alnFrac,
                    "strand": strandOf.get(chrom, "+") if chrom else "+",
                    "ref_start": startOf.get(chrom, 0) if chrom else 0, "note": note}

    refIndex = {c: i for i, c in enumerate(refOrder)}
    for chrom in refOrder:
        bbs = sorted((q for q, x in lab.items() if x["chrom"] == chrom and x["class"] == "backbone"),
                     key=lambda q: -meta[q]["length"])
        for n, q in enumerate(bbs, 1):
            lab[q]["name"] = f"{chrom}_s{n}"
        frs = sorted((q for q, x in lab.items() if x["chrom"] == chrom and x["class"] == "fragment"),
                     key=lambda q: (lab[q]["ref_start"], -meta[q]["length"]))
        for n, q in enumerate(frs, 1):
            lab[q]["name"] = f"{chrom}_f{n:03d}"
    n = 0
    for q, x in lab.items():
        if "name" not in x:
            n += 1
            x["name"] = f"unplaced_{n:05d}"

    # 5. mapping target: reuse an earlier one unchanged, or write a new one
    target = os.path.join(P, "mapping_target.fa")
    # labels.v2.tsv: written by an earlier run of this version; labels.tsv: by the first version
    oldLabels = next((os.path.join(P, f) for f in ("labels.v2.tsv", "labels.tsv")
                      if os.path.exists(os.path.join(P, f))), None)
    if os.path.exists(target):
        if oldLabels is None:
            die(f"{target} exists but no labels file it was written with (labels.v2.tsv or labels.tsv)")
        old, oldByID = read_old_labels(oldLabels)
        for q, x in lab.items():
            key = (meta[q]["file"], meta[q]["orig"])
            if key in old:
                x["target_id"], x["target_strand"] = old[key]
            elif len(oldByID.get(meta[q]["orig"], [])) == 1:
                x["target_id"], x["target_strand"] = oldByID[meta[q]["orig"]][0]
            else:
                die(f"cannot find {meta[q]['orig']} ({meta[q]['file']}) in {oldLabels}")
        log(f"reusing {target} and its sequence IDs, so the existing Hi-C mapping stays valid")
    else:
        for x in lab.values():
            x["target_id"] = x["name"]
            x["target_strand"] = x["strand"] if x["class"] != "unplaced" else "+"
        byQ = {q: x for q, x in lab.items()}
        seen = set()
        with open(target + ".part", "w") as out:
            for qid, seq in read_fasta(queryFa):
                x = byQ[qid]
                write_fasta(out, x["target_id"], revcomp(seq) if x["target_strand"] == "-" else seq)
                seen.add(x["target_id"])
            nDecoy = 0
            for path in a.decoys or []:
                for name, seq in read_fasta(path):
                    did = f"decoy_{name}" if not name.startswith("decoy_") else name
                    if did in seen:
                        die(f"duplicate decoy ID {did}")
                    seen.add(did)
                    write_fasta(out, did, seq)
                    nDecoy += 1
        os.replace(target + ".part", target)
        log(f"mapping target written with {nDecoy} decoy sequence(s)")

    # 6. homology among the assembly's own sequences, in target orientation. Every backbone is
    #    the minimap2 target for ALL other chromosome-assigned sequences, whatever their label,
    #    so homologs moved to another chromosome by a translocation are still found.
    H = os.path.join(P, "homology")
    os.makedirs(H, exist_ok=True)
    byTarget = {x["target_id"]: (q, x) for q, x in lab.items()}
    assigned = {x["target_id"] for x in lab.values() if x["class"] != "unplaced"}
    assignedFa = os.path.join(H, "assigned.fa")
    with open(assignedFa + ".part", "w") as allFh:
        for tid, seq in read_fasta(target):
            if tid not in assigned:
                continue
            q, x = byTarget[tid]
            if len(seq) != meta[q]["length"]:
                die(f"{tid} in {target} is {len(seq)} bp but {meta[q]['orig']} is {meta[q]['length']} bp")
            write_fasta(allFh, tid, seq)
            if x["class"] == "backbone":
                with open(os.path.join(H, f"bb_{tid}.fa"), "w") as fh:
                    write_fasta(fh, tid, seq)
    os.replace(assignedFa + ".part", assignedFa)

    lengthOf = {x["target_id"]: meta[q]["length"] for q, x in lab.items()}
    classOf = {x["target_id"]: x["class"] for x in lab.values()}
    # Homologous regions, per (backbone b, other sequence o): alignment blocks >= --mask_min_aln
    # chained along b and o (gaps <= --homolog_gap on both), kept when the chain aligns
    # >= min(--min_homolog_region, --partner_frac * len(o)) bp. Chaining in both coordinates
    # keeps scattered repeat hits from adding up to a false homolog.
    blocks = {}  # b -> o -> merged aligned blocks on b (coverage: collapse, testable, excess)
    spans = {}   # b -> backbone o -> merged chain spans on b (which contacts are homolog contacts)
    G = a.homolog_gap
    for q, x in lab.items():
        if x["class"] != "backbone":
            continue
        tid = x["target_id"]
        paf = os.path.join(H, f"bb_{tid}.v4.paf")
        if not os.path.exists(paf + ".ok"):
            run([minimap2, "-x", a.preset, "-t", str(a.threads), "--secondary=no", "-o", paf,
                 os.path.join(H, f"bb_{tid}.fa"), assignedFa], os.path.join(H, "minimap2.log"))
            touch(paf + ".ok")
        hits = collections.defaultdict(list)
        with open(paf) as fh:
            for line in fh:
                s = line.split("\t")
                if s[0] != tid and s[5] == tid and int(s[10]) >= a.mask_min_aln:
                    hits[s[0]].append((int(s[7]), int(s[8]), int(s[2]), int(s[3])))
        blocks[tid], spans[tid] = {}, {}
        for o, hs in hits.items():
            need_bp = min(a.min_homolog_region, a.partner_frac * lengthOf[o])
            keptBlocks, keptSpans, chain = [], [], None
            for ts, te, qs, qe in sorted(hs) + [(None,) * 4]:
                if chain is not None and (ts is None or ts - chain["te"] > G or
                                          max(qs, chain["qs"]) - min(qe, chain["qe"]) > G):
                    iv = merge(chain["iv"])
                    if total(iv) >= need_bp:
                        keptBlocks += iv
                        keptSpans.append((iv[0][0], iv[-1][1]))
                    chain = None
                if ts is None:
                    break
                if chain is None:
                    chain = {"te": te, "qs": qs, "qe": qe, "iv": [(ts, te)]}
                else:
                    chain["te"] = max(chain["te"], te)
                    chain["qs"], chain["qe"] = qs, qe  # next block must be near the last one on o
                    chain["iv"].append((ts, te))
            if keptBlocks:
                blocks[tid][o] = merge(keptBlocks)
                if classOf[o] == "backbone":
                    spans[tid][o] = merge(keptSpans)

    # 7. partners (any homologous region, with any chromosome label) and testable regions
    partners = collections.defaultdict(set)
    partnerRows = []
    for b, others in blocks.items():
        for o, iv in others.items():
            partners[b].add(o)
            partners[o].add(b)
            partnerRows.append((b, o, classOf[o], total(iv)))
    testable, excess, testBp, excessBp = {}, {}, 0, 0
    for b, others in blocks.items():
        L = lengthOf[b]
        bbIv = [iv for o, iv in others.items() if classOf[o] == "backbone"]
        depth = depth_at_least(list(others.values()), a.ploidy - 1, L)
        testable[b] = intersect(depth, depth_at_least(bbIv, 1, L))
        # extra copies: covered by more partner backbones than the ploidy allows (a duplicated or
        # unpurged copy, not a phase switch)
        excess[b] = [iv for iv in depth_at_least(bbIv, a.ploidy, L)  # drop slivers where blocks touch
                     if iv[1] - iv[0] >= a.mask_min_aln]
        testBp += total(testable[b])
        excessBp += total(excess[b])
    for fname, regions in (("testable.v4.bed", testable), ("excess.v4.bed", excess)):
        with open(os.path.join(P, fname), "w") as fh:
            for b in sorted(regions):
                for s, e in regions[b]:
                    fh.write(f"{b}\t{s}\t{e}\n")
    with open(os.path.join(P, "homology.v4.tsv"), "w") as fh:
        fh.write("#backbone\tstart\tend\thomologous_backbone\n")
        for b in sorted(spans):
            for o in sorted(spans[b]):
                for s, e in spans[b][o]:
                    fh.write(f"{b}\t{s}\t{e}\t{o}\n")
    with open(os.path.join(P, "partners.tsv"), "w") as fh:
        fh.write("#backbone\tpartner\tpartner_class\taligned_bp_on_backbone\n")
        for row in partnerRows:
            fh.write("\t".join(str(v) for v in row) + "\n")

    # 8. labels and summary
    nameOrder = lambda q: ((0, refIndex[lab[q]["chrom"]], lab[q]["name"]) if lab[q]["chrom"]
                           else (1, 0, lab[q]["name"]))
    with open(os.path.join(P, "labels.v2.tsv"), "w") as fh:
        fh.write("#" + "\t".join(LABEL_COLS) + "\n")
        for q in sorted(lab, key=nameOrder):
            x, m, tid = lab[q], meta[q], lab[q]["target_id"]
            vals = [x["name"], tid, m["length"], x["class"], x["chrom"] or ".", x["ref_start"],
                    x["target_strand"], f"{x['dom']:.3f}", f"{x['alnFrac']:.3f}",
                    ",".join(sorted(byTarget[p][1]["name"] for p in partners.get(tid, ()))) or ".",
                    total(testable.get(tid, [])), m["claimed_chrom"] or ".", m["claimed_hap"] or ".",
                    m["file"], m["orig"], ",".join(x["note"]) or "."]
            fh.write("\t".join(str(v) for v in vals) + "\n")

    count = collections.Counter(x["class"] for x in lab.values())
    bpOf = collections.Counter()
    for q, x in lab.items():
        bpOf[x["class"]] += meta[q]["length"]
    perChrom = collections.Counter(x["chrom"] for x in lab.values() if x["class"] == "backbone")
    chim = [x["name"] for x in lab.values() if "chimeric?" in x["note"]]
    claimed = [x["name"] for x in lab.values() if any(n.startswith("claimed_") for n in x["note"])]
    noPartner = [x["name"] for x in lab.values() if x["class"] == "backbone"
                 and not any(classOf[p] == "backbone" for p in partners.get(x["target_id"], ()))]
    summary = [("sequences", len(lab)), ("backbones", count["backbone"]), ("backbone_bp", bpOf["backbone"]),
               ("fragments", count["fragment"]), ("fragment_bp", bpOf["fragment"]),
               ("unplaced", count["unplaced"]), ("unplaced_bp", bpOf["unplaced"]),
               ("all_bp", sum(bpOf.values())), ("testable_bp", testBp), ("excess_copy_bp", excessBp),
               ("chroms_without_backbone", ",".join(c for c in refOrder if not perChrom[c]) or "."),
               ("backbones_per_chrom", ",".join(f"{c}:{perChrom[c]}" for c in refOrder)),
               ("backbones_without_partner", ",".join(noPartner) or "."),
               ("possible_chimeras", ",".join(chim) or "."),
               ("chrom_differs_from_claim", ",".join(claimed) or ".")]
    with open(os.path.join(P, "summary.v2.tsv"), "w") as fh:
        for k, v in summary:
            fh.write(f"{k}\t{v}\n")
    log(f"labelled: {count['backbone']} backbones ({bpOf['backbone'] / 1e6:.1f} Mb), "
        f"{count['fragment']} fragments ({bpOf['fragment'] / 1e6:.1f} Mb), {count['unplaced']} unplaced; "
        f"testable {testBp / 1e6:.1f} Mb ({100 * testBp / max(bpOf['backbone'], 1):.1f}% of backbones)")
    for c in refOrder:
        log(f"  {c}: {perChrom[c]} backbone(s)")
    if chim:
        log(f"  WARNING: possible chimeras: {', '.join(chim)}")
    if noPartner:
        log(f"  NOTE: backbones with no partner backbone (no homolog contacts can be scored): {', '.join(noPartner)}")
    if claimed:
        log(f"  NOTE: chromosome differs from the name's claim: {', '.join(claimed)}")
    if excessBp:
        log(f"  NOTE: {excessBp / 1e6:.2f} Mb of backbone has more copies than --ploidy allows (excess.v4.bed)")
    touch(os.path.join(P, "done.v4.ok"))

# ------------------------------------------------------------------------ map ----
PAIRS_AWK = r'''
function ok(   i) {
    if ($5 < 1) return 0
    if (maxnm > 0) for (i = 12; i <= NF; i++)
        if (substr($i, 1, 5) == "NM:i:") return (substr($i, 6) + 0 < maxnm)
    return 1
}
{
    good = ok()
    if ($1 == pname) {
        if (pgood && good) {
            if (int(pflag / 64) % 2) print pctg, ppos, pmq, $3, $4, $5
            else                     print $3, $4, $5, pctg, ppos, pmq
            n++
        }
        pname = ""
    } else { pname = $1; pflag = $2; pctg = $3; ppos = $4; pmq = $5; pgood = good }
}
END { print n + 0 " pairs written" > "/dev/stderr" }'''

def map_reads(a):
    P, Mp = os.path.join(a.outdir, "prepare"), os.path.join(a.outdir, "map")
    os.makedirs(Mp, exist_ok=True)
    if done(os.path.join(Mp, "done.ok")):
        log("map: already done")
        return
    if not done(os.path.join(P, "done.v4.ok")):
        die("run 'prepare' first")
    if not a.r1 or not a.r2:
        die("map needs -1 and -2")
    r1s, r2s = a.r1.split(","), a.r2.split(",")
    if len(r1s) != len(r2s):
        die("-1 and -2 need the same number of files")
    need("samtools")
    target = os.path.join(P, "mapping_target.fa")
    lg = os.path.join(Mp, "map.log")
    q = shlex.quote
    if a.aligner == "bwa-mem2":
        need("bwa-mem2")
        if not os.path.exists(target + ".bwt.2bit.64"):
            run(["bwa-mem2", "index", target], lg)
    else:
        need("bwa")
        if not os.path.exists(target + ".bwt"):
            run(["bwa", "index", target], lg)
    bam = os.path.join(Mp, "hic.namesorted.bam")
    tmp = os.path.join(Mp, "tmp")
    os.makedirs(tmp, exist_ok=True)
    if not done(bam + ".ok"):
        in1 = q(r1s[0]) if len(r1s) == 1 else "<(cat " + " ".join(q(x) for x in r1s) + ")"
        in2 = q(r2s[0]) if len(r2s) == 1 else "<(cat " + " ".join(q(x) for x in r2s) + ")"
        t, mem = a.threads, a.sort_mem
        run(f"{a.aligner} mem -t {t} -5SP -T0 {q(target)} {in1} {in2} 2>> {q(lg)}"
            f" | samtools fixmate -m --output-fmt bam,level=0 - -"
            f" | samtools sort -u -@ {t} -m {mem} -T {q(tmp + '/coord')} -"
            f" | samtools markdup -r -@ {t} -f {q(os.path.join(Mp, 'markdup.txt'))} --output-fmt bam,level=0 - -"
            f" | samtools view -u -F 0x904 -"
            f" | samtools sort -n -@ {t} -m {mem} -T {q(tmp + '/name')} -o {q(bam)} -", lg)
        touch(bam + ".ok")
    pairs = os.path.join(Mp, "pairs.tsv.gz")
    awk = os.path.join(Mp, "pairs.awk")
    with open(awk, "w") as fh:
        fh.write(PAIRS_AWK)
    run(f"samtools view -@ {a.threads} -F 0xF0C {q(bam)}"
        f" | awk -F'\\t' -v OFS='\\t' -v maxnm={a.max_nm} -f {q(awk)} | gzip -1 > {q(pairs + '.part')}", lg)
    os.replace(pairs + ".part", pairs)
    shutil.rmtree(tmp, ignore_errors=True)
    if a.remove_bam:
        os.remove(bam)
    touch(os.path.join(Mp, "done.ok"))

# ---------------------------------------------------------------------- score ----
def read_bed(path):
    out = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            n, s, e = line.split("\t")
            out[n].append((int(s), int(e)))
    return {n: merge(iv) for n, iv in out.items()}

def scan(keys, wrong, own, a):
    '''Flag bins that contact one other sequence much more than the assembly's typical bin does.
    wrong[k]: Counter of contacts per other sequence; own[k]: cis contacts >= --far away.
    Returns (baseline ratio, bins tested, [(key, other, contacts with other, own)]).'''
    vals = []
    for k in keys:
        other, t = max(wrong[k].items(), key=lambda x: x[1]) if wrong.get(k) else (None, 0)
        if t + own[k] >= a.min_hom:
            vals.append((k, other, t, own[k]))
    if not vals:
        return 0.0, 0, []
    ratios = sorted(t / max(o, 1) for k, x, t, o in vals)
    base = ratios[len(ratios) // 2]
    flagged = [v for v in vals if v[2] >= a.min_hom and v[2] > a.fold * base * max(v[3], 1)]
    return base, len(vals), flagged

def runs_of(flagged):
    '''Adjacent flagged bins of one backbone -> [name, first bin, last bin + 1, other, contacts, own]'''
    runs = []
    for (name, b), other, t, o in sorted(flagged, key=lambda v: v[0]):
        if runs and runs[-1][0] == name and runs[-1][2] == b:
            r = runs[-1]
            r[2] = b + 1
            r[3][other] += t
            r[4] += t
            r[5] += o
        else:
            runs.append([name, b, b + 1, collections.Counter({other: t}), t, o])
    return runs

def bin_intervals(keys, B, length):
    by = collections.defaultdict(list)
    for name, b in keys:
        by[name].append((b * B, min(length[name], (b + 1) * B)))
    return {n: merge(iv) for n, iv in by.items()}

def subtract(a, b):
    '''Merged intervals a minus merged intervals b.'''
    out, j = [], 0
    for s, e in a:
        while j < len(b) and b[j][1] <= s:
            j += 1
        k, cur = j, s
        while k < len(b) and b[k][0] < e:
            if b[k][0] > cur:
                out.append((cur, b[k][0]))
            cur = max(cur, b[k][1])
            k += 1
        if cur < e:
            out.append((cur, e))
    return out

CATS = [("C", "consistent_bp"), ("X", "conflict_bp"), ("D", "duplicated_bp"),
        ("U", "untestable_bp"), ("O", "off_backbone_bp")]

def headline(bp, denom):
    '''hic_score and the one-line breakdown from the category bp; O is whatever of denom is not on
    backbones (fragments, unplaced, or absent when --genome_size exceeds the assembly).'''
    bp = dict(bp)
    bp["off_backbone_bp"] = max(0, denom - sum(bp[k] for k in ("consistent_bp", "conflict_bp",
                                                                  "duplicated_bp", "untestable_bp")))
    pct = lambda x: f"{100 * x / denom:.2f}" if denom else "nan"
    line = (f"C:{pct(bp['consistent_bp'])}% X:{pct(bp['conflict_bp'])}% "
            f"[switch:{pct(bp['switch_bp'])}% misjoin:{pct(bp['misjoin_bp'])}%] "
            f"D:{pct(bp['duplicated_bp'])}% U:{pct(bp['untestable_bp'])}% O:{pct(bp['off_backbone_bp'])}% "
            f"n:{denom / 1e6:.1f}Mb")
    return pct(bp["consistent_bp"]), line, bp["off_backbone_bp"]

def score(a):
    P, Mp = os.path.join(a.outdir, "prepare"), os.path.join(a.outdir, "map")
    pairs = os.path.join(Mp, "pairs.tsv.gz")
    if not done(os.path.join(Mp, "done.ok")):
        die("run 'map' first")
    S = os.path.join(a.outdir, f"xtalk_q{a.mapq}")
    os.makedirs(S, exist_ok=True)

    labels = read_labels(a.outdir)
    byT = {r["target_id"]: r for r in labels}
    length = {r["target_id"]: r["length"] for r in labels}
    chromOrder = list(dict.fromkeys(r["chrom"] for r in labels if r["chrom"] != "."))
    backbones = [r["target_id"] for r in labels if r["class"] == "backbone"]
    homPath = os.path.join(P, "homology.v4.tsv")
    if not os.path.exists(homPath):
        die(f"{homPath} not found; rerun 'prepare' (it reuses the mapping target, so 'map' stays valid)")
    # homologous regions per ordered backbone pair: a contact between two backbones is a homolog
    # contact only when each end lies in the region of its backbone that is homologous to the other
    homIv = collections.defaultdict(list)
    with open(homPath) as fh:
        for line in fh:
            if not line.startswith("#"):
                b, s_, e_, o = line.rstrip("\n").split("\t")
                homIv[(b, o)].append((int(s_), int(e_)))
    homIv = {k: merge(iv) for k, iv in homIv.items()}
    homStarts = {k: [s_ for s_, e_ in iv] for k, iv in homIv.items()}
    homolog_of = collections.defaultdict(set)
    for b, o in homIv:
        homolog_of[b].add(o)

    def in_hom(b, o, x):  # x is 1-based
        st = homStarts.get((b, o))
        if not st:
            return False
        i = bisect.bisect_right(st, x - 1) - 1
        return i >= 0 and x - 1 < homIv[(b, o)][i][1]

    excess = read_bed(os.path.join(P, "excess.v4.bed"))
    # one phasing mask: testable and not an extra copy. Contacts in excess regions are left out of
    # numerator and denominator alike; their bp are reported as duplicated.
    testable = read_bed(os.path.join(P, "testable.v4.bed"))
    phase = {n: subtract(iv, excess.get(n, [])) for n, iv in testable.items()}
    pStarts = {n: [s for s, e in iv] for n, iv in phase.items()}
    pEnds = {n: [e for s, e in iv] for n, iv in phase.items()}

    def in_phase(n, x):  # x is 1-based
        st = pStarts.get(n)
        if not st:
            return False
        i = bisect.bisect_right(st, x - 1) - 1
        return i >= 0 and x - 1 < pEnds[n][i]

    B, FAR = a.bin, a.far
    c = collections.Counter()
    per = collections.defaultdict(collections.Counter)
    own = collections.Counter()      # far cis, phase mask applied (switch scan)
    ownAll = collections.Counter()   # far cis, no mask (misjoin scan)
    hom = collections.defaultdict(collections.Counter)    # bin -> partner backbone -> contacts
    trans = collections.defaultdict(collections.Counter)  # bin -> other chromosome -> contacts
    mapPairs = gzip.open(os.path.join(S, "pairs.map.tsv.gz"), "wt", compresslevel=1) if a.contact_map else None

    with gzip.open(pairs, "rt") as fh:
        for line in fh:
            ca, xa, qa, cb, xb, qb = line.rstrip("\n").split("\t")
            c["pairs_in"] += 1
            if int(qa) < a.mapq or int(qb) < a.mapq:
                c["dropped_mapq"] += 1
                continue
            if ca.startswith("decoy_") or cb.startswith("decoy_"):
                c["dropped_decoy"] += 1
                continue
            ra, rb = byT.get(ca), byT.get(cb)
            if ra is None or rb is None:
                die(f"pairs.tsv.gz has sequence {ca if ra is None else cb} that is not in labels.v2.tsv")
            if ra["class"] == "unplaced" or rb["class"] == "unplaced":
                c["dropped_unplaced"] += 1
                continue
            xa, xb = int(xa), int(xb)
            if mapPairs:  # maps: every chromosome-assigned pair, unmasked
                mapPairs.write(f"{ca}\t{xa}\t{cb}\t{xb}\n")
            if ra["class"] != "backbone" or rb["class"] != "backbone":
                c["fragment_pairs"] += 1  # fragments are drawn in the maps, never scored
                continue
            c["backbone_pairs"] += 1
            ka, kb = (ca, xa // B), (cb, xb // B)
            homologous = ca != cb and in_hom(ca, cb, xa) and in_hom(cb, ca, xb)
            if ra["chrom"] != rb["chrom"] and not homologous:
                # structure, not phasing: no mask, so misjoined segments stay visible. Homologous
                # regions with different labels (a translocation in one haplotype) are scored
                # as homologs below instead.
                c["interchrom"] += 1
                per[ca]["interchrom"] += 1
                per[cb]["interchrom"] += 1
                trans[ka][rb["chrom"]] += 1
                trans[kb][ra["chrom"]] += 1
                continue
            far = ca == cb and abs(xa - xb) >= FAR
            if far:
                ownAll[ka] += 1
                ownAll[kb] += 1
            if ca != cb and not homologous:
                c["unresolved"] += 1  # e.g. two pieces of one chromosome copy: not penalised
                per[ca]["unresolved"] += 1
                per[cb]["unresolved"] += 1
                continue
            if not (in_phase(ca, xa) and in_phase(cb, xb)):
                c["dropped_masked"] += 1
                continue
            if ca != cb:
                c["homolog"] += 1
                per[ca]["homolog"] += 1
                per[cb]["homolog"] += 1
                hom[ka][cb] += 1
                hom[kb][ca] += 1
            else:
                c["cis"] += 1
                per[ca]["cis"] += 1
                if far:
                    c["cis_far"] += 1
                    own[ka] += 1
                    own[kb] += 1
    if mapPairs:
        mapPairs.close()

    # bin scans: switch (partner backbones, phase-testable bins) and misjoin (other chromosomes, all bins)
    allBins, phaseBins = [], []
    for name in backbones:
        L, iv = length[name], phase.get(name, [])
        for b in range(-(-L // B)):
            bs, be = b * B, min(L, (b + 1) * B)
            allBins.append((name, b))
            if total(intersect(iv, [(bs, be)])) >= a.min_testable * (be - bs):
                phaseBins.append((name, b))
    swBase, swTested, swFlag = scan(phaseBins, hom, own, a)
    mjBase, mjTested, mjFlag = scan(allBins, trans, ownAll, a)
    swRuns, mjRuns = runs_of(swFlag), runs_of(mjFlag)
    nm = lambda t: byT[t]["name"]
    for fname, runs, otherName in (("switch_runs.tsv", swRuns, nm), ("misjoin_runs.tsv", mjRuns, str)):
        with open(os.path.join(S, fname), "w") as fh:
            fh.write("#backbone\tstart\tend\tbins\ttop_other\tcontacts_with_other\town_far_contacts\n")
            for name, b0, b1, others, t, o in runs:
                fh.write(f"{nm(name)}\t{b0 * B}\t{min(length[name], b1 * B)}\t{b1 - b0}\t"
                         f"{otherName(others.most_common(1)[0][0])}\t{t}\t{o}\n")

    # every backbone bp in exactly one category; misjoin > switch > duplicated > consistent > untestable
    mjIv = bin_intervals([v[0] for v in mjFlag], B, length)
    swIv = bin_intervals([v[0] for v in swFlag], B, length)
    bp = collections.Counter()
    perBp = {}
    for t in backbones:
        mj = mjIv.get(t, [])
        sw = subtract(swIv.get(t, []), mj)
        bad = merge(mj + sw)
        x = {"misjoin_bp": total(mj), "switch_bp": total(sw),
             "consistent_bp": total(subtract(phase.get(t, []), bad)),
             "duplicated_bp": total(subtract(excess.get(t, []), bad))}
        x["conflict_bp"] = x["misjoin_bp"] + x["switch_bp"]
        x["untestable_bp"] = length[t] - x["conflict_bp"] - x["consistent_bp"] - x["duplicated_bp"]
        perBp[t] = x
        bp.update(x)

    def frac(x, y):
        return f"{x / y:.5f}" if y else "nan"
    with open(os.path.join(S, "per_backbone.tsv"), "w") as fh:
        fh.write("#backbone\tchrom\tlength\tconsistent_bp\tswitch_bp\tmisjoin_bp\tduplicated_bp\tuntestable_bp\t"
                 "partner_backbones\tcis\thomolog\tinterchrom\tunresolved\thomolog_fraction\tsource_id\tnote\n")
        for t in backbones:
            r, g, x = byT[t], per[t], perBp[t]
            pb = sorted(nm(p) for p in homolog_of[t])
            fh.write(f"{r['name']}\t{r['chrom']}\t{r['length']}\t{x['consistent_bp']}\t{x['switch_bp']}\t"
                     f"{x['misjoin_bp']}\t{x['duplicated_bp']}\t{x['untestable_bp']}\t{','.join(pb) or '.'}\t"
                     f"{g['cis']}\t{g['homolog']}\t{g['interchrom']}\t{g['unresolved']}\t"
                     f"{frac(g['homolog'], g['cis'] + g['homolog'])}\t{r['source_id']}\t{r['note']}\n")

    prep = dict(l.rstrip("\n").split("\t", 1) for l in open(os.path.join(P, "summary.v2.tsv")))
    allBp = int(prep["all_bp"])
    denom = a.genome_size or allBp
    hscore, hline, offBp = headline(bp, denom)
    summary = [
        # headline: every bp of the assembly (or of --genome_size) in one category
        ("hic_score", hscore), ("hic_breakdown", hline), ("score_denominator_bp", denom),
        ("consistent_bp", bp["consistent_bp"]), ("conflict_bp", bp["conflict_bp"]),
        ("switch_bp", bp["switch_bp"]), ("misjoin_bp", bp["misjoin_bp"]),
        ("duplicated_bp", bp["duplicated_bp"]), ("untestable_bp", bp["untestable_bp"]),
        ("off_backbone_bp", offBp),
        # error rates and events
        ("homolog_fraction", frac(c["homolog"], c["cis"] + c["homolog"])),
        ("interchrom_fraction", frac(c["interchrom"], c["backbone_pairs"])),
        ("switch_runs", len(swRuns)), ("misjoin_runs", len(mjRuns)),
        # assembly composition
        ("all_bp", allBp), ("backbones", prep["backbones"]), ("backbone_bp", prep["backbone_bp"]),
        ("fragments", prep["fragments"]), ("fragment_bp", prep["fragment_bp"]),
        ("unplaced", prep["unplaced"]), ("unplaced_bp", prep["unplaced_bp"]),
        ("chroms_without_backbone", prep["chroms_without_backbone"]),
        ("backbones_without_partner", prep["backbones_without_partner"]),
        ("possible_chimeras", prep["possible_chimeras"]),
        # diagnostics
        ("switch_baseline", f"{swBase:.4f}"), ("switch_bins_tested", swTested),
        ("misjoin_baseline", f"{mjBase:.4f}"), ("misjoin_bins_tested", mjTested),
        ("mapq_min", a.mapq), ("pairs_in", c["pairs_in"]), ("dropped_mapq", c["dropped_mapq"]),
        ("dropped_decoy", c["dropped_decoy"]), ("dropped_unplaced", c["dropped_unplaced"]),
        ("fragment_pairs", c["fragment_pairs"]), ("backbone_pairs", c["backbone_pairs"]),
        ("interchrom", c["interchrom"]), ("unresolved", c["unresolved"]),
        ("dropped_masked", c["dropped_masked"]), ("cis", c["cis"]), ("cis_far", c["cis_far"]),
        ("homolog", c["homolog"]),
    ]
    with open(os.path.join(S, "summary.tsv"), "w") as fh:
        for k, v in summary:
            fh.write(f"{k}\t{v}\n")
    log(f"score (MAPQ >= {a.mapq}): hic_score {hscore} ({hline}); homolog_fraction "
        f"{frac(c['homolog'], c['cis'] + c['homolog'])}; interchrom_fraction "
        f"{frac(c['interchrom'], c['backbone_pairs'])}; {len(swRuns)} switch run(s), {len(mjRuns)} misjoin run(s)")

    if a.contact_map:
        # AGP turns mapping-target IDs into display names (they differ when a target was reused)
        agp = os.path.join(S, "names.agp")
        with open(agp, "w") as fh:
            for r in labels:
                if r["class"] != "unplaced":
                    L = r["length"]
                    fh.write(f"{r['name']}\t1\t{L}\t1\tW\t{r['target_id']}\t1\t{L}\t+\n")
        grp = os.path.join(S, "groups.tsv")
        with open(grp, "w") as fh:
            for chrom in chromOrder:
                bbs = sorted((r for r in labels if r["chrom"] == chrom and r["class"] == "backbone"),
                             key=lambda r: -r["length"])
                frs = sorted((r for r in labels if r["chrom"] == chrom and r["class"] == "fragment"),
                             key=lambda r: (int(r["ref_start"]), -r["length"]))
                fh.write(f"{chrom}\t{','.join(r['name'] for r in bbs + frs)}\n")
            if a.genome_map:
                fh.write("all_backbones\t" + ",".join(byT[t]["name"] for t in backbones) + "\n")
        mapsDir = os.path.join(S, "maps")
        os.makedirs(mapsDir, exist_ok=True)
        run([a.python, a.contact_map, "--pairs", os.path.join(S, "pairs.map.tsv.gz"), "--agp", agp,
             "--groups", grp, "--bin", str(a.map_bin), "--out", mapsDir], os.path.join(S, "contact_map.log"))

# -------------------------------------------------------------------- compare ----
def compare(a):
    cols, table = [], collections.OrderedDict()
    base = [os.path.basename(os.path.normpath(d)) for d in a.dirs]
    for d, b in zip(a.dirs, base):
        f = os.path.join(d, f"xtalk_q{a.mapq}", "summary.tsv")
        if not os.path.exists(f):
            die(f"{f} not found; run 'score --mapq {a.mapq}' there")
        # folder name as column, or the path when two folders share a name
        name = b if base.count(b) == 1 else os.path.normpath(d)
        cols.append(name)
        vals = dict(line.rstrip("\n").split("\t", 1) for line in open(f))
        if a.genome_size:  # one denominator for every assembly
            bp = {k: int(vals[k]) for k in ("consistent_bp", "conflict_bp", "switch_bp", "misjoin_bp",
                                            "duplicated_bp", "untestable_bp")}
            vals["hic_score"], vals["hic_breakdown"], vals["off_backbone_bp"] = headline(bp, a.genome_size)
            vals["score_denominator_bp"] = a.genome_size
        for k, v in vals.items():
            table.setdefault(k, {})[name] = str(v)
    out = sys.stdout if a.out is None else open(a.out, "w")
    out.write("#metric\t" + "\t".join(cols) + "\n")
    for k, vals in table.items():
        out.write(k + "\t" + "\t".join(vals.get(cn, ".") for cn in cols) + "\n")
    if a.out:
        out.close()

# ----------------------------------------------------------------------- main ----
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("-o", "--outdir", required=True, help="output directory for this assembly")
        p.add_argument("-t", "--threads", type=int, default=8)
        p.add_argument("--ploidy", type=int, default=2)

    def prep_args(p):
        p.add_argument("--reference", help="FASTA whose --ref_regex sequences name the chromosomes")
        p.add_argument("--ref_regex", default=r"(.+)_h1",
                       help="selects one haplotype of --reference; group 1 = chromosome name [%(default)s]")
        p.add_argument("--input", action="append",
                       help="[hapN=|backbone=]FASTA; repeat for each file of the assembly. Labels are only "
                            "recorded as claims (keep them as in earlier runs so their mapping is reused)")
        p.add_argument("--backbone_regex", default=r"(.+)_h(\d+)",
                       help="IDs of a backbone= file: group 1 = claimed chromosome, group 2 = haplotype [%(default)s]")
        p.add_argument("--decoys", nargs="+", help="FASTA(s) of organelles / set-aside sequence")
        p.add_argument("--annotarium", help="path to annotarium.py")
        p.add_argument("--minimap2", default=None)
        p.add_argument("--preset", default="asm20", choices=["asm5", "asm10", "asm20"])
        p.add_argument("--phase_min_aln", type=int, default=10000)
        p.add_argument("--scaffold_len", type=int, default=5000000,
                       help="sequences this long are backbones (chromosome-scale) [%(default)s]")
        p.add_argument("--min_dominance", type=float, default=0.9)
        p.add_argument("--min_aligned_frac", type=float, default=0.5)
        p.add_argument("--min_part_len", type=int, default=50000)
        p.add_argument("--mask_min_aln", type=int, default=10000,
                       help="minimum alignment block for homology between sequences [%(default)s]")
        p.add_argument("--partner_frac", type=float, default=0.2,
                       help="a homologous chain also counts when it aligns this fraction of the other "
                            "sequence (lets short fragments count) [%(default)s]")
        p.add_argument("--min_homolog_region", type=int, default=200000,
                       help="aligned bp a chain of blocks needs to count as a homologous region [%(default)s]")
        p.add_argument("--homolog_gap", type=int, default=100000,
                       help="largest gap, on either sequence, between chained alignment blocks [%(default)s]")

    def map_args(p):
        p.add_argument("-1", dest="r1", help="Hi-C R1 FASTQ(s), comma-separated")
        p.add_argument("-2", dest="r2", help="Hi-C R2 FASTQ(s), comma-separated")
        p.add_argument("--aligner", default="bwa-mem2", choices=["bwa-mem2", "bwa"])
        p.add_argument("--sort_mem", default="1G")
        p.add_argument("--max_nm", type=int, default=3)
        p.add_argument("--remove_bam", action="store_true", help="delete the name-sorted BAM once pairs exist")

    def score_args(p):
        p.add_argument("--mapq", type=int, default=20)
        p.add_argument("--bin", type=int, default=100000, help="bin size of the switch/misjoin scans [%(default)s]")
        p.add_argument("--far", type=int, default=2000000,
                       help="cis contacts this far apart are the scans' 'own' signal [%(default)s]")
        p.add_argument("--min_hom", type=int, default=30,
                       help="a bin is tested with >= this many contacts, flagged with >= this many "
                            "contacts with one other sequence [%(default)s]")
        p.add_argument("--fold", type=float, default=4.0,
                       help="flag a bin when its other/own ratio is > this many times the assembly's "
                            "median bin ratio [%(default)s]")
        p.add_argument("--min_testable", type=float, default=0.5,
                       help="fraction of a bin that must be phase-testable for the switch scan [%(default)s]")
        p.add_argument("--genome_size", type=int, default=None,
                       help="expected bp of all haplotypes, the hic_score denominator (default: the "
                            "assembly's own size); give the same value for every assembly")
        p.add_argument("--contact_map", help="agentssembler scripts/contact_map.py, to also draw maps")
        p.add_argument("--map_bin", type=int, default=250000)
        p.add_argument("--no_genome_map", dest="genome_map", action="store_false",
                       help="skip the all-backbones map")
        p.add_argument("--python", default=sys.executable, help="python with numpy+matplotlib for maps")

    p = sub.add_parser("prepare"); common(p); prep_args(p)
    p = sub.add_parser("map"); common(p); map_args(p)
    p = sub.add_parser("score"); common(p); score_args(p)
    p = sub.add_parser("all"); common(p); prep_args(p); map_args(p); score_args(p)
    p = sub.add_parser("compare")
    p.add_argument("dirs", nargs="+")
    p.add_argument("--mapq", type=int, default=20)
    p.add_argument("--genome_size", type=int, default=None,
                   help="recompute hic_score with this one denominator for every assembly")
    p.add_argument("--out", default=None)
    a = ap.parse_args()

    if a.cmd == "compare":
        compare(a)
        return
    os.makedirs(a.outdir, exist_ok=True)
    if a.cmd in ("prepare", "all"):
        prepare(a)
    if a.cmd in ("map", "all"):
        map_reads(a)
    if a.cmd in ("score", "all"):
        score(a)

if __name__ == "__main__":
    main()
