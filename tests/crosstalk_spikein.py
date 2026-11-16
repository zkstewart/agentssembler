#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
crosstalk_spikein.py -- does utilities/hic_crosstalk.py find errors that are known to be there?

The test puts known errors into an assembly, scores each spiked copy with hic_crosstalk.py, and
checks that every category of its score moves as it should, and only where it should.

Spiked copies ("mutants"), each in its own FASTA, with at most one event per sequence:

  neutral      Every second sequence reverse-complemented, all renamed, order reversed.
               Nothing is wrong, so the score must not move.
  switch       Homologous blocks exchanged between two haplotypes of a chromosome (one block
               size per chromosome), and one whole arm exchanged.     -> X [switch]
  misjoin      Segments moved from one chromosome to another (one size per pair).
                                                                       -> X [misjoin]
               Also one segment inverted in place. Order and orientation inside a chromosome
               are not scored by any category, so this event is reported, never required.
  composition  One segment duplicated as an extra sequence             -> D
               one sequence cut into pieces shorter than a backbone    -> O, and its homolog U
               one sequence removed                                    -> O, and its homolog U

Subcommands

  simulate  A small genome with known truth: random chromosomes, --ploidy haplotypes differing
            by SNPs, and Hi-C read pairs drawn within a haplotype with 1/distance decay, plus
            pairs between homologs and between chromosomes as noise. Sequences are named
            chrN_hK, so the genome is its own reference for hic_crosstalk.py (--ref_regex).
  mutate    Writes the four mutants of any assembly (simulated or real) and truth.tsv.
            Homologous breakpoints for the switches come from a minimap2 alignment of the two
            haplotypes, so real assemblies with indels between haplotypes are handled.
  evaluate  Reads the hic_crosstalk.py output of the baseline and of each mutant, with
            truth.tsv, and reports per event whether a flagged run lies at it, the flagged bp
            away from every event, and the change in each category. Exits 1 if a check fails.
  demo      simulate, mutate, hic_crosstalk.py on all five assemblies, evaluate. With the
            defaults this takes about half an hour on 16 threads and needs about 4 GB of disk.

On a real assembly (the stronger test, since repeats and a real Hi-C library are involved):

  crosstalk_spikein.py mutate --fasta final.fasta --out spike/
  # for X in baseline neutral switch misjoin composition, one job each, with the SAME
  # --reference (the unmodified assembly) and the SAME --genome_size for all five:
  hic_crosstalk.py all -o xt/X --reference final.fasta --input spike/X.fa ... --genome_size G
  crosstalk_spikein.py evaluate --truth spike/truth.tsv --genome_size G --run baseline xt/baseline \\
      --run neutral xt/neutral --run switch xt/switch --run misjoin xt/misjoin \\
      --run composition xt/composition --out spike/report.tsv
  (unplaced sequence can be given to every hic_crosstalk.py run as a second --input, unchanged)

No output file is ever overwritten. Python 3 standard library; mutate needs minimap2; demo needs
what hic_crosstalk.py needs (annotarium, minimap2, bwa-mem2 or bwa, samtools).
'''
import argparse, collections, math, os, random, re, shutil, subprocess, sys, time

COMP = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")
HERE = os.path.dirname(os.path.abspath(__file__))
CROSSTALK = os.path.join(HERE, "..", "utilities", "hic_crosstalk.py")
MUTANTS = ("neutral", "switch", "misjoin", "composition")

def log(msg):
    print(f"[{time.strftime('%F %T')}] {msg}", file=sys.stderr, flush=True)

def die(msg):
    sys.exit(f"ERROR: {msg}")

def revcomp(seq):
    return seq.translate(COMP)[::-1]

def read_fasta(path):
    name, chunks = None, []
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name is not None:
                    yield name, b"".join(chunks)
                name, chunks = line[1:].split()[0].decode(), []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, b"".join(chunks)

def write_fasta(path, records, width=80):
    with open(path, "xb") as fh:                      # "x": never overwrite
        for name, seq in records:
            fh.write(b">" + name.encode() + b"\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + b"\n")

# ------------------------------------------------------------------- simulate ----
def simulate(a):
    os.makedirs(a.out, exist_ok=True)
    rng = random.Random(a.seed)
    to_acgt = bytes(b"ACGT"[i & 3] for i in range(256))
    other = {65: b"CGT", 67: b"AGT", 71: b"ACT", 84: b"ACG"}
    haps = collections.OrderedDict()
    for c in range(1, a.chroms + 1):
        L = int(a.chrom_len * (1 + 0.05 * (c - 1)))   # unequal lengths, as in a real genome
        base = rng.randbytes(L).translate(to_acgt)
        for k in range(1, a.ploidy + 1):
            seq = bytearray(base)
            if k > 1:
                for p in rng.sample(range(L), int(L * a.divergence)):
                    seq[p] = rng.choice(other[seq[p]])
            haps[f"chr{c}_h{k}"] = bytes(seq)
    write_fasta(os.path.join(a.out, "genome.fa"), haps.items())
    log(f"genome: {len(haps)} sequences, {sum(map(len, haps.values())) / 1e6:.1f} Mb")

    names = list(haps)
    chromOf = {n: n.rsplit("_h", 1)[0] for n in names}
    rl, lo = a.read_len, math.log(1000)
    def read(name, pos):
        s = haps[name][pos:pos + rl]
        return revcomp(s) if rng.random() < 0.5 else s
    with open(os.path.join(a.out, "hic_R1.fa"), "xb") as f1, open(os.path.join(a.out, "hic_R2.fa"), "xb") as f2:
        for i in range(a.pairs):
            n1 = rng.choice(names)
            L1 = len(haps[n1]) - rl
            x = rng.randrange(L1)
            u = rng.random()
            if u < a.cis:                             # same molecule: 1/distance decay
                n2 = n1
                d = int(math.exp(rng.uniform(lo, math.log(max(x, L1 - x, 1001)))))
                sides = [p for p in (x + d, x - d) if 0 <= p < L1]
                y = rng.choice(sides) if sides else x
            else:
                if u < a.cis + a.homolog:             # noise between homologs
                    pool = [n for n in names if chromOf[n] == chromOf[n1] and n != n1]
                else:                                 # noise between chromosomes
                    pool = [n for n in names if chromOf[n] != chromOf[n1]]
                n2 = rng.choice(pool)
                y = rng.randrange(len(haps[n2]) - rl)
            f1.write(b">p%d\n%s\n" % (i, read(n1, x)))
            f2.write(b">p%d\n%s\n" % (i, read(n2, y)))
    log(f"Hi-C: {a.pairs} pairs of {rl} bp (cis {a.cis}, homolog noise {a.homolog}, "
        f"inter-chromosome noise {1 - a.cis - a.homolog:.2f})")

# --------------------------------------------------------------------- mutate ----
CIGAR = re.compile(rb"(\d+)([MIDNSHP=X])")

def align(minimap2, threads, qname, qseq, tname, tseq, tmp):
    '''minimap2 of one haplotype (query) against its homolog (target) -> '+' alignments, longest first'''
    qf, tf = os.path.join(tmp, "q.fa"), os.path.join(tmp, "t.fa")
    for f, n, s in ((qf, qname, qseq), (tf, tname, tseq)):
        if os.path.exists(f):
            os.remove(f)                              # scratch files of this script only
        write_fasta(f, [(n, s)])
    out = subprocess.run([minimap2, "-cx", "asm20", "--secondary=no", "-t", str(threads), tf, qf],
                         check=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL).stdout
    rows = []
    for line in out.splitlines():
        f = line.split(b"\t")
        if f[4] != b"+" or int(f[11]) < 20:
            continue
        cg = next((x[5:] for x in f[12:] if x.startswith(b"cg:Z:")), None)
        if cg:
            rows.append((int(f[2]), int(f[3]), int(f[7]), int(f[8]), cg))
    return sorted(rows, key=lambda r: r[0] - r[1])

def lift(rows, qpos):
    '''query position -> homologous target position, through the longest alignment covering it'''
    for qs, qe, ts, te, cg in rows:
        if not qs <= qpos < qe:
            continue
        q, t = qs, ts
        for m in CIGAR.finditer(cg):
            n, op = int(m.group(1)), m.group(2)
            if op in b"M=X":
                if q + n > qpos:
                    return t + (qpos - q)
                q += n; t += n
            elif op == b"I":
                if q + n > qpos:
                    return t
                q += n
            elif op in b"DN":
                t += n
        return te
    return None

def lift_near(rows, qpos, limit):
    '''as lift, moving right in 10 kb steps when the position is not aligned'''
    for step in range(200):
        p = qpos + step * 10000
        if p >= limit:
            break
        t = lift(rows, p)
        if t is not None:
            return p, t
    return None, None

def mutate(a):
    os.makedirs(a.out, exist_ok=True)
    for f in [m + ".fa" for m in ("baseline",) + MUTANTS] + ["truth.tsv"]:
        if os.path.exists(os.path.join(a.out, f)):
            die(f"{os.path.join(a.out, f)} exists; choose a new --out (nothing is overwritten)")
    minimap2 = a.minimap2 or shutil.which("minimap2") or die("minimap2 not found (--minimap2)")
    rng = random.Random(a.seed)
    seqs = collections.OrderedDict()
    for path in a.fasta:
        for name, seq in read_fasta(path):
            if name in seqs:
                die(f"duplicate sequence {name}")
            seqs[name] = seq
    nameRe = re.compile(a.name_regex)
    byChrom = collections.OrderedDict()
    for name, seq in seqs.items():
        m = nameRe.fullmatch(name)
        if m and len(seq) >= a.min_len:
            byChrom.setdefault(m.group(1), []).append(name)
    chroms = [c for c, v in byChrom.items() if len(v) >= 2]
    if len(chroms) < 3:
        die(f"need at least 3 chromosomes with two sequences >= --min_len matching --name_regex; found {len(chroms)}")
    log(f"{len(seqs)} sequences; {len(chroms)} chromosomes with at least two haplotypes of >= {a.min_len} bp")
    Mb = lambda x: int(float(x) * 1e6)
    truth = []            # mutant, event, type, size, sequence, start, end, expect
    def start_for(L, size):
        if size >= L:
            return None
        lo_, hi_ = int(0.25 * L), int(0.75 * L) - size
        return rng.randrange(lo_, hi_) if hi_ > lo_ else (L - size) // 2

    write_fasta(os.path.join(a.out, "baseline.fa"), seqs.items())

    # neutral: nothing is wrong
    out = []
    for i, (name, seq) in enumerate(seqs.items()):
        out.append((f"scaffold_{i + 1}", revcomp(seq) if i % 2 else seq))
    write_fasta(os.path.join(a.out, "neutral.fa"), out[::-1])

    # switch: one event per chromosome
    mut = dict(seqs)
    tmp = os.path.join(a.out, "tmp_align")
    os.makedirs(tmp, exist_ok=True)
    specs = a.switch_sizes.split(",")
    if len(specs) > len(chroms):
        log(f"WARNING: {len(specs)} switch sizes but {len(chroms)} chromosomes; using the first {len(chroms)}")
    for ev, (spec, chrom) in enumerate(zip(specs, chroms), 1):
        A, B = byChrom[chrom][:2]
        sa, sb = seqs[A], seqs[B]
        rows = align(minimap2, a.threads, A, sa, B, sb, tmp)
        arm = spec == "arm"
        size = None if arm else Mb(spec)
        ok = False
        for attempt in range(20):
            a1 = int(rng.uniform(0.55, 0.7) * len(sa)) if arm else start_for(len(sa), size)
            if a1 is None:
                break
            a1, b1 = lift_near(rows, a1, len(sa))
            if a1 is None:
                continue
            if arm:
                a2, b2 = len(sa), len(sb)
            else:
                a2, b2 = lift_near(rows, a1 + size, len(sa))
                if a2 is None:
                    continue
            if b2 > b1 and 0.8 <= (b2 - b1) / (a2 - a1) <= 1.25:
                ok = True
                break
        if not ok:
            log(f"WARNING: no collinear {spec} Mb block found between {A} and {B}; event skipped")
            continue
        mut[A] = sa[:a1] + sb[b1:b2] + sa[a2:]
        mut[B] = sb[:b1] + sa[a1:a2] + sb[b2:]
        truth.append(("switch", f"sw{ev}", "switch", a2 - a1, A, a1, a1 + (b2 - b1), "switch"))
        truth.append(("switch", f"sw{ev}", "switch", a2 - a1, B, b1, b1 + (a2 - a1), "switch"))
        log(f"switch sw{ev}: {A}:{a1}-{a2} <-> {B}:{b1}-{b2}" + (" (arm)" if arm else ""))
    shutil.rmtree(tmp)
    write_fasta(os.path.join(a.out, "switch.fa"), mut.items())

    # misjoin: segment j moves from chromosome j to chromosome j+1; then one inversion in place
    mut, used = dict(seqs), set()
    specs = a.misjoin_sizes.split(",")
    for ev, spec in enumerate(specs, 1):
        if ev >= len(chroms):
            log(f"WARNING: not enough chromosomes for misjoin size {spec}; skipped")
            continue
        size, h = Mb(spec), (ev - 1) % 2
        donor, recip = byChrom[chroms[ev - 1]][h], byChrom[chroms[ev]][h]
        s = start_for(len(seqs[donor]), size)
        if s is None:
            log(f"WARNING: {donor} is too short for a {spec} Mb segment; skipped")
            continue
        t = int(rng.uniform(0.3, 0.7) * len(seqs[recip]))
        seg = seqs[donor][s:s + size]
        mut[donor] = seqs[donor][:s] + seqs[donor][s + size:]
        mut[recip] = seqs[recip][:t] + seg + seqs[recip][t:]
        used.update((donor, recip))
        truth.append(("misjoin", f"mj{ev}", "misjoin_insert", size, recip, t, t + size, "misjoin"))
        truth.append(("misjoin", f"mj{ev}", "misjoin_excision", size, donor, s, s, "report"))
        log(f"misjoin mj{ev}: {donor}:{s}-{s + size} -> {recip}:{t}")
    free = [n for c in chroms for n in byChrom[c][:2] if n not in used]
    if free:
        name, size = free[-1], Mb(a.inversion_size)
        s = start_for(len(seqs[name]), size)
        if s is not None:
            mut[name] = seqs[name][:s] + revcomp(seqs[name][s:s + size]) + seqs[name][s + size:]
            truth.append(("misjoin", "inv1", "inversion", size, name, s, s + size, "report"))
            log(f"inversion inv1: {name}:{s}-{s + size} (not scored by design; reported only)")
    write_fasta(os.path.join(a.out, "misjoin.fa"), mut.items())

    # composition: duplication, fragmentation, removal, on three chromosomes
    mut = collections.OrderedDict(seqs)
    A, B = byChrom[chroms[0]][:2]
    size = max(int(1.2 * a.min_len), Mb(a.dup_size))      # long enough to be a backbone itself
    if len(seqs[A]) >= 1.5 * size:
        s = start_for(len(seqs[A]), size)
        mut[f"spike_dup_of_{A}"] = seqs[A][s:s + size]
        truth.append(("composition", "dup1", "duplication", size, A, s, s + size, "D"))
    else:
        log(f"WARNING: {A} is too short to duplicate {size} bp of it; duplication skipped")
    A, B = byChrom[chroms[1]][:2]
    piece = int(0.6 * a.min_len)
    parts = [seqs[A][i:i + piece] for i in range(0, len(seqs[A]), piece)]
    idx = list(mut).index(A)
    items = list(mut.items())
    items[idx:idx + 1] = [(f"{A}_part{i + 1}", p) for i, p in enumerate(parts)]
    mut = collections.OrderedDict(items)
    truth.append(("composition", "frag1", "fragmentation", len(seqs[A]), A, 0, len(seqs[A]), "O"))
    # the homologs lose their only partner backbone when there are two haplotypes; with three or
    # more, another backbone remains and the pieces still count as copies, so they stay testable
    if len(byChrom[chroms[1]]) == 2:
        truth.append(("composition", "frag1", "homolog_untestable", len(seqs[B]), B, 0, len(seqs[B]), "U"))
    A, B = byChrom[chroms[2]][:2]
    del mut[B]
    truth.append(("composition", "del1", "removal", len(seqs[B]), B, 0, len(seqs[B]), "O"))
    for other in byChrom[chroms[2]]:                  # a copy is missing for every haplotype left
        if other != B:
            truth.append(("composition", "del1", "homolog_untestable", len(seqs[other]), other, 0,
                          len(seqs[other]), "U"))
    write_fasta(os.path.join(a.out, "composition.fa"), mut.items())

    with open(os.path.join(a.out, "truth.tsv"), "x") as fh:
        fh.write("#mutant\tevent\ttype\tsize\tsequence\tstart\tend\texpect\n")
        for row in truth:
            fh.write("\t".join(map(str, row)) + "\n")
    total = sum(map(len, seqs.values()))
    log(f"wrote baseline.fa, {', '.join(m + '.fa' for m in MUTANTS)} and truth.tsv to {a.out}; "
        f"use --genome_size {total} for every run")

# ------------------------------------------------------------------- evaluate ----
def a_version():
    '''the score version of the hic_crosstalk.py beside this test, which names its output folder'''
    for line in open(CROSSTALK):
        if line.startswith("SCORE_VERSION"):
            return int(line.split("=")[1])
    die(f"SCORE_VERSION not found in {CROSSTALK}")

def read_run(d, mapq):
    S = os.path.join(d, f"xtalk_q{mapq}.v{a_version()}")
    summ = dict(l.rstrip("\n").split("\t", 1) for l in open(os.path.join(S, "summary.tsv")))
    runs = {}
    for kind in ("switch", "misjoin"):
        runs[kind] = collections.defaultdict(list)
        for line in open(os.path.join(S, f"{kind}_runs.tsv")):
            if not line.startswith("#"):
                f = line.split("\t")
                runs[kind][f[0]].append((int(f[1]), int(f[2])))
    labels = {}
    with open(os.path.join(d, "prepare", "labels.v2.tsv")) as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            labels[r["source_id"]] = r
    return summ, runs, labels

def categories(summ, denom):
    bp = {k: int(summ[k + "_bp"]) for k in ("consistent", "switch", "misjoin", "duplicated", "untestable")}
    bp["off_backbone"] = max(0, denom - sum(bp.values()))
    return bp

def evaluate(a):
    runs = collections.OrderedDict((n, read_run(d, a.mapq)) for n, d in a.run)
    if "baseline" not in runs:
        die("give the unmodified assembly as --run baseline DIR")
    truth = collections.defaultdict(list)
    with open(a.truth) as fh:
        header = fh.readline().lstrip("#").rstrip("\n").split("\t")
        for line in fh:
            r = dict(zip(header, line.rstrip("\n").split("\t")))
            for k in ("size", "start", "end"):
                r[k] = int(r[k])
            truth[r["mutant"]].append(r)
    denom = a.genome_size or int(runs["baseline"][0]["all_bp"])
    cats = {n: categories(s, denom) for n, (s, _, _) in runs.items()}
    pct = lambda x: 100 * x / denom
    out, checks = [], []
    w = out.append

    w(f"## Categories, % of {denom / 1e6:.1f} Mb (change from baseline in brackets)")
    w("run\tC\tX_switch\tX_misjoin\tD\tU\tO\tswitch_runs\tmisjoin_runs\tpairs_scored_per_Mb")
    keys = ("consistent", "switch", "misjoin", "duplicated", "untestable", "off_backbone")
    for n, (s, _, _) in runs.items():
        cells = [f"{pct(cats[n][k]):.2f}" + ("" if n == "baseline" else
                 f" ({pct(cats[n][k] - cats['baseline'][k]):+.2f})") for k in keys]
        scored = int(s["pairs_in"]) - int(s["dropped_mapq"])
        w("\t".join([n] + cells + [s["switch_runs"], s["misjoin_runs"], f"{scored / (denom / 1e6):.0f}"]))
    flagged = lambda n: cats[n]["switch"] + cats[n]["misjoin"]

    if a.expect_clean_baseline:
        checks.append(("baseline: nothing flagged in an error-free assembly", pct(flagged("baseline")) <= 0.1,
                       f"X = {pct(flagged('baseline')):.2f}%"))
    if "neutral" in runs:
        dC = pct(cats["neutral"]["consistent"] - cats["baseline"]["consistent"])
        dX = pct(flagged("neutral") - flagged("baseline"))
        checks.append(("neutral: score unchanged by renaming, reordering and reverse-complementing",
                       abs(dC) <= 0.5 and abs(dX) <= 0.1, f"C {dC:+.2f} pp, X {dX:+.2f} pp"))

    w("\n## Spiked events")
    w("mutant\tevent\ttype\tsize_Mb\tsequence\tbackbone\tstart\tend\texpected\tnearest_run_bp\tflagged_in_event_bp\tfound")
    for n in ("switch", "misjoin"):
        if n not in runs:
            continue
        s, rr, labels = runs[n]
        events, near = [], collections.defaultdict(list)   # backbone -> truth intervals
        for r in truth[n]:
            lab = labels.get(r["sequence"])
            if lab is None or lab["class"] != "backbone":
                w("\t".join(map(str, (n, r["event"], r["type"], f"{r['size'] / 1e6:.2f}", r["sequence"],
                                      "not_a_backbone", ".", ".", r["expect"], ".", ".", "no"))))
                events.append((r, False))
                continue
            L, st, en = int(lab["length"]), r["start"], r["end"]
            if lab["strand"] == "-":
                st, en = L - en, L - st
            near[lab["name"]].append((st, en))
            kinds = ("switch", "misjoin") if r["expect"] == "report" else (r["expect"],)
            dist, inside = None, 0
            for k in kinds:
                for b0, b1 in rr[k].get(lab["name"], []):
                    d = max(0, b0 - en, st - b1)
                    dist = d if dist is None else min(dist, d)
                    inside += max(0, min(b1, en) - max(b0, st))
            found = dist is not None and dist <= a.window
            events.append((r, found))
            w("\t".join(map(str, (n, r["event"], r["type"], f"{r['size'] / 1e6:.2f}", r["sequence"], lab["name"],
                                  st, en, r["expect"], "none" if dist is None else dist, inside,
                                  "yes" if found else "no"))))
        need = collections.defaultdict(list)          # event -> found, per haplotype it touches
        for r, f in events:
            if r["expect"] == n and r["size"] >= a.min_detect * 1e6:
                need[r["event"]].append(f)
        missed = sorted(e for e, fs in need.items() if not any(fs))
        both = sum(1 for fs in need.values() if all(fs))
        checks.append((f"{n}: every spiked {n} of >= {a.min_detect} Mb has a flagged {n} run at it", not missed,
                       f"{len(need) - len(missed)} of {len(need)} events found"
                       + (f" ({both} on every haplotype involved)" if n == "switch" else "")
                       + (f"; missed {', '.join(missed)}" if missed else "")))
        # flagged bp further than --far from every spiked event
        off = 0
        for k in ("switch", "misjoin"):
            for bb, ivs in rr[k].items():
                for b0, b1 in ivs:
                    if not any(b0 <= en + a.far and b1 >= st - a.far for st, en in near.get(bb, [])):
                        off += b1 - b0
        wrong = "misjoin" if n == "switch" else "switch"      # the other kind of conflict
        dW = pct(cats[n][wrong] - cats["baseline"][wrong])
        checks.append((f"{n}: spiked {n} events are not also reported as {wrong}", dW <= 0.1,
                       f"X [{wrong}] {dW:+.2f} pp"))
        allow = flagged("baseline") + max(1_000_000, int(0.002 * denom))
        checks.append((f"{n}: little is flagged away from the spiked events", off <= allow,
                       f"{off / 1e6:.2f} Mb flagged more than {a.far / 1e6:g} Mb from any event "
                       f"(baseline flags {flagged('baseline') / 1e6:.2f} Mb in all)"))

    if "composition" in runs:
        want = collections.Counter()
        for r in truth["composition"]:
            want[r["expect"]] += r["size"]
        d = {k: cats["composition"][k] - cats["baseline"][k] for k in keys}
        w("\n## Composition events: expected and observed change, Mb")
        w("category\texpected_at_least\tobserved")
        # O is what the fixed denominator leaves after the backbones, so the duplicate's own length,
        # which is extra backbone sequence, comes off it
        want["O"] -= want["D"]
        for cat, key, frac in (("D", "duplicated", 0.9), ("O", "off_backbone", 0.9), ("U", "untestable", 0.8)):
            if cat == "D" and not want["D"]:
                continue
            w(f"{cat}\t{frac * want[cat] / 1e6:.2f}\t{d[key] / 1e6:.2f}")
            checks.append((f"composition: {cat} rises by at least {frac:.0%} of the expected length",
                           d[key] >= frac * want[cat], f"{d[key] / 1e6:.2f} of {want[cat] / 1e6:.2f} Mb"))
        dX = pct(flagged("composition") - flagged("baseline"))
        checks.append(("composition: duplication, fragmentation and removal are not reported as conflicts",
                       abs(dX) <= 0.2, f"X {dX:+.2f} pp"))

    w("\n## Checks")
    for name, ok, detail in checks:
        w(f"{'PASS' if ok else 'FAIL'}\t{name}\t{detail}")
    text = "\n".join(out) + "\n"
    sys.stdout.write(text)
    if a.out:
        with open(a.out, "x") as fh:
            fh.write(text)
    if not all(ok for _, ok, _ in checks):
        sys.exit(1)

# ----------------------------------------------------------------------- demo ----
def demo(a):
    sim, spike, xt = (os.path.join(a.out, d) for d in ("sim", "spike", "xt"))
    if os.path.exists(a.out) and os.listdir(a.out):
        die(f"{a.out} is not empty; give a new directory (nothing is overwritten)")
    aligner = a.aligner or ("bwa-mem2" if shutil.which("bwa-mem2") else "bwa")
    a.out_root, a.out = a.out, sim
    simulate(a)
    a.fasta, a.out = [os.path.join(sim, "genome.fa")], spike
    mutate(a)
    G = sum(len(s) for _, s in read_fasta(os.path.join(sim, "genome.fa")))
    runs = []
    for name in ("baseline",) + MUTANTS:
        d = os.path.join(xt, name)
        log(f"hic_crosstalk.py all: {name}")
        cmd = [sys.executable, CROSSTALK, "all", "-o", d, "-t", str(a.threads), "--ploidy", str(a.ploidy),
               "--annotarium", a.annotarium, "--reference", os.path.join(sim, "genome.fa"),
               "--input", os.path.join(spike, name + ".fa"), "--aligner", aligner, "--remove_bam",
               "-1", os.path.join(sim, "hic_R1.fa"), "-2", os.path.join(sim, "hic_R2.fa"),
               "--genome_size", str(G)]
        os.makedirs(xt, exist_ok=True)
        with open(os.path.join(xt, name + ".log"), "x") as lg:
            if subprocess.run(cmd, stdout=lg, stderr=subprocess.STDOUT).returncode:
                die(f"hic_crosstalk.py failed on {name}; see {os.path.join(xt, name + '.log')}")
        runs.append([name, d])
    a.run, a.truth, a.genome_size = runs, os.path.join(spike, "truth.tsv"), G
    a.out, a.expect_clean_baseline = os.path.join(a.out_root, "report.tsv"), True
    evaluate(a)

# ----------------------------------------------------------------------- main ----
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def sim_args(p):
        p.add_argument("--chroms", type=int, default=5)
        p.add_argument("--chrom_len", type=int, default=12_000_000, help="length of chr1; each next is 5%% longer")
        p.add_argument("--ploidy", type=int, default=2)
        p.add_argument("--divergence", type=float, default=0.02, help="SNP rate of each further haplotype against h1")
        p.add_argument("--pairs", type=int, default=6_000_000)
        p.add_argument("--read_len", type=int, default=150)
        p.add_argument("--cis", type=float, default=0.85, help="share of pairs within one haplotype")
        p.add_argument("--homolog", type=float, default=0.03, help="share of pairs between homologs (noise)")
    def mut_args(p):
        p.add_argument("--name_regex", default=r"(.+)_h(\d+)", help="group 1 = chromosome; sequences of one "
                       "chromosome are its haplotypes [%(default)s]")
        p.add_argument("--min_len", type=int, default=5_000_000, help="hic_crosstalk.py's --scaffold_len")
        p.add_argument("--switch_sizes", default="0.2,0.5,1,3,arm", help="Mb, one per chromosome; 'arm' "
                       "exchanges everything beyond one breakpoint [%(default)s]")
        p.add_argument("--misjoin_sizes", default="0.2,0.5,1,3", help="Mb, one per chromosome pair [%(default)s]")
        p.add_argument("--inversion_size", default="3", help="Mb [%(default)s]")
        p.add_argument("--dup_size", default="6", help="Mb; at least 1.2 x --min_len is used [%(default)s]")
        p.add_argument("--minimap2", default=None)
    def eval_args(p):
        p.add_argument("--mapq", type=int, default=20, help="the --mapq hic_crosstalk.py scored with")
        p.add_argument("--window", type=int, default=200_000, help="an event is found when a flagged run of "
                       "the expected kind lies within this many bp of it [%(default)s]")
        p.add_argument("--far", type=int, default=2_000_000, help="flagged bp further than this from every "
                       "event count as off-target; hic_crosstalk.py's --far [%(default)s]")
        p.add_argument("--min_detect", type=float, default=0.5, help="events of at least this many Mb must be "
                       "found for the check to pass; smaller ones are reported [%(default)s]")
    def common(p):
        p.add_argument("--seed", type=int, default=1)
        p.add_argument("-t", "--threads", type=int, default=8)

    p = sub.add_parser("simulate"); p.add_argument("out"); sim_args(p); common(p)
    p = sub.add_parser("mutate")
    p.add_argument("--fasta", nargs="+", required=True, help="the chromosome-scale sequences of one assembly")
    p.add_argument("--out", required=True); mut_args(p); common(p)
    p = sub.add_parser("evaluate")
    p.add_argument("--truth", required=True)
    p.add_argument("--run", nargs=2, action="append", required=True, metavar=("NAME", "DIR"),
                   help="baseline, neutral, switch, misjoin or composition, and its hic_crosstalk.py -o directory")
    p.add_argument("--genome_size", type=int, default=None, help="the one given to every hic_crosstalk.py run")
    p.add_argument("--expect_clean_baseline", action="store_true", help="simulated data: the baseline has no errors")
    p.add_argument("--out", default=None, help="also write the report here"); eval_args(p)
    p = sub.add_parser("demo"); p.add_argument("out")
    p.add_argument("--annotarium", required=True, help="path to annotarium.py")
    p.add_argument("--aligner", default=None, choices=["bwa-mem2", "bwa"])
    sim_args(p); mut_args(p); eval_args(p); common(p)
    a = ap.parse_args()
    {"simulate": simulate, "mutate": mutate, "evaluate": evaluate, "demo": demo}[a.cmd](a)

if __name__ == "__main__":
    main()
