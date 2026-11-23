# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Driver: runs stages in order, stops at each decision gate, writes a report for the reader. See AGENTS.md."""
import argparse, datetime, fnmatch, glob, gzip, hashlib, json, os, shutil, subprocess, sys
import yaml
from .stages import STAGES, BY_ID, GUIDES
from . import report, compare

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULTS = dict(telomere_motif="CCCTAAA", telomere_window=20000, telomere_min_copies=50, link_bonus=20000, hops=12, restore_connectors=True, read_preset="lr:hq", raw_preset="lr:hq", hic_allele_rescue=False)
TOOLS = dict(minimap2="minimap2", samtools="samtools", seqkit="seqkit", bwa="bwa", samblaster="samblaster", filter_bam="filter_bam", juicer_tools="", compleasm="", plot_python="")

def now(): return datetime.datetime.now().isoformat(timespec="seconds")
def wd(a):
    w = os.path.abspath(a.work)
    if not os.path.exists(f"{w}/config.yaml"): sys.exit(f"{w} has no config.yaml: run `agentssembler init` first")
    return w
def base_dir(w):
    """Where `init` was run. Relative config paths resolve from here on every later command, whatever the current directory.
    Work directories made before this file existed fall back to the current directory, as before."""
    f = f"{w}/state/base_dir"
    return open(f).read().strip() if os.path.exists(f) else os.getcwd()
def cfg(w):
    c = yaml.safe_load(open(f"{w}/config.yaml")); b = base_dir(w)
    fix = lambda p: os.path.normpath(os.path.join(b, p)) if isinstance(p, str) and p and not os.path.isabs(p) else p
    i = c.get("inputs") or {}
    for k, v in i.items(): i[k] = [fix(p) for p in v] if isinstance(v, list) else fix(v)
    t = c.get("tools") or {}
    for k, v in t.items():
        if isinstance(v, str) and "/" in v: t[k] = fix(v)      # a bare name is looked up on PATH
    for sec, k in (("compleasm", "library"), ("compute", "guard_path")):
        if (c.get(sec) or {}).get(k): c[sec][k] = fix(c[sec][k])
    return c
def decisions(w):
    D = {}
    for f in sorted(glob.glob(f"{w}/decisions/D*.yaml")): D[os.path.basename(f)[:-5]] = yaml.safe_load(open(f))
    return D
def settings(w):
    """Config settings, then decision values on top: a decision file records why a setting has its value."""
    s = dict(DEFAULTS); c = cfg(w); s.update(c.get("settings") or {})
    for d in decisions(w).values(): s.update(d.get("set") or {})
    return s
SCALE_AT = 36e6      # mean chromosome length of the citrus genomes the fixed windows were chosen on: four builds measured 35.0 to 37.1 Mb
def length_scale(w, compute=False):
    """How long this sample's chromosomes are, relative to citrus. Seed length, chimera windows and bin sizes are chromosome-scale and are
    multiplied by it. Mean chromosome length is chain bp over chromosomes x ploidy, so it is known once stage 04 is done, and it is kept in
    state/length_scale.json from then on: a threshold must not move between stages of one build. `reset` to stage 04 or earlier clears it.
    Capped so that seeds are never longer than the chain N50, and floored at 0.02. `settings.length_scale: 1` pins the citrus values."""
    c = cfg(w); v = (c.get("settings") or {}).get("length_scale")
    if v not in (None, "auto"): return dict(scale=float(v), source="set in config.yaml")
    f = f"{w}/state/length_scale.json"
    if os.path.exists(f): return json.load(open(f))
    fai = f"{w}/08_final_chains/chains.fasta.fai"
    if not compute or not os.path.exists(fai): return None
    if os.path.exists(f"{w}/state/05_hic_chains.done"):      # a build already under way keeps the windows it began with
        d = dict(scale=1.0, source="fixed windows: stage 05 had run before the scale was measured"); json.dump(d, open(f, "w")); return d
    L = sorted((int(l.split("\t")[1]) for l in open(fai)), reverse=True); tot = sum(L); acc = 0
    for n50 in L:
        acc += n50
        if acc >= tot / 2: break
    S = tot / (int(c["chromosomes"]) * int(c["ploidy"])); by_len, by_n50 = S / SCALE_AT, n50 / 1e6
    d = dict(scale=max(round(min(by_len, by_n50), 2), 0.02), source="chains of stage 04", chain_bp=tot, mean_chromosome=round(S), chain_n50=n50,
             limited_by="chain N50" if by_n50 < by_len else "floor" if by_len < 0.02 else "mean chromosome length")
    json.dump(d, open(f, "w")); return d
MIN_WITHIN_OF_H = 0.13   # on the three citrus builds of ordinary Hi-C depth, the fixed 700 pairs/Mb was 0.125 to 0.149 of H
def hic_density(w, compute=False):
    """min_within, the within-chain Hi-C density below which a chain carries no placement information, as a share of what this library gives.
    Such chains are low because their reads rarely place uniquely, which is relative: a library four times deeper lifts satellite chains
    over any fixed count. H is the length-weighted median of within-chain pairs/Mb over chains of at least the seed length, from
    10_groups/chain_stats.tsv, so it is known once stage 06 is done, and it is kept in state/hic_density.json from then on. `reset` to
    stage 06 or earlier clears it. `settings.min_within: 700` pins the earlier fixed value."""
    v = settings(w).get("min_within")
    if v not in (None, "auto"): return dict(min_within=float(v), source="set in config.yaml or by a decision")
    f = f"{w}/state/hic_density.json"
    if os.path.exists(f): return json.load(open(f))
    st = f"{w}/10_groups/chain_stats.tsv"
    if not compute or not os.path.exists(st): return None
    if os.path.exists(f"{w}/state/07_group.done"):           # a build already under way keeps the threshold it grouped with
        d = dict(min_within=700.0, source="fixed value: stage 07 had run before the density was measured"); json.dump(d, open(f, "w")); return d
    ls = length_scale(w); seed = 1e6 * (ls["scale"] if ls else 1)
    R = [(float(r[2]), int(r[1])) for r in (l.split("\t") for l in open(st) if l[0] != "#" and l.strip())]
    use = sorted(r for r in R if r[1] >= seed) or sorted(R); tot = sum(n for _, n in use); acc = 0; H = 0.0
    for H, n in use:
        acc += n
        if acc >= tot / 2: break
    d = dict(min_within=float(max(round(MIN_WITHIN_OF_H * H), 1)), source="chain_stats.tsv of stage 06", hic_density=H, share=MIN_WITHIN_OF_H, chains_used=len(use), seed_len=round(seed))
    json.dump(d, open(f, "w")); return d
def raw_reads(c):
    """Raw basecalled read files: inputs.raw_reads, plus inputs.raw_bams, its earlier name."""
    i = c["inputs"]; return list(i.get("raw_reads") or []) + list(i.get("raw_bams") or [])
def read_format(p):
    """bam for BAM, CRAM or SAM (by extension); otherwise fasta or fastq from the first character, gzipped or not; None if neither."""
    if p.endswith((".bam", ".cram", ".sam")): return "bam"
    with open(p, "rb") as f: gz = f.read(2) == b"\x1f\x8b"
    with (gzip.open(p, "rb") if gz else open(p, "rb")) as f: c = f.read(1)
    return {b">": "fasta", b"@": "fastq"}.get(c)
def input_link(w, p, base, ext, old):
    """Path under WORK of the 00_inputs link to read file p: base plus ext[format], and .gz if gzipped. A missing or unrecognised file gets old,
    the name every input had before 0.3.2, as does a build already linked under it."""
    try:
        fmt = read_format(p)
        with open(p, "rb") as f: gz = f.read(2) == b"\x1f\x8b"
    except OSError: fmt, gz = None, False
    new = base + ext[fmt] + (".gz" if gz else "") if fmt in ext else old
    if not os.path.lexists(f"{w}/00_inputs/{new}") and os.path.lexists(f"{w}/00_inputs/{old}"): new = old
    return "00_inputs/" + new
def active(w, cond):
    if cond is None: return True
    if cond == "raw": return bool(raw_reads(cfg(w)))
    if cond == "collapsed": return str(settings(w).get("collapsed")).lower() == "true"
    raise ValueError(cond)
def stage_env(w):
    c = cfg(w); s = settings(w); i = c["inputs"]; t = dict(TOOLS); t.update(c.get("tools") or {}); comp = c.get("compute") or {}; bus = c.get("compleasm") or {}
    abspath = lambda p: os.path.abspath(p) if p else ""
    e = dict(os.environ); e.update(AS_WORK=w, AS_REPO=REPO, AS_SAMPLE=str(c["sample"]), AS_CHROMS=str(c["chromosomes"]), AS_PLOIDY=str(c["ploidy"]),
        AS_UNITIGS=abspath(i["unitigs_fasta"]), AS_GFA=abspath(i["unitigs_gfa"]), AS_READS=abspath(i["corrected_reads"]), AS_HIC_R1=abspath(i["hic_r1"]), AS_HIC_R2=abspath(i["hic_r2"]),
        AS_REFERENCE=abspath(i.get("reference")), AS_BASELINE_AGP=abspath(i.get("baseline_agp")), AS_RAW_READS=":".join(abspath(b) for b in raw_reads(c)),
        AS_THREADS=str(comp.get("threads", 16)), AS_MIN_FREE_GB=str(comp.get("min_free_gb", 40)), AS_GUARD_PATH=str(comp.get("guard_path", w)),
        AS_MOTIF=str(s["telomere_motif"]), AS_TEL_WINDOW=str(s["telomere_window"]), AS_TEL_MIN=str(s["telomere_min_copies"]), AS_LINK_BONUS=str(s["link_bonus"]), AS_HOPS=str(s["hops"]), AS_CONNECTORS=str(s["restore_connectors"]).lower(),
        AS_READ_PRESET=str(s["read_preset"]), AS_RAW_PRESET=str(s["raw_preset"]), AS_HIC_RESCUE=str(s["hic_allele_rescue"]).lower(),
        AS_READS_LINK=input_link(w, i["corrected_reads"], "ont_corrected_reads", dict(fasta=".fasta", fastq=".fastq"), "ont_corrected_reads.fasta"),
        AS_HIC_R1_LINK=input_link(w, i["hic_r1"], "hic_R1", dict(fasta=".fa", fastq=".fq"), "hic_R1.fq.gz"),
        AS_HIC_R2_LINK=input_link(w, i["hic_r2"], "hic_R2", dict(fasta=".fa", fastq=".fq"), "hic_R2.fq.gz"),
        AS_COLLAPSED=str(s.get("collapsed", "false")).lower(), AS_UNLOC_PURGE=str(s.get("unloc_purge", 0.2)), AS_HAP_FIRST=str(s.get("hap_first", c["ploidy"] == 2)).lower(),
        AS_LINEAGE=str(bus.get("lineage", "")), AS_BUSCO_LIB=abspath(bus.get("library")))
    for k, v in t.items(): e["AS_" + k.upper()] = str(v)
    if s.get("hap_depth"):
        h = float(s["hap_depth"]); r = lambda x: str(round(x))
        # Collapsed sequence can be carried by up to `ploidy` haplotypes, so every "above shared" bound sits that far up. At ploidy 2 these are the
        # multiples the citrus genomes were built with: 4.4, 3 and 2.68.
        P = max(int(c["ploidy"]), 2)
        e.update(AS_HAP=str(s["hap_depth"]), AS_REPEAT=str(s.get("repeat_depth") or round(max(4.4, P + 1.4) * h)), AS_TWO_X=r(1.5 * h), AS_MAX_DEPTH=r((P + 1) * h), AS_MIN_ALLELE=r(0.36 * h), AS_MAX_PAIR=r(2.68 * h),
                 AS_RESTORE_LO=r(1.54 * h), AS_RESTORE_HI=r((P + 0.68) * h), AS_ONE_LO=r(0.64 * h), AS_SHARED_LO=r(1.54 * h), AS_SHARED_HI=r((P + 0.68) * h))
    ls = length_scale(w, compute=os.path.exists(f"{w}/state/04_chains.done"))
    if ls:
        k = ls["scale"]; r = lambda x: str(max(int(round(x * k)), 1))
        e.update(AS_SCALE=str(k), AS_SEED_LEN=r(1e6), AS_BREAK_DIST=r(5e5), AS_BREAK_FLANK=r(1e6), AS_BREAK_MIN_SIDE=r(1e5), AS_MAP_BIN=r(2.5e5), AS_DOSAGE_BIN=r(5e5))
    hd = hic_density(w, compute=os.path.exists(f"{w}/state/06_evidence.done"))
    if hd: e["AS_MIN_WITHIN"] = f"{hd['min_within']:g}"
    return e
def done(w, sid): return os.path.exists(f"{w}/state/{sid}.done")
def pending(w):
    """First stage with an unmet gate: its stage is done, a decision is missing."""
    D = decisions(w)
    for st in STAGES:
        if not done(w, st.id): return None, []
        miss = [g for g, cond in st.gates if active(w, cond) and g not in D]
        if miss: return st, miss
    return None, []

def cmd_init(a):
    w = os.path.abspath(a.work)
    if os.path.exists(f"{w}/config.yaml"): sys.exit(f"{w}/config.yaml exists, so this work directory is already set up, and init will not replace it.\nTo change a setting, edit that file, then `agentssembler reset {w} STAGE` from the earliest stage the change affects.\nFor a different sample or different inputs, give init a new work directory.")
    c = yaml.safe_load(open(a.config))
    for k in ("sample", "chromosomes", "ploidy", "inputs"):
        if k not in c: sys.exit(f"config lacks {k}")
    for k in ("unitigs_fasta", "unitigs_gfa", "corrected_reads", "hic_r1", "hic_r2"):
        if not os.path.exists(c["inputs"].get(k, "")): sys.exit(f"input {k}: {c['inputs'].get(k)} not found")
    if os.path.getsize(c["inputs"]["corrected_reads"]) and read_format(c["inputs"]["corrected_reads"]) not in ("fasta", "fastq"):
        sys.exit(f"input corrected_reads: {c['inputs']['corrected_reads']} is not FASTA or FASTQ (gzipped or not). Convert a BAM with `samtools fastq`.")
    for p in raw_reads(c):
        if not os.path.exists(p): sys.exit(f"raw reads not found: {p}")
        if read_format(p) is None: sys.exit(f"raw reads: {p} is not BAM, CRAM, SAM, FASTQ or FASTA (gzipped or not)")
        if ":" in p: sys.exit(f"raw reads: {p}: a path with ':' cannot be passed to the stages")
    for d in ("state", "decisions", "edits", "reports", "logs"): os.makedirs(f"{w}/{d}", exist_ok=True)
    shutil.copyfile(a.config, f"{w}/config.yaml")
    open(f"{w}/state/base_dir", "w").write(os.getcwd() + "\n")
    if not os.path.exists(f"{w}/README.md"): shutil.copyfile(f"{REPO}/config/handover_templates/README.md", f"{w}/README.md")
    src = f"{REPO}/config/edits_templates"
    for f in os.listdir(src):
        if not os.path.exists(f"{w}/edits/{f}"): shutil.copyfile(f"{src}/{f}", f"{w}/edits/{f}")
    print(f"Work directory ready: {w}\nNext: agentssembler run {w}")

def cmd_status(a):
    w = wd(a); D = decisions(w)
    for st in STAGES:
        skip = st.when and not active(w, st.when)
        mark = "skip" if skip else "done" if done(w, st.id) else "todo"
        gates = " ".join(f"{g}{'+' if g in D else '?'}" for g, cond in st.gates if active(w, cond))
        print(f"{mark:5} {st.id:16} {st.title:38} {gates}")
    st, miss = pending(w)
    if st: print(f"\nWaiting on {', '.join(miss)} after {st.id}. Read reports/{st.id}.md and docs/decisions/{GUIDES[miss[0]]}.")

def run_stage(w, st):
    env = stage_env(w); sh = f"{REPO}/stages/{st.id}.sh"; logf = f"{w}/logs/{st.id}.log"
    print(f"== {st.id}: {st.title}  {('(' + st.hours + ')') if st.hours else ''}", flush=True)
    with open(logf, "a") as lf:
        lf.write(f"\n=== {now()}\n"); lf.flush()
        p = subprocess.Popen(["bash", sh], env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, cwd=w)
        for line in p.stdout: lf.write(line); lf.flush()
        rc = p.wait()
    if rc:
        tail = open(logf).read().splitlines()[-15:]
        sys.exit(f"{st.id} failed with exit {rc}. Last lines of {logf}:\n" + "\n".join(tail))
    if st.id == "04_chains": length_scale(w, compute=True)
    if st.id == "06_evidence": hic_density(w, compute=True)
    report.write(w, st.id, cfg(w), settings(w), lambda cond: active(w, cond))
    json.dump(dict(finished=now(), settings=settings(w)), open(f"{w}/state/{st.id}.done", "w"), default=str)

def check_inputs(w):
    """Inputs can vanish between stages when disk space is freed. Stop now, not hours later at the stage that needs them."""
    i = cfg(w)["inputs"]; miss = []
    for k, v in i.items():
        for p in (v if isinstance(v, list) else [v]):
            if not p: continue
            q = os.path.abspath(p)
            if not os.path.exists(q): miss.append(f"{k}: {q}")
            elif os.path.getsize(q) == 0: print(f"WARNING: input {k} is empty: {q}", file=sys.stderr)
    if miss: sys.exit("Inputs missing. Restore them, or correct the paths in config.yaml:\n  " + "\n  ".join(miss))

def cmd_run(a):
    w = wd(a); ids = [s.id for s in STAGES]; check_inputs(w)
    if a.to and a.to not in ids: sys.exit(f"unknown stage {a.to}")
    if a.only and a.only not in ids: sys.exit(f"unknown stage {a.only}")
    for st in STAGES:
        if a.to and ids.index(st.id) > ids.index(a.to): print(f"Stopped after {a.to}, as asked."); return
        if st.when and not active(w, st.when):
            continue
        if done(w, st.id):
            if a.only != st.id: continue
        if a.only and a.only != st.id: continue
        gst, miss = pending(w)
        if gst and not a.only:
            print(f"\nSTOP: decision needed: {', '.join(miss)}\nReport: {w}/reports/{gst.id}.md"); [print(f"Guide:  {REPO}/docs/decisions/{GUIDES[m]}") for m in miss]
            print(f"Record: agentssembler decide {w} {miss[0]} --by NAME --reason '...'"); return
        run_stage(w, st)
        print(f"   report: {w}/reports/{st.id}.md", flush=True)
        if a.only: break
    gst, miss = pending(w)
    if gst:
        print(f"\nSTOP: decision needed: {', '.join(miss)}\nReport: {w}/reports/{gst.id}.md"); [print(f"Guide:  {REPO}/docs/decisions/{GUIDES[m]}") for m in miss]
        print(f"Record: agentssembler decide {w} {miss[0]} --by NAME --reason '...'")
    elif all(done(w, s.id) or (s.when and not active(w, s.when)) for s in STAGES):
        print("All stages finished and all decisions recorded.")

def cmd_decide(a):
    w = wd(a)
    if a.decision not in GUIDES: sys.exit(f"unknown decision {a.decision}: {', '.join(GUIDES)}")
    if not a.reason.strip(): sys.exit("a decision needs a reason")
    vals = {}
    for x in a.set or []:
        k, v = x.split("=", 1); vals[k] = yaml.safe_load(v)
    if a.decision == "D01":
        s = dict(settings(w)); s.update(vals)
        for k in ("hap_depth", "repeat_depth", "collapsed"):
            if s.get(k) is None: sys.exit(f"D01 must set {k}: --set {k}=VALUE")
    rec = dict(decision=a.decision, by=a.by, when=now(), reason=a.reason, set=vals or None, edits={})
    for f in sorted(glob.glob(f"{w}/edits/*")):
        if os.path.getsize(f): rec["edits"][os.path.basename(f)] = hashlib.md5(open(f, "rb").read()).hexdigest()
    yaml.safe_dump(rec, open(f"{w}/decisions/{a.decision}.yaml", "w"), sort_keys=False)
    print(f"Recorded {a.decision}. Next: agentssembler run {w}")

def cmd_reset(a):
    """Later results rest on this stage, so they go too; decisions made from them become stale and must be taken again."""
    w = wd(a); ids = [s.id for s in STAGES]
    if a.stage not in BY_ID: sys.exit(f"unknown stage {a.stage}")
    later = STAGES[ids.index(a.stage):]; arch = f"{w}/archive/{datetime.datetime.now():%Y%m%d_%H%M%S}"
    for st in later:
        # a mapping to a target an earlier stage builds is stale once that stage is reset, --remap or not
        remap = a.remap or bool(st.target_from and ids.index(a.stage) <= ids.index(st.target_from))
        for pat in st.outputs + (st.keep if remap else []):
            for p in glob.glob(f"{w}/{pat}"):
                if st.keep and not remap and any(fnmatch.fnmatch(os.path.relpath(p, w), k) for k in st.keep): continue
                d = f"{arch}/{os.path.relpath(p, w)}"; os.makedirs(os.path.dirname(d), exist_ok=True); shutil.move(p, d)
        for f in (f"{w}/state/{st.id}.done", f"{w}/reports/{st.id}.md") + ((f"{w}/state/length_scale.json",) if st.id == "04_chains" else (f"{w}/state/hic_density.json",) if st.id == "06_evidence" else ()):
            if os.path.exists(f): os.remove(f)
        for g, _ in st.gates if st.id != a.stage or not a.keepDecision else []:
            f = f"{w}/decisions/{g}.yaml"
            if os.path.exists(f):
                os.makedirs(f"{arch}/decisions", exist_ok=True); shutil.move(f, f"{arch}/decisions/{g}.yaml")
    print(f"Reset from {a.stage}. Earlier outputs and stale decisions are in {arch}. Delete it when you are sure.")

def cmd_report(a):
    w = wd(a); ids = [a.stage] if a.stage else [s.id for s in STAGES if done(w, s.id)]
    for i in ids:
        f = f"{w}/reports/{i}.md"
        if os.path.exists(f): print(open(f).read())

def cmd_env(a):
    w = wd(a); e = stage_env(w)
    for k in sorted(e):
        if k.startswith("AS_"): print(f"export {k}={json.dumps(e[k])}")

def cmd_compare(a):
    wa, wb = os.path.abspath(a.work), os.path.abspath(a.other)
    for w in (wa, wb):
        if not os.path.exists(f"{w}/config.yaml"): sys.exit(f"{w} has no config.yaml")
    print(compare.table(wa, wb, str(cfg(wa)["sample"]), str(cfg(wb)["sample"])))

def main():
    p = argparse.ArgumentParser(prog="agentssembler", description=__doc__); sp = p.add_subparsers(dest="cmd", required=True)
    x = sp.add_parser("init", help="create a work directory from a config"); x.add_argument("config"); x.add_argument("work"); x.set_defaults(fn=cmd_init)
    x = sp.add_parser("status", help="stages and decisions"); x.add_argument("work"); x.set_defaults(fn=cmd_status)
    x = sp.add_parser("run", help="run until the next decision"); x.add_argument("work"); x.add_argument("--to", help="stop after this stage"); x.add_argument("--only", help="run just this stage, even if done"); x.set_defaults(fn=cmd_run)
    x = sp.add_parser("decide", help="record a decision and its reason"); x.add_argument("work"); x.add_argument("decision"); x.add_argument("--by", required=True); x.add_argument("--reason", required=True); x.add_argument("--set", action="append"); x.set_defaults(fn=cmd_decide)
    x = sp.add_parser("reset", help="undo a stage and all later ones, moving outputs to archive/"); x.add_argument("work"); x.add_argument("stage"); x.add_argument("--remap", action="store_true", help="also move aside costly mappings"); x.add_argument("--keepDecision", action="store_true", help="keep this stage's own decisions; later ones still go"); x.set_defaults(fn=cmd_reset)
    x = sp.add_parser("report", help="print reports"); x.add_argument("work"); x.add_argument("stage", nargs="?"); x.set_defaults(fn=cmd_report)
    x = sp.add_parser("env", help="print the environment a stage gets, to run a stage script by hand or on an HPC"); x.add_argument("work"); x.set_defaults(fn=cmd_env)
    x = sp.add_parser("compare", help="side-by-side metrics of two finished work directories, for D13"); x.add_argument("work"); x.add_argument("other"); x.set_defaults(fn=cmd_compare)
    a = p.parse_args(); a.fn(a)
