"""Driver: runs stages in order, stops at each decision gate, writes a report for the reader. See AGENTS.md."""
import argparse, datetime, fnmatch, glob, hashlib, json, os, shutil, subprocess, sys
import yaml
from .stages import STAGES, BY_ID, GUIDES
from . import report

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEFAULTS = dict(telomere_motif="CCCTAAA", min_within=700, link_bonus=20000, hops=8)
TOOLS = dict(minimap2="minimap2", samtools="samtools", seqkit="seqkit", bwa="bwa", samblaster="samblaster", filter_bam="filter_bam", juicer_tools="", compleasm="", plot_python="")

def now(): return datetime.datetime.now().isoformat(timespec="seconds")
def wd(a):
    w = os.path.abspath(a.work)
    if not os.path.exists(f"{w}/config.yaml"): sys.exit(f"{w} has no config.yaml: run `agentssembler init` first")
    return w
def cfg(w): return yaml.safe_load(open(f"{w}/config.yaml"))
def decisions(w):
    D = {}
    for f in sorted(glob.glob(f"{w}/decisions/D*.yaml")): D[os.path.basename(f)[:-5]] = yaml.safe_load(open(f))
    return D
def settings(w):
    """Config settings, then decision values on top: a decision file records why a setting has its value."""
    s = dict(DEFAULTS); c = cfg(w); s.update(c.get("settings") or {})
    for d in decisions(w).values(): s.update(d.get("set") or {})
    return s
def active(w, cond):
    if cond is None: return True
    if cond == "raw": return bool(cfg(w)["inputs"].get("raw_bams"))
    if cond == "collapsed": return str(settings(w).get("collapsed")).lower() == "true"
    raise ValueError(cond)
def stage_env(w):
    c = cfg(w); s = settings(w); i = c["inputs"]; t = dict(TOOLS); t.update(c.get("tools") or {}); comp = c.get("compute") or {}; bus = c.get("compleasm") or {}
    abspath = lambda p: os.path.abspath(p) if p else ""
    e = dict(os.environ); e.update(AS_WORK=w, AS_REPO=REPO, AS_SAMPLE=str(c["sample"]), AS_CHROMS=str(c["chromosomes"]), AS_PLOIDY=str(c["ploidy"]),
        AS_UNITIGS=abspath(i["unitigs_fasta"]), AS_GFA=abspath(i["unitigs_gfa"]), AS_READS=abspath(i["corrected_reads"]), AS_HIC_R1=abspath(i["hic_r1"]), AS_HIC_R2=abspath(i["hic_r2"]),
        AS_REFERENCE=abspath(i.get("reference")), AS_BASELINE_AGP=abspath(i.get("baseline_agp")), AS_RAW_BAMS=":".join(abspath(b) for b in (i.get("raw_bams") or [])),
        AS_THREADS=str(comp.get("threads", 16)), AS_MIN_FREE_GB=str(comp.get("min_free_gb", 40)), AS_GUARD_PATH=str(comp.get("guard_path", w)),
        AS_MOTIF=str(s["telomere_motif"]), AS_MIN_WITHIN=str(s["min_within"]), AS_LINK_BONUS=str(s["link_bonus"]), AS_HOPS=str(s["hops"]),
        AS_COLLAPSED=str(s.get("collapsed", "false")).lower(), AS_HAP_FIRST=str(s.get("hap_first", c["ploidy"] == 2)).lower(),
        AS_LINEAGE=str(bus.get("lineage", "")), AS_BUSCO_LIB=abspath(bus.get("library")))
    for k, v in t.items(): e["AS_" + k.upper()] = str(v)
    if s.get("hap_depth"):
        h = float(s["hap_depth"]); r = lambda x: str(round(x))
        e.update(AS_HAP=str(s["hap_depth"]), AS_REPEAT=str(s.get("repeat_depth") or round(4.4 * h)), AS_TWO_X=r(1.5 * h), AS_MAX_DEPTH=r(3 * h), AS_MIN_ALLELE=r(0.36 * h), AS_MAX_PAIR=r(2.68 * h),
                 AS_RESTORE_LO=r(1.54 * h), AS_RESTORE_HI=r(2.68 * h), AS_ONE_LO=r(0.64 * h), AS_SHARED_LO=r(1.54 * h), AS_SHARED_HI=r(2.68 * h))
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
    w = os.path.abspath(a.work); os.makedirs(w, exist_ok=True)
    for d in ("state", "decisions", "edits", "reports", "logs"): os.makedirs(f"{w}/{d}", exist_ok=True)
    c = yaml.safe_load(open(a.config))
    for k in ("sample", "chromosomes", "ploidy", "inputs"):
        if k not in c: sys.exit(f"config lacks {k}")
    for k in ("unitigs_fasta", "unitigs_gfa", "corrected_reads", "hic_r1", "hic_r2"):
        if not os.path.exists(c["inputs"].get(k, "")): sys.exit(f"input {k}: {c['inputs'].get(k)} not found")
    for p in c["inputs"].get("raw_bams") or []:
        if not os.path.exists(p): sys.exit(f"raw BAM not found: {p}")
    if os.path.exists(f"{w}/config.yaml") and not a.force: sys.exit(f"{w}/config.yaml exists: use --force to replace it")
    shutil.copyfile(a.config, f"{w}/config.yaml")
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
    report.write(w, st.id, cfg(w), settings(w), lambda cond: active(w, cond))
    json.dump(dict(finished=now(), settings=settings(w)), open(f"{w}/state/{st.id}.done", "w"), default=str)

def cmd_run(a):
    w = wd(a); ids = [s.id for s in STAGES]
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
        for pat in st.outputs + (st.keep if a.remap else []):
            for p in glob.glob(f"{w}/{pat}"):
                if st.keep and not a.remap and any(fnmatch.fnmatch(os.path.relpath(p, w), k) for k in st.keep): continue
                d = f"{arch}/{os.path.relpath(p, w)}"; os.makedirs(os.path.dirname(d), exist_ok=True); shutil.move(p, d)
        for f in (f"{w}/state/{st.id}.done", f"{w}/reports/{st.id}.md"):
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

def main():
    p = argparse.ArgumentParser(prog="agentssembler", description=__doc__); sp = p.add_subparsers(dest="cmd", required=True)
    x = sp.add_parser("init", help="create a work directory from a config"); x.add_argument("config"); x.add_argument("work"); x.add_argument("--force", action="store_true"); x.set_defaults(fn=cmd_init)
    x = sp.add_parser("status", help="stages and decisions"); x.add_argument("work"); x.set_defaults(fn=cmd_status)
    x = sp.add_parser("run", help="run until the next decision"); x.add_argument("work"); x.add_argument("--to", help="stop after this stage"); x.add_argument("--only", help="run just this stage, even if done"); x.set_defaults(fn=cmd_run)
    x = sp.add_parser("decide", help="record a decision and its reason"); x.add_argument("work"); x.add_argument("decision"); x.add_argument("--by", required=True); x.add_argument("--reason", required=True); x.add_argument("--set", action="append"); x.set_defaults(fn=cmd_decide)
    x = sp.add_parser("reset", help="undo a stage and all later ones, moving outputs to archive/"); x.add_argument("work"); x.add_argument("stage"); x.add_argument("--remap", action="store_true", help="also move aside costly mappings"); x.add_argument("--keepDecision", action="store_true", help="keep this stage's own decisions; later ones still go"); x.set_defaults(fn=cmd_reset)
    x = sp.add_parser("report", help="print reports"); x.add_argument("work"); x.add_argument("stage", nargs="?"); x.set_defaults(fn=cmd_report)
    x = sp.add_parser("env", help="print the environment a stage gets, to run a stage script by hand or on an HPC"); x.add_argument("work"); x.set_defaults(fn=cmd_env)
    a = p.parse_args(); a.fn(a)
