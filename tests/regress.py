#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Regression test on real data: seed a work directory with an earlier build's costly alignments and hand decisions, run the driver through every stage, and compare outputs with that build.
A spec (see tests/specs/) names the seeds, decisions and comparisons, with paths relative to --saved.
Usage: regress.py spec.yaml --saved BUILD_DIR --work NEW_DIR [--to STAGE] [--missingReadsOk]
       regress.py spec.yaml --saved BUILD_DIR --work NEW_DIR --compareOnly
--missingReadsOk: when the saved build's read files have since been removed, stand in an empty read file and drop raw BAMs. Stage 08 chooses bridges with the reads, so the saved 13_scaffold/bridges/choices.tsv is seeded, and its comparison then shows nothing. Use --to 11_build: stage 12 onward needs reads, and those comparisons report MISSING.
'''
import argparse, filecmp, gzip, os, shutil, subprocess, sys, yaml
HERE = os.path.dirname(os.path.abspath(__file__)); CLI = os.path.join(HERE, "..", "bin", "agentssembler")
ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("spec"); ap.add_argument("--saved", required=True); ap.add_argument("--work", required=True); ap.add_argument("--to"); ap.add_argument("--compareOnly", action="store_true"); ap.add_argument("--threads", type=int, default=16); ap.add_argument("--missingReadsOk", action="store_true")
a = ap.parse_args(); sp = yaml.safe_load(open(a.spec)); saved = os.path.abspath(a.saved); w = os.path.abspath(a.work)
def run(*x):
    r = subprocess.run([sys.executable, CLI, *x], capture_output=True, text=True)
    if r.returncode: sys.exit(f"agentssembler {x[0]} failed:\n{r.stdout}{r.stderr}")
    return r
def fix(p): return os.path.join(saved, p)
if not a.compareOnly:
    cfg = sp["config"]
    for k, v in list(cfg["inputs"].items()):
        cfg["inputs"][k] = [fix(x) for x in v] if isinstance(v, list) else (fix(v) if v else v)
    cfg.setdefault("compute", {})["threads"] = a.threads
    os.makedirs(w, exist_ok=True)
    if a.missingReadsOk:
        if not a.to: sys.exit("--missingReadsOk needs --to, since later stages read the reads")
        if not os.path.exists(cfg["inputs"]["corrected_reads"]):
            print("NOTE corrected reads are missing; using an empty stand-in:", cfg["inputs"]["corrected_reads"], flush=True)
            cfg["inputs"]["corrected_reads"] = os.path.join(w, "missing_reads.fasta"); open(cfg["inputs"]["corrected_reads"], "w").close()
        kept = [b for b in cfg["inputs"].get("raw_bams") or [] if os.path.exists(b)]
        if len(kept) != len(cfg["inputs"].get("raw_bams") or []): print("NOTE raw BAMs are missing and were dropped", flush=True)
        cfg["inputs"]["raw_bams"] = kept
    c = os.path.join(w, "regress.config.yaml"); yaml.safe_dump(cfg, open(c, "w"))
    if not os.path.exists(f"{w}/config.yaml"): run("init", c, w)
    if a.missingReadsOk: sp.setdefault("seed", {})["13_scaffold/bridges/choices.tsv"] = "13_scaffold/bridges/choices.tsv"
    for dst, src in (sp.get("seed") or {}).items():
        d = f"{w}/{dst}"; os.makedirs(os.path.dirname(d), exist_ok=True)
        if not os.path.lexists(d): os.symlink(fix(src), d)
    for dst, src in (sp.get("edits") or {}).items():
        shutil.copyfile(fix(src), f"{w}/edits/{dst}")
    for dst, text in (sp.get("edits_text") or {}).items():
        open(f"{w}/edits/{dst}", "w").write(text)
    for _ in range(40):
        r = subprocess.run([sys.executable, CLI, "run", w] + (["--to", a.to] if a.to else []), capture_output=True, text=True)
        print(r.stdout, flush=True)
        if r.returncode: print(r.stderr); sys.exit(1)
        if "decision needed" not in r.stdout: break
        ids = [l for l in r.stdout.splitlines() if l.startswith("STOP: decision needed")][-1].split(":")[-1].replace(" ", "").split(",")
        for d in ids:
            cmd = ["decide", w, d, "--by", "regression", "--reason", f"Replaying the saved build's decision: {sp.get('reasons', {}).get(d, 'as saved')}"]
            for k, v in (sp.get("decisions", {}).get(d) or {}).items(): cmd += ["--set", f"{k}={v}"]
            run(*cmd)
def norm(p, mode):
    op = gzip.open if p.endswith(".gz") else open
    with op(p, "rt") as f: L = [l for l in f if not l.startswith("##")]
    return sorted(L) if mode == "sorted" else L
ok = bad = 0
for row in sp["compare"]:
    new, old, mode = row[0], row[1], row[2] if len(row) > 2 else "exact"
    pn, po = f"{w}/{new}", fix(old)
    if not os.path.exists(pn) or not os.path.exists(po): print(f"MISSING  {new}  (new {os.path.exists(pn)}, saved {os.path.exists(po)})"); bad += 1; continue
    same = filecmp.cmp(pn, po, shallow=False) if mode == "exact" else norm(pn, mode) == norm(po, mode)
    if same: ok += 1; print(f"same     {new}")
    else:
        x, y = norm(pn, "sorted"), norm(po, "sorted"); sx, sy = set(x), set(y); bad += 1
        print(f"DIFFERS  {new}  lines new {len(x)} saved {len(y)}; only new {len(sx - sy)}, only saved {len(sy - sx)}")
print(f"\n{ok} same, {bad} different or missing")
sys.exit(1 if bad else 0)
