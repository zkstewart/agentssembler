#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''Free disk in a finished work directory: delete large files that no later stage reads and that can be remade, and remake them on demand (docs/CLEANUP.md).
list: what would go, and what stays and why. Deletes nothing.
check: as delete, rebuilding and checksumming, but deletes nothing.
delete: before deleting, rebuild every cheap file in scratch and require a byte-for-byte match, and checksum the rest. Then delete, record each file in WORK/CLEANED.tsv, and add a line to WORK/README.md.
restore: remake deleted files from CLEANED.tsv and check them against it. A Hi-C BAM is checked through the pairs remade from it.
Covers WORK and WORK/archive/*/; --noArchive leaves the archive alone. --rules limits to some rules (names as list prints them); --skip leaves some out.
Usage: clean_work.py {list|check|delete|restore} WORK [--rules a,b] [--skip a,b] [--noArchive] [--tmp DIR]'''
import argparse, datetime, glob, gzip, hashlib, json, os, shlex, shutil, subprocess, sys, tempfile
R = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); S = f"{R}/scripts"
LINEAGE_FILES = ("miniprot_output.gff", "translated_protein.fasta")
HIC_INDEX = ("targets.fasta", "targets.fasta.amb", "targets.fasta.ann", "targets.fasta.bwt", "targets.fasta.pac", "targets.fasta.sa")

def md5(p, chunk=1 << 24):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""): h.update(b)
    return h.hexdigest()
def md5_stream(cmd, env):
    h = hashlib.md5(); p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, env=env, executable="/bin/bash")
    for b in iter(lambda: p.stdout.read(1 << 24), b""): h.update(b)
    if p.wait(): sys.exit(f"failed: {cmd}")
    return h.hexdigest()
def gz_md5(p):
    h = hashlib.md5()
    with gzip.open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 24), b""): h.update(b)
    return h.hexdigest()
def sh(cmd, env):
    r = subprocess.run(cmd, shell=True, env=env, executable="/bin/bash")
    if r.returncode: sys.exit(f"failed ({r.returncode}): {cmd}")
q = shlex.quote

def stage_env(w):
    r = subprocess.run([sys.executable, f"{R}/bin/agentssembler", "env", w], capture_output=True, text=True)
    if r.returncode: sys.exit(f"agentssembler env failed:\n{r.stderr}")
    e = dict(os.environ)
    for l in r.stdout.splitlines():
        if l.startswith("export "): k, v = l[7:].split("=", 1); e[k] = json.loads(v)
    # as stages/lib.sh: the tools' folders on PATH, and compleasm's, where it finds miniprot and hmmsearch
    e["PATH"] = ":".join([os.path.dirname(e[k]) for k in ("AS_MINIMAP2", "AS_SAMTOOLS", "AS_SEQKIT") if e.get(k)] + [e["PATH"]] + ([os.path.dirname(e["AS_COMPLEASM"])] if e.get("AS_COMPLEASM") else []))
    return e
def guard(e):
    '''The stages' disk guard (compute.guard_path, compute.min_free_gb), before anything is remade.'''
    g, lim = e.get("AS_GUARD_PATH"), int(e.get("AS_MIN_FREE_GB") or 0)
    if g and lim and shutil.disk_usage(g).free / 2**30 < lim: sys.exit(f"Only {shutil.disk_usage(g).free / 2**30:.0f} GB free on {g} (compute.min_free_gb {lim}): not restoring")

# Each rule: the files it covers in one base directory (WORK, or an archive/<time>/ inside it), and how to remake each.
# A file is a candidate only if no later stage reads it, or (archive=True) only inside archive/, where no stage runs.
# kind: "exact" files are rebuilt in scratch and must match before deletion; "hic" and "compleasm" are checksummed only.
class Ctx:
    def __init__(s, w, base, env): s.w, s.b, s.e, s.sm = w, base, env, env["AS_SAMPLE"]
    def p(s, *x): return os.path.join(s.b, *x)
    def final_dir(s): return "15_final_raw" if os.path.exists(s.p("15_final_raw", f"{s.sm}.final.agp")) else "15_final"

def cat_rule(name, made, target, parts, archive=False):
    def find(c): return [target] if os.path.exists(c.p(target)) else []
    def needs(c, path): return [c.p(x.format(s=c.sm)) for x in parts]
    def build(c, path, out): sh(f"cat {' '.join(q(x) for x in needs(c, path))} > {q(out)}", c.e)
    return dict(name=name, made=made, kind="exact", archive=archive, find=find, needs=needs, build=build, how="cat " + " + ".join(parts))

def r_candidate():
    t = "17_raw_reads/fill/candidate.fasta"
    def needs(c, path): return [c.p("17_raw_reads/fill/candidate.agp"), c.w + "/00_inputs/unitigs.fasta", c.p("14_readfill/fills.fasta"), c.p("17_raw_reads/fill/fills.fasta")]
    def build(c, path, out):
        n = needs(c, path); sh(f"python3 {q(S)}/agp_to_fasta.py {q(n[0])} {q(','.join(n[1:]))} {q(out)}", c.e)
    return dict(name="candidate", made="stage 13 (agp_to_fasta.py); stage 14 reads candidate.agp, not this", kind="exact", archive=False,
                find=lambda c: [t] if os.path.exists(c.p(t)) else [], needs=needs, build=build, how="agp_to_fasta.py candidate.agp")

def r_subset(name, made, pattern, names_of, src, archive):
    def find(c): return [os.path.relpath(x, c.b) for x in sorted(glob.glob(c.p(pattern)))]
    def needs(c, path): return [c.p(names_of(path)), c.p(src)]
    def build(c, path, out):
        n = needs(c, path); sh(f"{q(c.e['AS_SEQKIT'])} grep -f {q(n[0])} {q(n[1])} -w 0 > {q(out)} 2>/dev/null", c.e)
    return dict(name=name, made=made, kind="exact", archive=archive, find=find, needs=needs, build=build, how=f"seqkit grep -f names.txt {src}")

def r_mask():
    pat = "14_readfill/target.h*.fasta"
    def needs(c, path): return [c.p("13_scaffold", f"{c.sm}.graph.agp"), c.p("13_scaffold/target.fasta")]
    def build(c, path, out):
        n = needs(c, path); h = path.split(".")[-2]; d = tempfile.mkdtemp(dir=os.path.dirname(out))
        try: sh(f"python3 {q(S)}/mask_shared.py {q(n[0])} {q(n[1])} {q(d)}/target > /dev/null", c.e); shutil.move(f"{d}/target.{h}.fasta", out)
        finally: shutil.rmtree(d, ignore_errors=True)
    return dict(name="readfill_masked", made="stage 12 (mask_shared.py, collapsed samples); stage 13 maps raw reads to it, so archive only", kind="exact", archive=True,
                find=lambda c: [os.path.relpath(x, c.b) for x in sorted(glob.glob(c.p(pat)))], needs=needs, build=build, how="mask_shared.py graph.agp target.fasta")

def r_bridge_reads():
    t = "13_scaffold/bridges/reads.fasta"
    def needs(c, path): return [c.p("13_scaffold/bridges/ids.txt"), c.e["AS_READS"]]
    def build(c, path, out):
        n = needs(c, path); k = q(c.e["AS_SEQKIT"]); sh(f"set -o pipefail; {k} grep -f {q(n[0])} {q(n[1])} 2>/dev/null | {k} fq2fa -w 0 -o {q(out)} 2>/dev/null", c.e)
    return dict(name="bridge_reads", made="stage 08 (choose_bridges.py), which deletes it itself; only older or archived runs have it", kind="exact", archive=False,
                find=lambda c: [t] if os.path.exists(c.p(t)) else [], needs=needs, build=build, how="seqkit grep -f ids.txt corrected reads | fq2fa", slow=True)

def hic_ok(c):
    '''The Hi-C mapping can go only if both pair tables exist, were made from these chains, and the reads and decoy list are still there.'''
    d = c.p("09_hic_map"); why = []
    for f in ("pairs.tsv.gz", "pairs.rescued.tsv.gz", "decoys.txt", "chains.fasta", "alleles08.pairs.tsv"):
        if not os.path.exists(f"{d}/{f}"): why.append(f"no 09_hic_map/{f}")
    for k in ("AS_HIC_R1", "AS_HIC_R2"):
        if not os.path.exists(c.e.get(k, "")): why.append(f"Hi-C reads missing: {c.e.get(k)}")
    if not why and os.path.exists(f"{d}/pairs.tsv.gz.target") and open(f"{d}/pairs.tsv.gz.target").read().strip() != md5(f"{d}/chains.fasta"):
        why.append("pairs.tsv.gz.target does not match chains.fasta")
    return why
def hic_targets(c, out):
    d = c.p("09_hic_map"); sh(f"cp {q(d)}/chains.fasta {q(out)}", c.e)
    if os.path.getsize(f"{d}/decoys.txt"): sh(f"{q(c.e['AS_SAMTOOLS'])} faidx {q(c.w)}/00_inputs/unitigs.fasta -r {q(d)}/decoys.txt | sed 's/^>/>decoy_/' >> {q(out)}", c.e)
def r_hic_targets():
    return dict(name="hic_targets", made="stage 05 (chains plus decoy unitigs); read only to map Hi-C, and pairs.tsv.gz is kept", kind="exact", archive=False,
                find=lambda c: ["09_hic_map/targets.fasta"] if os.path.exists(c.p("09_hic_map/targets.fasta")) else [], cond=hic_ok,
                needs=lambda c, path: [c.p("09_hic_map/chains.fasta"), c.p("09_hic_map/decoys.txt")], build=lambda c, path, out: hic_targets(c, out), how="chains.fasta + decoy_ unitigs")
def r_hic_map():
    def find(c): return [f"09_hic_map/{f}" for f in ("HiC.bam",) + HIC_INDEX[1:] if os.path.exists(c.p("09_hic_map", f))]
    return dict(name="hic_map", made="stage 05 (bwa index, bwa mem, samblaster); only pairs.tsv.gz and pairs.rescued.tsv.gz are read later", kind="hic", archive=False,
                find=find, cond=hic_ok, how="bwa index + bwa mem -5SP | samblaster | samtools view -F 3340; about 1-2 h")
def compleasm_runs(c):
    '''compleasm runs whose two big intermediates can go: their tables were kept, and their input is known.'''
    L = c.e.get("AS_LINEAGE", ""); runs = []
    for d in [c.p("01_genes")] + sorted(glob.glob(c.p("16_validate/compleasm/*"))):
        if os.path.isdir(f"{d}/{L}") and os.path.exists(f"{d}/{L}/full_table_busco_format.tsv") and compleasm_sources(c, os.path.relpath(d, c.b)):
            runs.append(os.path.relpath(d, c.b))
    return runs
def compleasm_sources(c, run):
    '''The files a compleasm run's input is made from, or None when the run is not one the stages make (it is then left alone).'''
    sm, F = c.sm, c.final_dir(); n = os.path.basename(run)
    if run == "01_genes": src = [f"{c.w}/00_inputs/unitigs.fasta"]
    elif n == "unplaced": src = [c.p(F, f"{sm}.final.unplaced.fasta")]
    elif n.startswith("h") and n[1:].isdigit(): src = [c.p(F, f"{sm}.final.fasta")]
    elif n.endswith("_unloc") and n[1:-6].isdigit(): src = [c.p(F, f"{sm}.final.fasta"), c.p("16_validate/unloc", f"{sm}.unloc.fasta")]
    else: return None
    return src if all(os.path.exists(x) for x in src) else None
def compleasm_input(c, run, out):
    sm, F = c.sm, c.final_dir(); n = os.path.basename(run); k = q(c.e["AS_SEQKIT"])
    if run == "01_genes": return f"{c.w}/00_inputs/unitigs.fasta"
    if n == "unplaced": return c.p(F, f"{sm}.final.unplaced.fasta")
    if n.startswith("h") and n[1:].isdigit(): sh(f"{k} grep -r -p {q('_' + n + '$')} {q(c.p(F, sm + '.final.fasta'))} > {q(out)} 2>/dev/null", c.e); return out
    if n.endswith("_unloc") and n[1:-6].isdigit():
        h = n[:-6]; sh(f"{{ {k} grep -r -p {q('_' + h + '$')} {q(c.p(F, sm + '.final.fasta'))}; {k} grep -r -p {q('_' + h + '_unloc_[0-9]+$')} {q(c.p('16_validate/unloc', sm + '.unloc.fasta'))}; }} > {q(out)} 2>/dev/null", c.e); return out
    return None
def r_compleasm():
    def find(c):
        L = c.e.get("AS_LINEAGE", "")
        return [f"{r}/{L}/{f}" for r in compleasm_runs(c) for f in LINEAGE_FILES if os.path.exists(c.p(r, L, f))]
    def cond(c, path=None):
        why = [] if os.path.exists(c.e.get("AS_COMPLEASM", "") or "/nonexistent") else ["compleasm not set or missing"]
        return why
    return dict(name="compleasm", made="compleasm (stages 01 and 15); stages read full_table*.tsv and summary.txt, which stay", kind="compleasm", archive=False,
                find=find, cond=cond, how="compleasm run on the same input in scratch; copy the two files back; about 15 min a run")

RULES = [r_hic_map(), r_hic_targets(), r_compleasm(),
         cat_rule("pre_target", "stage 10 (cat), which remakes it when missing", "12_shared/pre/target.fasta", ["12_shared/pre/{s}.fasta", "12_shared/pre/{s}.unplaced.fasta"]),
         cat_rule("scaffold_target", "stage 11 (cat); stages 12 and 13 read it, so archive only", "13_scaffold/target.fasta", ["13_scaffold/{s}.graph.fasta", "13_scaffold/{s}.graph.unplaced.fasta"], archive=True),
         cat_rule("rawfill_in", "stage 13 (cat), read only inside stage 13", "17_raw_reads/fill/in.fasta", ["15_final/{s}.final.fasta", "15_final/{s}.final.unplaced.fasta"]),
         r_subset("rawfill_hap_in", "stage 13 (seqkit, collapsed samples), read only inside stage 13", "17_raw_reads/fill/h*/in.fasta",
                  lambda p: os.path.join(os.path.dirname(p), "names.txt"), "17_raw_reads/fill/in.fasta", archive=False),
         r_subset("readfill_hap_in", "stage 12 (seqkit, collapsed samples); stage 13 reads it, so archive only", "14_readfill/h*/in.fasta",
                  lambda p: os.path.join(os.path.dirname(p), "names.txt"), "13_scaffold/target.fasta", archive=True),
         r_mask(), r_candidate(), r_bridge_reads()]
ORDER = ["scaffold_target", "pre_target", "rawfill_in", "readfill_hap_in", "rawfill_hap_in", "readfill_masked", "candidate", "bridge_reads", "hic_targets", "hic_map", "compleasm"]
BYNAME = {r["name"]: r for r in RULES}

def bases(w, noArchive):
    b = [(w, False)]
    if not noArchive: b += [(d, True) for d in sorted(glob.glob(f"{w}/archive/*/")) if os.path.isdir(d)]
    return [(os.path.normpath(x), a) for x, a in b]

def candidates(w, env, rules, noArchive):
    '''(ctx, rule, path rel to base, size) for files that may go; and (rel to WORK, why) for files kept.'''
    go, keep = [], []
    for base, isArch in bases(w, noArchive):
        c = Ctx(w, base, env)
        for r in rules:
            if r["archive"] and not isArch: continue
            files = r["find"](c)
            if not files: continue
            why = r.get("cond", lambda c: [])(c)
            for f in files:
                rel = os.path.relpath(c.p(f), w); miss = [x for x in (r["needs"](c, f) if "needs" in r else []) if not os.path.exists(x)]
                if os.path.islink(c.p(f)): keep.append((rel, "a symlink")); continue
                if why or miss: keep.append((rel, "; ".join(why + [f"missing {os.path.relpath(x, w) if x.startswith(w) else x}" for x in miss]))); continue
                go.append((c, r, f, os.path.getsize(c.p(f))))
    return go, keep

def G(n): return f"{n / 1e9:.2f} GB"

def cmd_list(w, env, rules, a):
    go, keep = candidates(w, env, rules, a.noArchive); tot = 0
    for name in ORDER:
        rows = [x for x in go if x[1]["name"] == name]
        if not rows: continue
        s = sum(x[3] for x in rows); tot += s; r = BYNAME[name]
        print(f"\n{name}: {len(rows)} files, {G(s)}. Made by {r['made']}. Remade by: {r['how']}{' (slow)' if r.get('slow') else ''}")
        for c, r, f, n in rows: print(f"  {G(n):>9}  {os.path.relpath(c.p(f), w)}")
    print(f"\nTotal: {len(go)} files, {G(tot)}")
    if keep: print("\nKept, although a rule covers them:"); [print(f"  {p}: {why}") for p, why in keep]
    return go

MAN = "CLEANED.tsv"; HDR = "#path\tbytes\tmd5\trule\tcheck\tdeleted\n"
def read_manifest(w):
    m = {}
    if os.path.exists(f"{w}/{MAN}"):
        for l in open(f"{w}/{MAN}"):
            if l[0] != "#": s = l.rstrip("\n").split("\t"); m[s[0]] = s
    return m

def cmd_check(w, env, rules, a): cmd_delete(w, env, rules, a, dry=True)
def cmd_delete(w, env, rules, a, dry=False):
    go = cmd_list(w, env, rules, a)
    if not go: return
    tmp = tempfile.mkdtemp(prefix="clean_", dir=a.tmp or w); rows = []; bad = []; when = datetime.date.today().isoformat()
    try:
        for c, r, f, n in go:
            path = c.p(f); rel = os.path.relpath(path, w); print(f"checking {rel}", flush=True); h = md5(path); check = "md5"
            if r["kind"] == "exact":
                out = f"{tmp}/rebuilt"; r["build"](c, f, out); hb = md5(out); os.remove(out)
                if hb != h: bad.append(rel); print(f"  DIFFERS when rebuilt: kept"); continue
                check = "rebuilt and matched"
            elif r["kind"] == "hic" and f.endswith("HiC.bam"):
                d = c.p("09_hic_map"); h = "-"; check = f"pairs={gz_md5(d + '/pairs.tsv.gz')};rescued={gz_md5(d + '/pairs.rescued.tsv.gz')}"
            elif r["kind"] == "compleasm": check = "md5, for information (compleasm is not deterministic)"
            rows.append((path, rel, n, h, r["name"], check))
    finally: shutil.rmtree(tmp, ignore_errors=True)
    if dry:
        print(f"\nCheck only: {len(rows)} files would go, {G(sum(x[2] for x in rows))}" + (f"; {len(bad)} did not rebuild identically and would stay: {', '.join(bad)}" if bad else "; every rebuilt file matched."))
        return
    m = read_manifest(w)   # a file deleted again (after a restore) keeps one row, the latest
    for path, rel, n, h, rule, check in rows: m[rel] = [rel, str(n), h, rule, check, when]
    with open(f"{w}/{MAN}.part", "w") as o: o.write(HDR); o.writelines("\t".join(x) + "\n" for x in m.values())
    os.replace(f"{w}/{MAN}.part", f"{w}/{MAN}")
    for path, *_ in rows: os.remove(path)
    s = sum(x[2] for x in rows); print(f"\nDeleted {len(rows)} files, {G(s)}; recorded in {MAN}." + (f" Kept {len(bad)} that did not rebuild identically: {', '.join(bad)}" if bad else ""))
    if rows and os.path.exists(f"{w}/README.md"):
        with open(f"{w}/README.md", "a") as o:
            o.write(f"\n**Disk cleanup ({when}):** {len(rows)} files, {G(s)}, deleted with `scripts/cleanup/clean_work.py` (docs/CLEANUP.md). "
                    f"Listed with checksums in `{MAN}`; `clean_work.py restore` remakes them.\n")

def restore_hic(c, env, log):
    d = c.p("09_hic_map"); e = env; t = e["AS_THREADS"]
    if not os.path.exists(f"{d}/targets.fasta"): hic_targets(c, f"{d}/targets.fasta")
    sh(f"cd {q(d)} && {q(e['AS_BWA'])} index targets.fasta > bwa_index.restore.log 2>&1", e)
    if not os.path.exists(f"{d}/HiC.bam"):
        r1, r2 = (os.path.join(c.w, e[k]) for k in ("AS_HIC_R1_LINK", "AS_HIC_R2_LINK"))
        sh(f"cd {q(d)} && set -o pipefail && {q(e['AS_BWA'])} mem -5SP -t {t} targets.fasta {q(r1)} {q(r2)} 2> bwa.restore.log | {q(e['AS_SAMBLASTER'])} 2> samblaster.restore.log | "
           f"{q(e['AS_SAMTOOLS'])} view - -@ 4 -S -h -b -F 3340 -o HiC.bam.part && mv HiC.bam.part HiC.bam", e)
    p = md5_stream(f"{q(e['AS_SAMTOOLS'])} view -@ 4 {q(d)}/HiC.bam | python3 {q(S)}/hic_pairs.py --minMapq 1 --maxNM 3 2>/dev/null", e)
    sh(f"{q(e['AS_SAMTOOLS'])} view -@ 2 {q(d)}/HiC.bam | python3 {q(S)}/rescue/hic_rescue_pairs.py {q(d)}/alleles08.pairs.tsv {q(d)}/rescued.restore.tsv.gz 2>/dev/null", e)
    rr = gz_md5(f"{d}/rescued.restore.tsv.gz"); os.remove(f"{d}/rescued.restore.tsv.gz")
    return f"pairs={p};rescued={rr}"

def cmd_restore(w, env, rules, a):
    m = read_manifest(w); names = {r["name"] for r in rules}; L = env.get("AS_LINEAGE", "")
    todo = [s for p, s in m.items() if s[3] in names and not os.path.exists(f"{w}/{p}")]
    if not todo: print("Nothing to restore."); return
    guard(env)
    todo.sort(key=lambda s: ORDER.index(s[3])); bad = []; doneRuns = set(); doneHic = set()
    for s in todo:
        rel = s[0]; path = f"{w}/{rel}"; r = BYNAME[s[3]]
        base = f"{w}/archive/{rel.split('/')[1]}" if rel.startswith("archive/") else w
        c = Ctx(w, base, env); f = os.path.relpath(path, base)
        if os.path.exists(path): continue
        print(f"restoring {rel}", flush=True)
        if r["kind"] == "exact":
            r["build"](c, f, path + ".part"); os.rename(path + ".part", path); ok = md5(path) == s[2]
        elif r["kind"] == "hic":
            if base in doneHic: continue
            got = restore_hic(c, env, None); doneHic.add(base)
            bam = [x for x in m.values() if x[0] == os.path.relpath(c.p("09_hic_map/HiC.bam"), w)]
            ok = not bam or got == bam[0][4]; print(f"  pairs from the remade BAM {'match' if ok else 'DIFFER from'} the recorded ones")
            for x in m.values():
                if x[3] == "hic_map" and x[2] != "-" and x[0].startswith(os.path.relpath(c.p("09_hic_map"), w)) and md5(f"{w}/{x[0]}") != x[2]: bad.append(x[0])
        else:
            run = os.path.dirname(os.path.dirname(f))
            if (base, run) in doneRuns: continue
            doneRuns.add((base, run)); t = tempfile.mkdtemp(prefix="restore_", dir=a.tmp or w)
            try:
                inp = compleasm_input(c, run, f"{t}/in.fasta")
                if not inp: print(f"  no known input for {run}: not restored"); bad.append(rel); continue
                sh(f"{q(env['AS_COMPLEASM'])} run -a {q(inp)} -o {q(t)}/out -l {q(L)} -L {q(env['AS_BUSCO_LIB'])} -t {env['AS_THREADS']} > {q(t)}/log 2>&1", env)
                for x in LINEAGE_FILES:
                    dst = c.p(run, L, x)
                    if not os.path.exists(dst): shutil.move(f"{t}/out/{L}/{x}", dst)
                same = open(f"{t}/out/summary.txt").read() == open(c.p(run, "summary.txt")).read() if os.path.exists(c.p(run, "summary.txt")) else None
            finally: shutil.rmtree(t, ignore_errors=True)
            hashes = all(md5(c.p(run, L, x)) == m.get(os.path.relpath(c.p(run, L, x), w), [None] * 3)[2] for x in LINEAGE_FILES if os.path.relpath(c.p(run, L, x), w) in m)
            # compleasm is not deterministic, so a new hash is expected and is not a failure; the summary says whether the result moved
            print(f"  {run}: checksums {'match' if hashes else 'differ (expected: compleasm is not deterministic)'}; the rerun's summary.txt "
                  + ("is the same as the kept one" if same else "DIFFERS from the kept one: the kept tables are what the stages used" if same is False else "has no kept copy to compare"))
            ok = True
        if not ok: bad.append(rel)
    print(f"\nRestored {len(todo)} files" + (f"; differing from {MAN}: {', '.join(bad)}" if bad else "; all match the record."))
    if bad: sys.exit(1)

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("cmd", choices=["list", "check", "delete", "restore"]); ap.add_argument("work")
    ap.add_argument("--rules"); ap.add_argument("--skip"); ap.add_argument("--noArchive", action="store_true"); ap.add_argument("--tmp")
    a = ap.parse_args(); w = os.path.abspath(a.work)
    if not os.path.exists(f"{w}/config.yaml"): sys.exit(f"{w} is not a work directory (no config.yaml)")
    sel = set(a.rules.split(",")) if a.rules else set(BYNAME); sel -= set(a.skip.split(",")) if a.skip else set()
    if sel - set(BYNAME): sys.exit(f"unknown rules: {sorted(sel - set(BYNAME))}; known: {', '.join(ORDER)}")
    env = stage_env(w); rules = [BYNAME[n] for n in ORDER if n in sel]
    if a.cmd != "restore" and not os.path.exists(f"{w}/state/15_validate.done"): sys.exit("The run has not finished (state/15_validate.done is missing): clean only a finished work directory")
    {"list": cmd_list, "check": cmd_check, "delete": cmd_delete, "restore": cmd_restore}[a.cmd](w, env, rules, a)

if __name__ == "__main__": main()
