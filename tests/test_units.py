#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Fast checks with no data: every script parses and states its purpose, no hard-coded paths or samples, and the edit-applying scripts behave. Run: python3 tests/test_units.py"""
import numpy as np
import ast, glob, gzip, os, re, shutil, subprocess, sys, tempfile, textwrap, unittest
R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = glob.glob(f"{R}/scripts/*.py") + glob.glob(f"{R}/scripts/*/*.py")
def run(script, *a): return subprocess.run([sys.executable, f"{R}/scripts/{script}", *a], capture_output=True, text=True)
class T(unittest.TestCase):
    def test_parse_and_doc(self):
        for f in SCRIPTS:
            t = ast.parse(open(f).read())
            if not f.endswith("chainlib.py"): self.assertTrue(ast.get_docstring(t), f)
    def test_scripts_doc_current(self):
        # docs/SCRIPTS.md is generated: a hand edit, or a script folder the generator does not read, is lost at the next run
        o = subprocess.run([sys.executable, f"{R}/tools/gen_scripts_doc.py", "--check"], capture_output=True, text=True)
        self.assertEqual(o.returncode, 0, o.stderr)
    def test_no_hard_coding(self):
        pat = re.compile(r"/home/|research/genomes|13Q010|WatsonsChoice|00_inputs|16_validate/|\bBusk\b")
        for f in SCRIPTS + glob.glob(f"{R}/stages/*.sh"):
            if f.endswith("lib.sh") or "/scripts/cleanup/" in f: continue   # these know the work directory layout, as the stages do
            for i, l in enumerate(open(f), 1):
                if f.endswith(".sh"): continue
                self.assertIsNone(pat.search(l), f"{f}:{i}: {l.strip()[:80]}")
    def test_sample_settings_required(self):
        for s, flags in (("copy_number.py", ["--hap", "--repeat", "--fai"]), ("round2.py", ["--hap", "--repeat"]), ("make_tails.py", ["--twoX", "--maxDepth", "--hap", "--ploidy"]), ("assign_arms2.py", ["--twoX", "--minAllele", "--maxPair", "--hap", "--ploidy"]), ("restore_arms.py", ["--lo", "--hi"]), ("dosage.py", ["--ploidy"])):
            r = run(s); self.assertNotEqual(r.returncode, 0, s)
            for x in flags: self.assertIn(x, r.stderr, f"{s}: {x} should be required")
    def test_overlap_longer_than_unitig(self):
        # Chain A = G[0:100], chain B = G[70:85] + G[80:200]. A's overlap with B is 30, longer than B's first unitig. Either walk direction must give G[0:200].
        import random; random.seed(1); G = "".join(random.choice("ACGT") for _ in range(200))
        rc = lambda x: x[::-1].translate(str.maketrans("ACGT", "TGCA"))
        U = {"u1": G[0:100], "u2": G[70:85], "u3": G[80:200]}
        for keep in (["u1", "u2", "u3"], ["u3", "u2", "u1"]):
            with tempfile.TemporaryDirectory() as d:
                open(f"{d}/k", "w").write("\n".join(keep) + "\n")
                open(f"{d}/f", "w").write("".join(f"{u}\t{len(x)}\n" for u, x in U.items()))
                open(f"{d}/g", "w").write("L\tu1\t+\tu2\t+\t30M\nL\tu2\t+\tu3\t+\t5M\n")
                r = run("build_chains.py", "--keep", f"{d}/k", "--gfa", f"{d}/g", "--fai", f"{d}/f", "--out", f"{d}/c"); self.assertEqual(r.returncode, 0, r.stderr)
                seq = ""
                for l in open(f"{d}/c.agp"):
                    if l[0] == "#": continue
                    x = l.split(); part = U[x[5]][int(x[6]) - 1:int(x[7])]; seq += part if x[8] == "+" else rc(part)
                self.assertIn(seq, (G, rc(G)), keep)
    def test_connectors(self):
        # a and b end at d1, a dropped unitig. c also touches d2, which leads nowhere; r is repeat class and must not count.
        with tempfile.TemporaryDirectory() as t:
            open(f"{t}/agp", "w").write("".join(f"{c}\t1\t10\t1\tW\t{u}\t1\t10\t+\n" for c, u in [("ca", "a"), ("cb", "b"), ("cc", "c")]))
            open(f"{t}/gfa", "w").write("L\ta\t+\td1\t+\t5M\nL\td1\t+\tb\t+\t5M\nL\tc\t+\td2\t+\t5M\nL\tc\t-\tr\t+\t5M\nL\tr\t+\tb\t-\t5M\n".replace("L\tc\t-\tr\t+", "L\tc\t-\tr\t+"))
            open(f"{t}/cn", "w").write("".join(f"{u}\t10\t30\t0\t0\t{k}\t1\t1\tnormal\tredundant\n" for u, k in [("d1", 60), ("d2", 60), ("r", 900)]))
            open(f"{t}/rep", "w").write("")
            r = run("restore_connectors.py", f"{t}/agp", f"{t}/gfa", f"{t}/cn", f"{t}/rep", f"{t}/out", "--repeat", "270"); self.assertEqual(r.returncode, 0, r.stderr)
            got = {l.split("\t")[0] for l in open(f"{t}/out") if l[0] != "#"}
            self.assertEqual(got, {"d1"})
    def test_mask_shared_any_ploidy(self):
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/a.agp", "w").write("".join(f"{s}\t{x}\t{y}\t{i}\tW\t{u}\t1\t{y-x+1}\t+\n" for s, x, y, i, u in [("c_h1",1,4,1,"u1"),("c_h1",5,8,2,"u2"),("c_h2",1,4,1,"u1"),("c_h2",5,8,2,"u3"),("c_h3",1,4,1,"u4"),("c_h3",5,8,2,"u2")]))
            open(f"{d}/t.fa", "w").write(">c_h1\nAAAACCCC\n>c_h2\nAAAAGGGG\n>c_h3\nTTTTCCCC\n>unplaced_x\nACGT\n")
            r = run("mask_shared.py", f"{d}/a.agp", f"{d}/t.fa", f"{d}/m"); self.assertEqual(r.returncode, 0, r.stderr)
            # keeping h1: h2's shared u1 and h3's shared u2 are masked; h1 whole
            self.assertEqual(open(f"{d}/m.h1.fasta").read(), ">c_h1\nAAAACCCC\n>c_h2\nNNNNGGGG\n>c_h3\nTTTTNNNN\n>unplaced_x\nACGT\n")
            self.assertEqual(open(f"{d}/m.h3.fasta").read(), ">c_h1\nAAAANNNN\n>c_h2\nAAAAGGGG\n>c_h3\nTTTTCCCC\n>unplaced_x\nACGT\n")
    def test_break_chain_multiple(self):
        # chain c = u1 (1-10) + u2 (11-20, 4 bp overlap trimmed) + u3 reversed (21-30, 4 bp trimmed); two breaks give c, cb, cc, each later part regaining its overlap.
        import random; random.seed(2); U = {u: "".join(random.choice("ACGT") for _ in range(14)) for u in ("u1", "u2", "u3")}
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/u.fa", "w").write("".join(f">{u}\n{s}\n" for u, s in U.items()))
            open(f"{d}/in.agp", "w").write("c\t1\t10\t1\tW\tu1\t1\t10\t+\nc\t11\t20\t2\tW\tu2\t5\t14\t+\nc\t21\t30\t3\tW\tu3\t1\t10\t-\nd\t1\t14\t1\tW\tu1\t1\t14\t+\n")
            with gzip.open(f"{d}/p.gz", "wt") as o: o.write("c\t5\tc\t15\nc\t25\td\t3\n")
            for brk, names in (("c 10\n", ["c", "cb"]), ("c 20\nc 10\n", ["c", "cb", "cc"])):
                open(f"{d}/b", "w").write(brk)
                r = run("break_chain.py", f"{d}/b", f"{d}/in.agp", f"{d}/p.gz", f"{d}/o", f"{d}/u.fa"); self.assertEqual(r.returncode, 0, r.stderr)
                rows = [l.split("\t") for l in open(f"{d}/o/chains.agp") if l[0] != "#"]
                self.assertEqual(sorted({x[0] for x in rows if x[0] != "d"}), names)
                pairs = [l.rstrip("\n").split("\t") for l in gzip.open(f"{d}/o/pairs.tsv.gz", "rt")]
            # two breaks: cb is all of u2, cc is u3 with its trimmed 4 bp back; c:15 lands on cb:9, c:25 on cc:9
            self.assertEqual([(x[0], x[5], x[6], x[7]) for x in rows if x[0] != "d"], [("c", "u1", "1", "10"), ("cb", "u2", "1", "14"), ("cc", "u3", "1", "14")])
            self.assertEqual(pairs, [["c", "5", "cb", "9"], ["cc", "9", "d", "3"]])
            # --pairsOnly lifts the same pairs from the unitig lengths alone, writing no AGP or FASTA
            open(f"{d}/u.fa.fai", "w").write("".join(f"{u}\t14\n" for u in U))
            r = run("break_chain.py", f"{d}/b", f"{d}/in.agp", f"{d}/p.gz", f"{d}/o3", f"{d}/u.fa", "--pairsOnly"); self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual([l.rstrip("\n").split("\t") for l in gzip.open(f"{d}/o3/pairs.tsv.gz", "rt")], pairs)
            self.assertEqual(sorted(os.listdir(f"{d}/o3")), ["pairs.tsv.gz"])
            open(f"{d}/b", "w").write("c 15\n")
            self.assertNotEqual(run("break_chain.py", f"{d}/b", f"{d}/in.agp", f"{d}/p.gz", f"{d}/o2", f"{d}/u.fa").returncode, 0)
    def test_weak_placed_copy_names(self):
        # collapsed order names a shared chain's copy x.g1; groups, links and pairs name the chain x. y is linked to x's copy, z is not.
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/o", "w").write("#group\trank\tchain\torientation\ng1\t1\tx.g1\t+\ng1\t2\ty\t+\ng1\t3\tz\t-\n")
            open(f"{d}/e", "w").write("#end_a\tend_b\treads\nx+\ty-\t3\n"); open(f"{d}/c", "w").write("#end_a\tend_b\tchosen\n")
            open(f"{d}/g", "w").write("".join(f"{c}\t2000000\tg1\n" for c in "xyz"))
            with gzip.open(f"{d}/p.gz", "wt") as o: o.write("z\t1\tx\t1\nz\t1\ty\t1\n")
            r = run("review/weak_placed.py", f"{d}/o", f"{d}/e", f"{d}/c", f"{d}/g", f"{d}/p.gz"); self.assertEqual(r.returncode, 0, r.stderr)
            rows = [l.split("\t") for l in r.stdout.splitlines() if l[0] != "#"]
            self.assertEqual([(x[0], x[4]) for x in rows], [("z", "2")])
            # the edited groups file names the copies too (x.g1, x.g2)
            open(f"{d}/g", "w").write("x.g1\t2000000\tg1\nx.g2\t2000000\tg2\ny\t2000000\tg1\nz\t2000000\tg1\n")
            r = run("review/weak_placed.py", f"{d}/o", f"{d}/e", f"{d}/c", f"{d}/g", f"{d}/p.gz"); self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual([l.split("\t")[0] for l in r.stdout.splitlines() if l[0] != "#"], ["z"])
    def test_agp_to_fasta(self):
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/u.fa", "w").write(">a\nAACC\n>b\nGGTT\n")
            open(f"{d}/a.agp", "w").write("s\t1\t4\t1\tW\ta\t1\t4\t+\ns\t5\t7\t2\tU\t3\tscaffold\tyes\tx\ns\t8\t11\t3\tW\tb\t1\t4\t-\n")
            r = run("agp_to_fasta.py", f"{d}/a.agp", f"{d}/u.fa", f"{d}/o.fa"); self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(open(f"{d}/o.fa").read(), ">s\nAACCNNNAACC\n")
            # a fill name defined twice with different sequence (one pass's fill overwriting another's) must stop the build
            open(f"{d}/f1.fa", "w").write(">rf_s_1\nAAAA\n"); open(f"{d}/f2.fa", "w").write(">rf_s_1\nCCCC\n")
            r = run("agp_to_fasta.py", f"{d}/a.agp", f"{d}/u.fa,{d}/f1.fa,{d}/f2.fa", f"{d}/o.fa"); self.assertNotEqual(r.returncode, 0)
            r = run("agp_to_fasta.py", f"{d}/a.agp", f"{d}/u.fa,{d}/f1.fa,{d}/f1.fa", f"{d}/o.fa"); self.assertEqual(r.returncode, 0, r.stderr)
    def test_group_edits_need_reason(self):
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/g.tsv", "w").write("#chain\tlength\tgroup\nc1\t10\tg01\nc2\t10\tg02\nc3\t10\tshared\n")
            open(f"{d}/g.shared.tsv", "w").write("#chain\tlength\tchrom\tlinks\tgroups\tscore_ratio\nc3\t10\tc0\t5\tg01,g02\t3\n")
            open(f"{d}/g.chroms.tsv", "w").write("g01\tc0\ng02\tc0\n")
            open(f"{d}/e.tsv", "w").write("c1\tunassigned\t.\n")
            self.assertNotEqual(run("apply_group_edits.py", f"{d}/g", f"{d}/e.tsv", f"{d}/o").returncode, 0)
            open(f"{d}/e.tsv", "w").write("c1\tshared\tc0\tdepth 2x\nc3\tgroup\tg02\tgraph link\n")
            r = run("apply_group_edits.py", f"{d}/g", f"{d}/e.tsv", f"{d}/o"); self.assertEqual(r.returncode, 0, r.stderr)
            g = dict(l.split("\t")[0::2] for l in open(f"{d}/o.tsv").read().splitlines() if l[0] != "#")
            self.assertEqual(g, {"c1": "shared", "c2": "g02", "c3": "g02"})
            self.assertIn("c1\t10\tc0", open(f"{d}/o.shared.tsv").read()); self.assertNotIn("c3", open(f"{d}/o.shared.tsv").read())
    def test_group_edit_with_chain(self):
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/g.tsv", "w").write("#chain\tlength\tgroup\na\t5\tg01\nb\t5\tg02\nx\t5\tunassigned\ny\t5\tunassigned\n")
            open(f"{d}/g.chroms.tsv", "w").write("g01\tc1\ng02\tc0\n")
            open(f"{d}/e", "w").write("x\tgroup\twith:b\tr\ny\tshared\twith:a\tr\n")
            r = run("apply_group_edits.py", f"{d}/g", f"{d}/e", f"{d}/o"); self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("x\t5\tg02", open(f"{d}/o.tsv").read()); self.assertIn("y\t5\tc1", open(f"{d}/o.shared.tsv").read())
            open(f"{d}/e", "w").write("x\tgroup\twith:y\tr\n")
            self.assertNotEqual(run("apply_group_edits.py", f"{d}/g", f"{d}/e", f"{d}/o").returncode, 0)
    def test_shared_in_across_translocation(self):
        # g01/g02 are chromosome c0, g03/g04 c1. Hub s is shared by g01 (c0) and g03 (c1), as across a reciprocal translocation.
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/g.tsv", "w").write("#chain\tlength\tgroup\na\t9\tg01\nb\t9\tg02\nc\t9\tg03\ns\t9\tunassigned\n")
            open(f"{d}/g.chroms.tsv", "w").write("g01\tc0\ng02\tc0\ng03\tc1\ng04\tc1\n")
            open(f"{d}/e", "w").write("s\tshared_in\twith:a,with:c\thub on both partners\n")
            r = run("apply_group_edits.py", f"{d}/g", f"{d}/e", f"{d}/f"); self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("s\t9\tc0+c1", open(f"{d}/f.shared.tsv").read()); self.assertIn("g01,g03", open(f"{d}/f.dosage.tsv").read())
            open(f"{d}/agp", "w").write("".join(f"{x}\t1\t9\t1\tW\tu{x}\t1\t9\t+\n" for x in "abcs"))
            open(f"{d}/fai", "w").write("".join(f"{x}\t9\n" for x in "abcs"))
            with gzip.open(f"{d}/p.gz", "wt") as o: o.write("s\t1\ta\t1\ns\t2\tc\t2\ns\t3\tb\t3\n")
            open(f"{d}/l", "w").write("#h\ns+\ta-\t3\tNA\t9000\ns-\tc+\t3\tNA\t9000\n")
            for x in ("st", "te"): open(f"{d}/{x}", "w").write("")
            r = run("make_shared_copies.py", f"{d}/f", f"{d}/agp", f"{d}/fai", f"{d}/p.gz", f"{d}/l", f"{d}/st", f"{d}/te", f"{d}/o"); self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("s.g01\t9\tg01\ns.g03\t9\tg03", open(f"{d}/o/groups.tsv").read())
            self.assertEqual(sorted(gzip.open(f"{d}/o/pairs.tsv.gz", "rt").read().splitlines()), ["s.g01\t1\ta\t1", "s.g03\t2\tc\t2"])   # b's group g02 has no copy
            self.assertEqual(open(f"{d}/o/end_links.tsv").read().splitlines()[1:], ["s.g01+\ta-\t3\tNA\t9000", "s.g03-\tc+\t3\tNA\t9000"])
    def test_hap_first_any_ploidy(self):
        # Two chromosomes of three haplotypes. a1 and b1 align over more bp than any true allele pair, but at 9% divergence:
        # identity weighting keeps them apart, and plain bp (--deWeight 0) does not.
        with tempfile.TemporaryDirectory() as d:
            N = [f"{c}{i}" for c in "ab" for i in (1, 2, 3)]
            open(f"{d}/fai", "w").write("".join(f"{n}\t{2000000 + 1000 * i}\t0\t0\t0\n" for i, n in enumerate(N)))
            open(f"{d}/links", "w").write("".join(f"{n}\t.\t1000\n" for n in N))
            al = [(x, y, 1500000, 0.01) for c in "ab" for x, y in ((f"{c}1", f"{c}2"), (f"{c}1", f"{c}3"), (f"{c}2", f"{c}3"))] + [("a1", "b1", 1900000, 0.09)]
            open(f"{d}/al", "w").write("".join(f"{x}\t{y}\t{n}\t{n}\t{de}\t3\n" for x, y, n, de in al))
            def chroms(*extra):
                r = run("hic_group.py", f"{d}/links", f"{d}/al", f"{d}/fai", f"{d}/o", "--chroms", "2", "--ploidy", "3", "--hapFirst", *extra); self.assertEqual(r.returncode, 0, r.stderr)
                g = {l.split("\t")[0]: l.rstrip().split("\t")[2] for l in open(f"{d}/o.tsv") if l[0] != "#"}; ch = dict(l.split() for l in open(f"{d}/o.chroms.tsv"))
                return sorted(sorted(n for n in N if ch[g[n]] == c) for c in set(ch.values())), r.stdout
            got, log = chroms(); self.assertEqual(got, [["a1", "a2", "a3"], ["b1", "b2", "b3"]], log); self.assertNotIn("WARNING", log)
            self.assertNotEqual(chroms("--deWeight", "0")[0], got)
    def test_shared_dosage(self):
        # Two chromosomes of four haplotypes. s is at 2x with Hi-C to a1 and a3; t is at 2x with graph links to a2 and a4; u is at 4x.
        with tempfile.TemporaryDirectory() as d:
            N = [f"{c}{i}" for c in "ab" for i in (1, 2, 3, 4)]; L = {n: 2800000 - 100000 * i for i, n in enumerate(N)}; L.update(s=400000, t=300000, u=200000)
            g = {n: f"g{i + 1:02d}" for i, n in enumerate(N)}                          # groups are numbered by size
            open(f"{d}/fai", "w").write("".join(f"{n}\t{x}\t0\t0\t0\n" for n, x in L.items()))
            lk = [(n, ".", 1000) for n in L] + [(x, y, 50 if x[0] == y[0] else 1) for i, x in enumerate(N) for y in N[i + 1:]]
            lk += [("s", n, 200 if n in ("a1", "a3") else 10 if n[0] == "a" else 1) for n in N] + [(c, n, 50 if n[0] == "a" else 1) for c in "tu" for n in N]
            open(f"{d}/links", "w").write("".join(f"{x}\t{y}\t{n}\n" for x, y, n in lk))
            open(f"{d}/al", "w").write("".join(f"{x}\t{y}\t1500000\t1500000\t0.01\t3\n" for i, x in enumerate(N) for y in N[i + 1:] if x[0] == y[0]))
            open(f"{d}/sh", "w").write("s\nt\nu\n")
            open(f"{d}/cls", "w").write("#chain\tlength\tmedian_depth\n" + "".join(f"{n}\t{L[n]}\t{30}\n" for n in N) + "s\t400000\t61\nt\t300000\t59\nu\t200000\t118\n")
            open(f"{d}/el", "w").write("#end_a\tend_b\treads\tq_gap_median\tgraph_overlap\nt+\ta2-\t5\t0\t6000\nt-\ta4+\t5\t0\t6000\n")
            r = run("hic_group.py", f"{d}/links", f"{d}/al", f"{d}/fai", f"{d}/g", "--chroms", "2", "--ploidy", "4", "--shared", f"{d}/sh", "--classDepth", f"{d}/cls", "--hap", "30", "--endLinks", f"{d}/el")
            self.assertEqual(r.returncode, 0, r.stderr)
            D = {l.split("\t")[0]: l.rstrip().split("\t") for l in open(f"{d}/g.dosage.tsv") if l[0] != "#"}
            self.assertEqual({k: (v[4], v[5], v[8]) for k, v in D.items()}, {"s": ("2", "g01,g03", "hic"), "t": ("2", "g02,g04", "graph")})
            # without depth, as at ploidy 2, every shared chain is in every haplotype and no dosage file is written
            r = run("hic_group.py", f"{d}/links", f"{d}/al", f"{d}/fai", f"{d}/n", "--chroms", "2", "--ploidy", "4", "--shared", f"{d}/sh"); self.assertEqual(r.returncode, 0, r.stderr)
            self.assertFalse(os.path.exists(f"{d}/n.dosage.tsv"))
            # a hand edit names the haplotypes; hic_group's call for t is carried over; a plain shared edit means every haplotype
            open(f"{d}/e", "w").write("s\tshared_in\twith:a1,g02\tgraph position\n")
            r = run("apply_group_edits.py", f"{d}/g", f"{d}/e", f"{d}/f"); self.assertEqual(r.returncode, 0, r.stderr)
            F = {l.split("\t")[0]: l.rstrip().split("\t") for l in open(f"{d}/f.dosage.tsv") if l[0] != "#"}
            self.assertEqual((F["s"][5], F["s"][8], F["t"][5]), ("g01,g02", "edit", "g02,g04"))
            open(f"{d}/e2", "w").write("s\tshared_in\tg01,g01\trepeated group\n"); self.assertNotEqual(run("apply_group_edits.py", f"{d}/g", f"{d}/e2", f"{d}/x").returncode, 0)
            open(f"{d}/e2", "w").write("s\tshared_in\tg01,g99\tunknown group\n"); self.assertNotEqual(run("apply_group_edits.py", f"{d}/g", f"{d}/e2", f"{d}/x").returncode, 0)
            open(f"{d}/e3", "w").write("t\tshared\twith:a1\tall four\n"); r = run("apply_group_edits.py", f"{d}/g", f"{d}/e3", f"{d}/y"); self.assertEqual(r.returncode, 0, r.stderr)
            self.assertNotIn("\nt\t", open(f"{d}/y.dosage.tsv").read())
            # copies: only in the named haplotypes; a pair to a haplotype without a copy is dropped
            open(f"{d}/agp", "w").write("".join(f"{n}\t1\t{x}\t1\tW\tu_{n}\t1\t{x}\t+\n" for n, x in L.items()))
            with gzip.open(f"{d}/pairs.gz", "wt") as o: o.write("s\t5\ta1\t9\ns\t5\ta3\t9\ns\t5\ts\t90\nu\t5\ta4\t9\ns\t7\tt\t8\n")
            for x in ("st", "tel"): open(f"{d}/{x}", "w").write("".join(f"{n}\t0\n" for n in L))
            os.mkdir(f"{d}/out"); r = run("make_shared_copies.py", f"{d}/f", f"{d}/agp", f"{d}/fai", f"{d}/pairs.gz", f"{d}/el", f"{d}/st", f"{d}/tel", f"{d}/out"); self.assertEqual(r.returncode, 0, r.stderr)
            cp = sorted(l.split("\t")[0] for l in open(f"{d}/out/groups.tsv") if "." in l.split("\t")[0])
            self.assertEqual(cp, ["s.g01", "s.g02", "t.g02", "t.g04", "u.g01", "u.g02", "u.g03", "u.g04"])
            self.assertEqual(gzip.open(f"{d}/out/pairs.tsv.gz", "rt").read(), "s.g01\t5\ta1\t9\ns.g01\t5\ts.g01\t90\ns.g02\t5\ts.g02\t90\nu.g04\t5\ta4\t9\ns.g02\t7\tt.g02\t8\n")
            self.assertEqual(sorted(l.split("\t")[0] for l in open(f"{d}/out/chains.agp")).count("s.g03"), 0)
            self.assertIn("t.g02+\ta2-", open(f"{d}/out/end_links.tsv").read())
    def test_reset_drops_stale_scaffold_mappings(self):
        """Read mappings to the scaffolds (stages 12, 13) are kept by a reset of their own stage, and archived by a reset of 11_build or earlier, which rebuilds their target."""
        sys.path.insert(0, R); from agentssembler import cli
        import argparse
        files = ["14_readfill/ont_vs_target.h1.paf.gz", "14_readfill/target.h1.fasta", "14_readfill/target.masked.tsv", "17_raw_reads/raw_vs_target.h1.paf.gz", "17_raw_reads/raw.h1.paf.gz", "06_chains/ont_vs_chains.paf.gz"]
        for stage, remap, gone in (("12_fill", False, {"14_readfill/target.h1.fasta", "14_readfill/target.masked.tsv"}),
                                   ("09_order", False, set(files[:5])), ("12_fill", True, set(files[:5]))):
            with tempfile.TemporaryDirectory() as w:
                for d in ("state", "decisions", "reports"): os.makedirs(f"{w}/{d}")
                open(f"{w}/config.yaml", "w").write("sample: s\n")
                for f in files: os.makedirs(os.path.dirname(f"{w}/{f}"), exist_ok=True); open(f"{w}/{f}", "w").close()
                cli.cmd_reset(argparse.Namespace(work=w, stage=stage, remap=remap, keepDecision=False))
                self.assertEqual({f for f in files if not os.path.exists(f"{w}/{f}")}, gone, (stage, remap))
    def test_read_formats(self):
        """init accepts FASTA or FASTQ reads, gzipped or not, and raw BAM/CRAM/SAM; raw_bams is read as well as raw_reads."""
        sys.path.insert(0, R); from agentssembler import cli
        d = tempfile.mkdtemp()
        try:
            fa, fq = ">r1\nACGT\n", "@r1\nACGT\n+\nIIII\n"
            for n, t in (("a.fa", fa), ("a.fq", fq), ("a.txt", "not reads\n")): open(f"{d}/{n}", "w").write(t)
            for n, t in (("a.fa.gz", fa), ("a.fq.gz", fq)):
                with gzip.open(f"{d}/{n}", "wt") as o: o.write(t)
            want = {"a.fa": "fasta", "a.fa.gz": "fasta", "a.fq": "fastq", "a.fq.gz": "fastq", "a.txt": None, "x.bam": "bam", "x.cram": "bam"}
            for n, k in want.items(): self.assertEqual(cli.read_format(f"{d}/{n}"), k, n)
            self.assertEqual(cli.raw_reads({"inputs": {"raw_reads": ["a"], "raw_bams": ["b"]}}), ["a", "b"])
            self.assertEqual(cli.raw_reads({"inputs": {}}), [])
            # 00_inputs links take the extension of what they hold; unknown or missing files, and builds linked under the old names, keep the old names
            w = f"{d}/work"; os.makedirs(f"{w}/00_inputs")
            C = lambda p: cli.input_link(w, p, "ont_corrected_reads", dict(fasta=".fasta", fastq=".fastq"), "ont_corrected_reads.fasta")
            Hc = lambda p: cli.input_link(w, p, "hic_R1", dict(fasta=".fa", fastq=".fq"), "hic_R1.fq.gz")
            want = {"a.fa": ".fasta", "a.fa.gz": ".fasta.gz", "a.fq": ".fastq", "a.fq.gz": ".fastq.gz", "a.txt": ".fasta", "gone.fq": ".fasta"}
            for n, x in want.items(): self.assertEqual(C(f"{d}/{n}"), "00_inputs/ont_corrected_reads" + x, n)
            for n, x in {"a.fq.gz": ".fq.gz", "a.fq": ".fq", "a.fa": ".fa"}.items(): self.assertEqual(Hc(f"{d}/{n}"), "00_inputs/hic_R1" + x, n)
            os.symlink(f"{d}/a.fq.gz", f"{w}/00_inputs/ont_corrected_reads.fasta"); self.assertEqual(C(f"{d}/a.fq.gz"), "00_inputs/ont_corrected_reads.fasta")
            os.symlink(f"{d}/a.fq.gz", f"{w}/00_inputs/ont_corrected_reads.fastq.gz"); self.assertEqual(C(f"{d}/a.fq.gz"), "00_inputs/ont_corrected_reads.fastq.gz")
        finally: shutil.rmtree(d)
    def test_length_scale(self):
        # Chromosome-scale windows follow mean chromosome length (chain bp over chromosomes x ploidy), anchored at 36 Mb = 1.
        sys.path.insert(0, R); from agentssembler import cli
        def env(chroms, ploidy, chains, extra=""):
            d = tempfile.mkdtemp(); os.makedirs(f"{d}/state"); os.makedirs(f"{d}/08_final_chains"); os.makedirs(f"{d}/decisions")
            with open(f"{d}/config.yaml", "w") as o: o.write(f"sample: x\nchromosomes: {chroms}\nploidy: {ploidy}\ninputs: {{unitigs_fasta: a, unitigs_gfa: b, corrected_reads: c, hic_r1: d, hic_r2: e}}\n{extra}")
            with open(f"{d}/08_final_chains/chains.fasta.fai", "w") as o: o.write("".join(f"c{i}\t{n}\t0\t0\t0\n" for i, n in enumerate(chains)))
            return d
        d = env(9, 2, [9000000] * 72); self.assertNotIn("AS_SEED_LEN", cli.stage_env(d))              # stage 04 not finished: nothing to measure yet
        open(f"{d}/state/04_chains.done", "w").close(); e = cli.stage_env(d)
        self.assertEqual([e[k] for k in ("AS_SCALE", "AS_SEED_LEN", "AS_BREAK_DIST", "AS_BREAK_FLANK", "AS_BREAK_MIN_SIDE", "AS_MAP_BIN", "AS_DOSAGE_BIN")], ["1.0", "1000000", "500000", "1000000", "100000", "250000", "500000"])
        with open(f"{d}/08_final_chains/chains.fasta.fai", "a") as o: o.write("more\t500000000\t0\t0\t0\n")
        self.assertEqual(cli.stage_env(d)["AS_SEED_LEN"], "1000000")                                        # kept for the whole build
        for chroms, ploidy, chains, extra, seed, why in ((16, 1, [500000] * 16, "", "20000", "floor"), (23, 2, [30000000] * 207, "", "3750000", "mean chromosome length"),
                                                        (21, 1, [2000000] * 7000, "", "2000000", "chain N50"), (16, 1, [750000] * 16, "settings: {length_scale: 1}\n", "1000000", None)):
            d = env(chroms, ploidy, chains, extra); open(f"{d}/state/04_chains.done", "w").close()
            self.assertEqual(cli.stage_env(d)["AS_SEED_LEN"], seed)
            if why: self.assertEqual(cli.length_scale(d)["limited_by"], why)
    def test_hic_density(self):
        # min_within is 0.13 x H, H the length-weighted median within-chain density of chains at least the seed length long.
        sys.path.insert(0, R); from agentssembler import cli
        def env(stats, extra="", done=("04_chains", "06_evidence")):
            d = tempfile.mkdtemp(); [os.makedirs(f"{d}/{x}") for x in ("state", "08_final_chains", "10_groups", "decisions")]
            with open(f"{d}/config.yaml", "w") as o: o.write(f"sample: x\nchromosomes: 1\nploidy: 2\ninputs: {{unitigs_fasta: a, unitigs_gfa: b, corrected_reads: c, hic_r1: d, hic_r2: e}}\n{extra}")
            with open(f"{d}/08_final_chains/chains.fasta.fai", "w") as o: o.write("".join(f"c{i}\t{n}\t0\t0\t0\n" for i, (n, _) in enumerate(stats)))
            with open(f"{d}/10_groups/chain_stats.tsv", "w") as o: o.write("#chain\tlength\twithin_per_mb\n" + "".join(f"c{i}\t{n}\t{x}\n" for i, (n, x) in enumerate(stats)))
            for s in done: open(f"{d}/state/{s}.done", "w").close()
            return d
        # 72 Mb of chains over 2 haplotypes: scale 1, seeds of 1 Mb. Short dense chains and a long sparse one must not move H off 5000.
        S = [(10000000, 5000)] * 3 + [(20000000, 6000)] + [(12000000, 4000)] + [(5000000, 100)] + [(500000, 90000)] * 6
        e = cli.stage_env(env(S)); self.assertEqual(e["AS_MIN_WITHIN"], "650")
        self.assertNotIn("AS_MIN_WITHIN", cli.stage_env(env(S, done=("04_chains",))))                    # stage 06 not finished: nothing to measure yet
        self.assertEqual(cli.stage_env(env(S, "settings: {min_within: 700}\n"))["AS_MIN_WITHIN"], "700")  # a number pins it
        d = env(S); cli.stage_env(d)
        with open(f"{d}/10_groups/chain_stats.tsv", "a") as o: o.write("big\t900000000\t99999\n")
        self.assertEqual(cli.stage_env(d)["AS_MIN_WITHIN"], "650")                                         # kept for the whole build
        self.assertEqual(cli.stage_env(env(S, done=("04_chains", "05_hic_chains", "06_evidence", "07_group")))["AS_MIN_WITHIN"], "700")   # a build that grouped with 700 keeps it
    def test_join_reversal(self):
        with tempfile.TemporaryDirectory() as d:
            # pre: a(1-100) gap b(1-100,-) gap c(1-100). cand closes both gaps with trims. reverse the first only: a and b's facing ends return; b's other trim stays
            pre = [("s",1,100,1,"W","a",1,100,"+"),("s",101,110,2,"U",10,"scaffold","yes","p"),("s",111,210,3,"W","b",1,100,"-"),("s",211,220,4,"U",10,"scaffold","yes","p"),("s",221,320,5,"W","c",1,100,"+")]
            cand = [("s",1,95,1,"W","a",1,95,"+"),("s",96,175,2,"W","b",11,90,"-"),("s",176,275,3,"W","c",1,100,"+")]
            w = lambda p, R: open(p, "w").write("".join("\t".join(map(str, r)) + "\n" for r in R))
            w(f"{d}/pre.agp", pre); w(f"{d}/cand.agp", cand); open(f"{d}/dec.tsv", "w").write("s\ta\tb\treopen, array\n")
            r = run("apply_join_decisions.py", f"{d}/cand.agp", f"{d}/pre.agp", f"{d}/dec.tsv", f"{d}/o.agp"); self.assertEqual(r.returncode, 0, r.stderr)
            o = [l.split("\t") for l in open(f"{d}/o.agp").read().splitlines() if l[0] != "#"]
            self.assertEqual([x[5] for x in o], ["a", "10", "b", "c"])
            self.assertEqual(o[0][6:8], ["1", "100"]); self.assertEqual(o[2][6:8], ["11", "100"]); self.assertEqual(o[3][6:8], ["1", "100"])
            self.assertEqual(o[-1][2], str(100 + 10 + 90 + 100))
            # the same scaffold written the other way round (edits/flip.txt) after pre.agp: the join is still found and the same unitig ends restored
            w(f"{d}/cand.agp", [("s",1,100,1,"W","c",1,100,"-"),("s",101,180,2,"W","b",11,90,"+"),("s",181,275,3,"W","a",1,95,"-")])
            r = run("apply_join_decisions.py", f"{d}/cand.agp", f"{d}/pre.agp", f"{d}/dec.tsv", f"{d}/o.agp"); self.assertEqual(r.returncode, 0, r.stderr)
            o = [l.split("\t") for l in open(f"{d}/o.agp").read().splitlines() if l[0] != "#"]
            self.assertEqual([x[5] for x in o], ["c", "b", "10", "a"])
            self.assertEqual(o[1][6:8], ["11", "100"]); self.assertEqual(o[3][6:8], ["1", "100"])
    def test_homolog_and_tail_alarms(self):
        sys.path.insert(0, R)
        from agentssembler import report as r
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/f.fai", "w").write("c1_h1\t30000000\t0\t60\t61\nc1_h2\t25000000\t0\t60\t61\nc2_h1\t30000000\t0\t60\t61\nc2_h2\t29500000\t0\t60\t61\n")
            open(f"{d}/a.agp", "w").write("c1_h2\t1\t100\t1\tW\tu\t1\t100\t+\nc1_h2\t101\t200\t2\tN\t100\tscaffold\tyes\tna\n")
            A = r.homolog_lengths(f"{d}/f.fai", f"{d}/a.agp")
            self.assertEqual(len(A), 1); self.assertIn("c1_h2", A[0]); self.assertIn("1 gap:", A[0])
            open(f"{d}/t.tsv", "w").write("#h\nc1_h2\t9\tright\t14\tunplaced_a\t14\t1\t1\t2\t+\t60\nc1_h2\t9\tleft\t10\tunplaced_b\t2\t1\t1\t2\t+\t60\nc1_h2\t9\tleft\t4\tc1_h2\t4\t1\t1\t2\t+\t60\nc1_h2\t7\tleft\t4\tc1_h1\t4\t1\t1\t2\t+\t60\nc2_h1\t5\tleft\t3\t.\t0\t0\t.\t.\t.\t.\n"
                                       "c2_h2\t8\tright\t4\tunplaced_c\t4\t1\t1\t2\t+\t0\t4\n")
            T = r.tail_alarms(f"{d}/t.tsv")
            self.assertEqual(len(T), 3); self.assertIn("14 of 14", T[0]); self.assertIn("duplicates placed", T[1]); self.assertIn("c2_h2:8 right -> c", T[1]); self.assertIn("1 gap edges", T[2])
    def test_rescue_diploid_only(self):
        # the Hi-C allele rescue is refused above ploidy 2 (13Q010 trial), and the stage 05 report does not suggest it there
        sys.path.insert(0, R)
        from agentssembler import cli, report as r
        self.assertIsNone(cli.rescue_refused(2, True)); self.assertIsNone(cli.rescue_refused(3, False)); self.assertIsNone(cli.rescue_refused(4, "false"))
        self.assertIn("diploids only", cli.rescue_refused(3, True)); self.assertIn("diploids only", cli.rescue_refused(4, "true"))
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(f"{d}/09_hic_map")
            open(f"{d}/09_hic_map/rescue_signal.txt", "w").write("Hi-C-blind: 41 chains\nALARM: 8.4 Mb of allele-paired chains have almost no unique Hi-C.\n")
            o, a = r.r_hic(d, {"ploidy": 2}, {}); self.assertEqual(len(a), 1)
            o, a = r.r_hic(d, {"ploidy": 3}, {}); self.assertEqual(a, []); self.assertTrue(any("diploids only" in x for x in o))
    def test_cleanup_rules(self):
        # clean_work.py finds exactly the rebuildable files (archive-only ones only in archive/), never a regression seed or compared file
        sys.path.insert(0, f"{R}/scripts/cleanup"); import clean_work as cw, yaml
        L = "lin_odb12"; env = {"AS_SAMPLE": "S", "AS_LINEAGE": L, "AS_COMPLEASM": sys.executable, "AS_HIC_R1": __file__, "AS_HIC_R2": __file__}
        go = ["09_hic_map/HiC.bam", "09_hic_map/targets.fasta"] + [f"09_hic_map/targets.fasta.{x}" for x in ("amb", "ann", "bwt", "pac", "sa")] + \
             [f"{d}/{L}/{f}" for d in ("01_genes", "16_validate/compleasm/h1", "16_validate/compleasm/h1_unloc", "16_validate/compleasm/unplaced") for f in cw.LINEAGE_FILES] + \
             ["12_shared/pre/target.fasta", "17_raw_reads/fill/in.fasta", "17_raw_reads/fill/h1/in.fasta", "17_raw_reads/fill/candidate.fasta", "13_scaffold/bridges/reads.fasta"]
        archOnly = ["13_scaffold/target.fasta", "14_readfill/h1/in.fasta", "14_readfill/target.h1.fasta"]
        keep = ["09_hic_map/pairs.tsv.gz", "09_hic_map/pairs.rescued.tsv.gz", "09_hic_map/decoys.txt", "09_hic_map/chains.fasta", "09_hic_map/alleles08.pairs.tsv",
                "00_inputs/unitigs.fasta", "12_shared/pre/S.fasta", "12_shared/pre/S.unplaced.fasta", "13_scaffold/S.graph.fasta", "13_scaffold/S.graph.unplaced.fasta", "13_scaffold/S.graph.agp",
                "15_final/S.final.fasta", "15_final/S.final.unplaced.fasta", "15_final_raw/S.final.agp", "15_final_raw/S.final.fasta", "15_final_raw/S.final.unplaced.fasta", "16_validate/unloc/S.unloc.fasta", "17_raw_reads/fill/candidate.agp", "17_raw_reads/fill/fills.fasta",
                "14_readfill/fills.fasta", "17_raw_reads/fill/h1/names.txt", "14_readfill/h1/names.txt", "13_scaffold/bridges/ids.txt", "16_validate/compleasm/tmp/hap1.fasta",
                "16_validate/compleasm/odd/" + L + "/miniprot_output.gff", "16_validate/compleasm/odd/" + L + "/full_table_busco_format.tsv"] + \
               [f"{d}/{L}/full_table_busco_format.tsv" for d in ("01_genes", "16_validate/compleasm/h1", "16_validate/compleasm/h1_unloc", "16_validate/compleasm/unplaced")]
        for sp in glob.glob(f"{R}/tests/specs/*.yaml"):
            with open(sp) as fh: y = yaml.safe_load(fh)
            keep += list((y.get("seed") or {})) + [r[0] for r in y.get("compare", [])]
        keep = [k.replace("13Q010", "S").replace("WatsonsChoice", "S") for k in keep]
        with tempfile.TemporaryDirectory() as w:
            arch = f"{w}/archive/20260101_000000"
            for b in (w, arch):
                for f in set(go + archOnly + keep):
                    os.makedirs(os.path.dirname(f"{b}/{f}"), exist_ok=True)
                    with open(f"{b}/{f}", "w") as fh: fh.write(">x\nA\n")
            env["AS_READS"] = f"{w}/00_inputs/unitigs.fasta"
            found, kept = cw.candidates(w, env, cw.RULES, False)
            got = sorted(os.path.relpath(c.p(f), w) for c, r, f, n in found)
            want = sorted(go + [f"archive/20260101_000000/{f}" for f in go + archOnly])
            self.maxDiff = None; self.assertEqual(sorted(set(want) - set(got)), []); self.assertEqual(sorted(set(got) - set(want)), [])
            for k in keep: self.assertNotIn(k, got)
            os.remove(f"{w}/09_hic_map/decoys.txt")   # a missing source keeps what it would remake
            found, kept = cw.candidates(w, env, cw.RULES, True)
            self.assertNotIn("09_hic_map/HiC.bam", [f for c, r, f, n in found]); self.assertIn("09_hic_map/HiC.bam", [p for p, why in kept])
    def test_hic_report_hub_alarm(self):
        # stage 05 report with one breakable and one inside hub core, and no other stage 05 files
        sys.path.insert(0, R)
        from agentssembler import report as r
        with tempfile.TemporaryDirectory() as d:
            os.makedirs(f"{d}/10_groups")
            open(f"{d}/10_groups/hub_cores.d04.tsv", "w").write("#chain\n"
                "cA\t500000\tspecific\t\t100000\t300000\t200\t150\t150\t0.2\t120000\t310000\tunitig/unitig\t0.90\tyes\n"
                "cB\t500000\tspecific\t\t100000\t300000\t200\t150\t150\t0.2\t.\t.\tinside/inside\t1.00\tno\n")
            o, a = r.r_hic(d, None, {})
            self.assertEqual(len(a), 1); self.assertIn("cA at 120000/310000", a[0]); self.assertTrue(any("not breakable here: 1" in x for x in o))
    def test_unplaced_components(self):
        # u1 placed on s1; u2 - u3 unplaced, joined to each other and to u1 in the graph; u4 unplaced, alone
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/p.tsv", "w").write("#h\nc1\t1000000\t.\ts1\tplaced\nc2\t400000\t.\t.\tunplaced\nc3\t400000\t.\t.\tunplaced\nc4\t500000\t.\t.\tunplaced\n")
            open(f"{d}/c.agp", "w").write("".join(f"c{i}\t1\t10\t1\tW\tu{i}\t1\t10\t+\n" for i in range(1, 5)))
            open(f"{d}/g.gfa", "w").write("L\tu1\t+\tu2\t+\t0M\nL\tu2\t+\tu3\t+\t0M\n")
            open(f"{d}/l.tsv", "w").write("")
            o = run("review/unplaced_components.py", f"{d}/p.tsv", f"{d}/c.agp", f"{d}/g.gfa", f"{d}/l.tsv")
            self.assertEqual(o.returncode, 0, o.stderr)
            self.assertRegex(o.stdout, r"2 chains, 0\.80 Mb; graph contacts: \{'s1': 1\}")
            self.assertRegex(o.stdout, r"1 chains, 0\.50 Mb; graph contacts: none")
    @unittest.skipUnless(os.environ.get("MINIMAP2") or shutil.which("minimap2"), "needs minimap2")
    def test_gap_tails(self):
        # A read aligned to the right flank's start, whose unaligned part lies in an unplaced sequence: its tail must be traced there.
        import random; random.seed(2); seq = lambda n: "".join(random.choice("ACGT") for _ in range(n))
        L, Rt, U = seq(20000), seq(20000), seq(20000); read = U[5000:12000] + Rt[:6000]
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/a.agp", "w").write("s\t1\t20000\t1\tW\tu1\t1\t20000\t+\ns\t20001\t20100\t2\tN\t100\tscaffold\tyes\tna\ns\t20101\t40100\t3\tW\tu2\t1\t20000\t+\n")
            open(f"{d}/all.fa", "w").write(f">s\n{L}{'N' * 100}{Rt}\n>un\n{U}\n")
            # the reads as FASTA and FASTQ, gzipped or not; a quality line starting with '>' once passed for a FASTA header
            reads = {"r.fa": f">r1\n{read}\n", "r.fq": f"@r1\n{read}\n+\n{'>' * len(read)}\n"}
            for n, t in reads.items():
                open(f"{d}/{n}", "w").write(t)
                with gzip.open(f"{d}/{n}.gz", "wt") as o: o.write(t)
            with gzip.open(f"{d}/m.paf.gz", "wt") as o: o.write("r1\t13000\t7000\t13000\t+\ts\t40100\t20100\t26100\t6000\t6000\t60\n")
            mm, sk = os.environ.get("MINIMAP2") or shutil.which("minimap2"), os.environ.get("SEQKIT") or shutil.which("seqkit")
            if not (mm and sk): self.skipTest("needs minimap2 and seqkit: on PATH, or MINIMAP2= and SEQKIT=")
            for n in ("r.fa", "r.fa.gz", "r.fq", "r.fq.gz"):
                x = run("review/gap_tails.py", f"{d}/a.agp", f"{d}/{n}", f"{d}/all.fa", f"{d}/o.tsv", "--paf", f"{d}/m.paf.gz", "--minimap2", mm, "--seqkit", sk, "--threads", "2")
                self.assertEqual(x.returncode, 0, x.stderr)
                F = [l.split("\t") for l in open(f"{d}/o.tsv") if l[0] != "#"]
                self.assertEqual(len(F), 1, n); self.assertEqual(F[0][:6], ["s", "20000", "right", "1", "un", "1"], n)
    def test_complement_genes(self):
        # c1_h1 lacks gene b, carried by c2_h2: numbering decides whether a set is complete. Genes a and d have all copies; gene e has one copy short.
        rows = {"h1": [("a", "c1_h1"), ("d", "c2_h1"), ("e", "c1_h1")], "h2": [("a", "c1_h2"), ("b", "c1_h2"), ("b", "c2_h2"), ("d", "c2_h2")]}
        with tempfile.TemporaryDirectory() as d:
            for h, v in rows.items():
                with open(f"{d}/{h}.tsv", "w") as o:
                    for g, s in v: o.write(f"{g}\tComplete\t{s}\t1\t2\t+\t1\t1\n")
            x = run("validate/complement_genes.py", "2", f"{d}/h1.tsv", f"{d}/h2.tsv")
            self.assertEqual(x.returncode, 0, x.stderr)
            self.assertIn("not complete in all sets 2.", x.stdout.splitlines()[0])
            self.assertIn("not complete in all sets 1.", x.stdout.splitlines()[1])
            self.assertIn("Gene copies short of 2 per gene, on any scaffold: 1 of 8", x.stdout)
    def test_chain_forks_lost_route(self):
        # Chain X = a, h; h is a 2x hub whose left side also links to dropped d, which reaches the end of chain Y. Both chains are in no group:
        # that must not read as the same group. With placement (X placed, Y unplaced at 1x) the row is LOST, with the break just before h.
        with tempfile.TemporaryDirectory() as d:
            p = lambda n, s: open(f"{d}/{n}", "w").write(textwrap.dedent(s).lstrip())
            p("c.agp", "X\t1\t1000\t1\tW\ta\t1\t1000\t+\nX\t1001\t2000\t2\tW\th\t1\t1000\t+\nY\t1\t100000\t1\tW\tt\t1\t100000\t+\n")
            p("g.gfa", "L\ta\t+\th\t+\t0M\nL\td\t+\th\t+\t0M\nL\tt\t+\td\t+\t0M\n")
            p("cn.tsv", "#u\tlen\tprim\tuq\tsec\tclass\nh\t1000\t30\t0\t0\t60\na\t1000\t30\t0\t0\t30\nd\t1000\t30\t0\t0\t30\nt\t100000\t30\t0\t0\t30\n")
            p("gr.tsv", "#chain\tlength\tgroup\nX\t2000\tunassigned\nY\t100000\tunassigned\n")
            p("ch.tsv", "g01\tc1\n")
            p("pl.tsv", "#chain\tlength\tgroup\tname\tstatus\nX\t2000\tunassigned\t.\tplaced:bridge\nY\t100000\tunassigned\t.\tunplaced:no_group\n")
            p("cls.tsv", "#chain\tlength\tmedian_depth\nX\t2000\t45\nY\t100000\t30\n")
            base = [f"{d}/c.agp", f"{d}/g.gfa", f"{d}/cn.tsv", f"{d}/gr.tsv", f"{d}/ch.tsv", "--hap", "30"]
            r = run("review/chain_forks.py", *base); self.assertEqual(r.returncode, 0, r.stderr)
            row = r.stdout.splitlines()[1].split("\t")
            self.assertEqual((row[0], row[1], row[9], row[11], row[12]), (".", "X", "Y", "ungrouped", "1000"))
            r = run("review/chain_forks.py", *base, "--placement", f"{d}/pl.tsv", "--chainClass", f"{d}/cls.tsv", "--oneLo", "19", "--oneHi", "46")
            self.assertEqual(r.returncode, 0, r.stderr); self.assertEqual(r.stdout.splitlines()[1].split("\t")[0], "LOST")
            p("tails.tsv", "#end\tgroup\tpath\tlabel\tpairs\tbp\nZ+\tg01\tY-\tonly\t1\t100000\n")
            r = run("review/chain_forks.py", *base, "--placement", f"{d}/pl.tsv", "--chainClass", f"{d}/cls.tsv", "--oneLo", "19", "--oneHi", "46", "--tails", f"{d}/tails.tsv")
            self.assertEqual(r.stdout.splitlines()[1].split("\t")[0], ".", "a chain on a tail is in use, not lost")
            # with the chain's own depth over the hub: under 2x (1.54 x 30) is an allele bubble, BUBBLE; 2x keeps LOST
            for depth, flag in (("31", "BUBBLE"), ("60", "LOST")):
                p("w.tsv", f"#unitig\tstart\tend\ttotal\tunique\tsecondary\nX\t1000\t2000\t{depth}\t{depth}\t0\n")
                r = run("review/chain_forks.py", *base, "--placement", f"{d}/pl.tsv", "--chainClass", f"{d}/cls.tsv", "--oneLo", "19", "--oneHi", "46", "--windows", f"{d}/w.tsv")
                self.assertEqual(r.returncode, 0, r.stderr); self.assertEqual(r.stdout.splitlines()[1].split("\t")[0], flag, depth)
    def test_gene_evidence(self):
        # genes lifted from unitigs to chains: a short gene with a copy on an unassigned chain, made whole by a rescued chain; tail genes
        with tempfile.TemporaryDirectory() as d:
            p = lambda n, x: open(f"{d}/{n}", "w").write(x)
            p("agp", "cA\t1\t1000\t1\tW\tu1\t1\t1000\t+\ncB\t1\t1000\t1\tW\tu2\t1\t1000\t-\ncC\t1\t500\t1\tW\tu3\t1\t500\t+\n")
            p("genes", "# Busco id\tStatus\tSequence\tGene Start\tGene End\n"
                       "G1\tDuplicated\tu1\t100\t200\nG1\tDuplicated\tu3\t100\t200\nG2\tDuplicated\tu1\t300\t400\nG2\tDuplicated\tu2\t300\t400\n"
                       "G3\tFragmented\tu3\t1\t50\nG4\tMissing\t\t\t\n")
            p("groups", "#chain\tlength\tgroup\ncA\t1000\tg01\ncB\t1000\tg02\ncC\t500\tunassigned\n")
            p("resc", "#chain\tlength\tchrom\tgroup\tcall\trescued_links\tratio_to_next\tround\ncC\t500\tc1\tg02\tallele\t12\t3.00\t1\n")
            r = run("review/genes_placed.py", "--genes", f"{d}/genes", "--agp", f"{d}/agp", "--groups", f"{d}/groups", "--ploidy", "2", "--out", f"{d}/out", "--rescued", f"{d}/resc")
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("short of 2 copies in the groups: 1; of these, with a copy on an unassigned chain: 1", r.stdout)
            self.assertIn("short genes made whole: 1 of 1", r.stdout)
            p("pl", "#chain\tlength\tgroup\tname\tstatus\ncA\t1000\tg01\tchr1_h1\tplaced\ncB\t1000\tg02\tchr1_h2\tplaced\ncC\t500\t.\t.\tunplaced:no_group\n")
            p("tails", "#end\tgroup\tpath\tlabel\tpairs\tbp\ncA+\tg01\tcC+\tonly\t3\t500\n")
            r = run("review/tail_genes.py", "--genes", f"{d}/genes", "--agp", f"{d}/agp", "--placement", f"{d}/pl", "--tails", f"{d}/tails", "--ploidy", "2", "--out", f"{d}/tg")
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(open(f"{d}/tg").read().splitlines()[1].split("\t")[4:], ["1", "1"], "tail cC carries G1, short of a placed copy")
            # unplaced chains off every tail: cC with no tail is listed, with the LOST fork that leads to it; on a tail it is not
            p("forks", "#flag\tchain\tgroup\tstart\tend\tunitig\tprimary\tclass\tvia\tto_chain\nLOST\tcA\tg01\t1\t10\tu1\t1\t2\tx\tcC\n")
            r = run("review/tail_genes.py", "--genes", f"{d}/genes", "--agp", f"{d}/agp", "--placement", f"{d}/pl", "--ploidy", "2", "--unplacedOut", f"{d}/ug", "--forks", f"{d}/forks")
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(open(f"{d}/ug").read().splitlines()[1:], ["cC\t500\tunplaced:no_group\t1\t1\tcA"])
            r = run("review/tail_genes.py", "--genes", f"{d}/genes", "--agp", f"{d}/agp", "--placement", f"{d}/pl", "--ploidy", "2", "--unplacedOut", f"{d}/ug", "--tails", f"{d}/tails")
            self.assertEqual(open(f"{d}/ug").read().splitlines()[1:], [])
    def test_apply_order_edits_unplace(self):
        # an unplaced chain leaves the order and becomes unassigned in the groups; blank lines are skipped
        with tempfile.TemporaryDirectory() as d:
            p = lambda n, x: open(f"{d}/{n}", "w").write(x)
            p("e", "#chain\taction\tgroup\tx\torientation\treason\n\ny\tunplace\t.\t.\t.\tsatellite, no unique pairs\n")
            p("o", "#group\trank\tchain\to\tlen\tlinks\tollr\tpllr\tend_link\ng1\t1\tx\t+\t10\t5\t1\t1\t.\ng1\t2\ty\t+\t10\t5\t1\t1\t.\ng1\t3\tz\t+\t10\t5\t1\t1\t.\n")
            p("g", "#chain\tlength\tgroup\nx\t10\tg1\ny\t10\tg1\nz\t10\tg1\n")
            r = run("apply_order_edits.py", f"{d}/e", f"{d}/o", f"{d}/g", f"{d}/oo", f"{d}/go"); self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual([l.split("\t")[2] for l in open(f"{d}/oo") if l[0] != "#"], ["x", "z"])
            self.assertIn("y\t10\tunassigned", open(f"{d}/go").read())
    def test_choose_bridges_none(self):
        # no bridges: no read extraction or mapping (the tools do not exist here), a header-only choices table, no blank id line
        with tempfile.TemporaryDirectory() as d:
            for n in ("b", "e"): open(f"{d}/{n}", "w").write("#header\n")
            open(f"{d}/c.fa", "w").write(">c1\nACGT\n")
            with gzip.open(f"{d}/p.gz", "wt") as o: o.write("")
            r = run("choose_bridges.py", f"{d}/b", f"{d}/e", f"{d}/c.fa", f"{d}/p.gz", f"{d}/missing.fa", f"{d}/o", "--minimap2", "/nonexistent", "--seqkit", "/nonexistent")
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(open(f"{d}/o/choices.tsv").read().count("\n"), 1); self.assertEqual(open(f"{d}/o/ids.txt").read(), "")
    def test_phase_scan_calibrated(self):
        # 20 bins on each of two homologs: a background of own far pairs and a little homolog crosstalk; bin 5 of a_h1 sends many pairs to
        # a_h2 (a phase switch). Bin 15 has lost its own far pairs but has no more homolog pairs than the rest: not flagged (the share test).
        import random; random.seed(4); B = 100000; rows = []
        for h, o in (("a_h1", "a_h2"), ("a_h2", "a_h1")):
            for b in range(20):
                if not (h == "a_h1" and b == 15):
                    rows += [f"{h}\t{b * B + 10}\t{h}\t{((b + 10) % 20) * B + 10}\n"] * 100   # own pairs 1 Mb apart
                rows += [f"{h}\t{b * B + 20}\t{h}\t{b * B + 30}\n"] * 200
                rows += [f"{h}\t{b * B + 40}\t{o}\t{b * B + 40}\n"] * (300 if (h == "a_h1" and b == 5) else 5)
        with tempfile.TemporaryDirectory() as d:
            with gzip.open(f"{d}/p.gz", "wt") as fh: fh.writelines(rows)
            open(f"{d}/agp", "w").write("a_h1\t1\t2000000\t1\tW\tu1\t1\t2000000\t+\na_h2\t1\t2000000\t1\tW\tu2\t1\t2000000\t+\n")
            r = run("validate/phase_scan.py", f"{d}/p.gz", f"{d}/agp", "--far", "500000"); self.assertEqual(r.returncode, 0, r.stderr)
            runs = [l.split() for l in r.stdout.splitlines()[2:]]
            self.assertIn(["a_h1", "0.5", "0.6"], [x[:3] for x in runs]); self.assertNotIn("1.5", [x[1] for x in runs if x[0] == "a_h1"])
    def test_rescue_signal(self):
        # chain X: allele of Y, no unique pairs, many ambiguous ends: Hi-C-blind; Y has unique pairs
        with tempfile.TemporaryDirectory() as d:
            p = lambda n, x: open(f"{d}/{n}", "w").write(x)
            p("fai", "X\t100000\nY\t300000\nZ\t2000000\n"); p("al", "#a\tb\nX\tY\t50000\t50000\n")
            with gzip.open(f"{d}/u.gz", "wt") as o: o.write("".join(f"Y\t{i}\tZ\t{i}\n" for i in range(50)))
            with gzip.open(f"{d}/r.gz", "wt") as o: o.write("".join(f"X\t{i}\tZ\t{i}\tau\n" for i in range(30)))
            r = run("rescue/rescue_signal.py", f"{d}/u.gz", f"{d}/r.gz", f"{d}/fai", f"{d}/al", f"{d}/out")
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("Hi-C-blind", r.stdout); self.assertIn("ALARM", r.stdout)
            self.assertEqual([l.split("\t")[0] for l in open(f"{d}/out") if l.endswith("yes\n")], ["X"])
    def test_order_patience(self):
        # four chains in a line, Hi-C decaying with distance: restarts stop after --patience failures, the order is found,
        # and a second run gives the same file
        import random; random.seed(3)
        with tempfile.TemporaryDirectory() as d:
            L = 200000; names = ["cA", "cB", "cC", "cD"]
            open(f"{d}/fai", "w").write("".join(f"{c}\t{L}\n" for c in names))
            open(f"{d}/groups", "w").write("#chain\tlength\tgroup\n" + "".join(f"{c}\t{L}\tg01\n" for c in names))
            with gzip.open(f"{d}/pairs.gz", "wt") as o:
                for _ in range(20000):
                    x = random.randrange(4 * L); dist = int(random.expovariate(1 / 60000)); y = x + dist
                    if y >= 4 * L: continue
                    o.write(f"{names[x // L]}\t{x % L}\t{names[y // L]}\t{y % L}\n")
            args = [f"{d}/groups", f"{d}/pairs.gz", f"{d}/fai"]
            r = run("hic_order.py", *args, f"{d}/o1", "--minPairs", "1"); self.assertEqual(r.returncode, 0, r.stderr)
            m = re.search(r"g01: restarts (\d+) of at most 40 \(patience 3\), improved (\d+)", r.stderr); self.assertTrue(m, r.stderr)
            if m.group(2) == "0": self.assertEqual(m.group(1), "3", "no restart improved, so three failures stop the search")
            got = [l.split("\t")[2] for l in open(f"{d}/o1") if l[0] != "#"]
            self.assertIn(got, (names, names[::-1]))
            run("hic_order.py", *args, f"{d}/o2", "--minPairs", "1")
            self.assertEqual(open(f"{d}/o1").read(), open(f"{d}/o2").read())
    def test_kept_mappings_keep_target(self):
        # a kept .paf.gz must take its .target with it, or a reset --remap leaves a stale checksum that stops the next copied-in mapping
        sys.path.insert(0, R); from agentssembler.stages import STAGES
        for st in STAGES:
            for k in st.keep:
                if ".paf.gz" in k: self.assertTrue(k.endswith("*"), f"{st.id}: {k}")
    def test_check_target(self):
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/t.fa", "w").write(">a\nACGT\n"); open(f"{d}/m.paf.gz", "w").write("x")
            sh = lambda: subprocess.run(["bash", "-c", f'source {R}/stages/lib.sh; check_target m.paf.gz t.fa; echo ok'], capture_output=True, text=True,
                                    env=dict(os.environ, AS_WORK=d, AS_REPO=R, AS_MINIMAP2="/bin/true", AS_SAMTOOLS="/bin/true", AS_SEQKIT="/bin/true"))
            self.assertEqual(sh().stdout.strip(), "ok")                               # no .target: a mapping made before checksums
            subprocess.run(["bash", "-c", f"md5sum < {d}/t.fa | cut -d' ' -f1 > {d}/m.paf.gz.target"])
            self.assertEqual(sh().stdout.strip(), "ok")
            open(f"{d}/t.fa", "w").write(">a\nACGTA\n")
            r = sh(); self.assertEqual(r.returncode, 4); self.assertIn("--remap", r.stderr)
    def test_tail_joins_filled_gap(self):
        # flank | tail | gap | other flank, the gap since filled: sequence on both sides of the tail, so the join is found by the flank's unitigs
        W = lambda c, k, x0, x1, u: f"{c}\t{x0 + 1}\t{x1}\t{k}\tW\t{u}\t1\t{x1 - x0}\t+\n"
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/c.agp", "w").write(W("cF", 1, 0, 10000, "uF") + W("cT", 1, 0, 5000, "uT") + W("cB", 1, 0, 10000, "uB"))
            open(f"{d}/n.tsv", "w").write("#g\tname\ng1\ts_h1\n")
            open(f"{d}/t.tsv", "w").write("#end\tgroup\tpath\tlabel\tpairs\tbp\ncF+\tg1\tcT+\tonly\t5\t5000\n")
            for agp, j in ((W("s_h1", 1, 0, 10000, "uF") + W("s_h1", 2, 10000, 15000, "uT") + W("s_h1", 3, 15000, 15500, "fill1") + W("s_h1", 4, 15500, 25500, "uB"), 10000),
                           (W("s_h1", 1, 0, 10000, "uB") + W("s_h1", 2, 10000, 10500, "fill1") + W("s_h1", 3, 10500, 15500, "uT") + W("s_h1", 4, 15500, 25500, "uF"), 15500)):
                open(f"{d}/s.agp", "w").write(agp)
                with gzip.open(f"{d}/r.paf.gz", "wt") as o: o.write(f"r1\t9000\t0\t9000\t+\ts_h1\t25500\t{j - 4000}\t{j + 4000}\t8000\t8000\t60\n")
                o = run("rawreads/tail_joins.py", "--tails", f"{d}/t.tsv", "--chainsAgp", f"{d}/c.agp", "--agp", f"{d}/s.agp", "--names", f"{d}/n.tsv", "--paf", f"{d}/r.paf.gz", "--out", f"{d}/o.tsv")
                self.assertEqual(o.returncode, 0, o.stderr)
                self.assertEqual([l.rstrip("\n").split("\t")[4:] for l in open(f"{d}/o.tsv") if l[0] != "#"], [[str(j), "1", "1"]])
    def test_gfa_overlap(self):
        sys.path.insert(0, f"{R}/scripts"); from chainlib import gfa_overlap
        self.assertEqual([gfa_overlap(x) for x in ("12345M", "10M1I5M", "10M2D5M", "*", "0M\n")], [12345, 16, 15, 0, 0])
    def test_chain_names_past_8000(self):
        # chain8* and chain9* mark connectors and restored arms: past 7999 chains the main chains take six digits
        with tempfile.TemporaryDirectory() as d:
            for n, want in ((12, "chain0012"), (8000, "chain008000")):
                open(f"{d}/fai", "w").write("".join(f"u{i}\t{100000 - i}\n" for i in range(n))); open(f"{d}/keep", "w").write("".join(f"u{i}\n" for i in range(n)))
                open(f"{d}/gfa", "w").write("")
                r = run("build_chains.py", "--keep", f"{d}/keep", "--gfa", f"{d}/gfa", "--fai", f"{d}/fai", "--out", f"{d}/c"); self.assertEqual(r.returncode, 0, r.stderr)
                names = {l.split("\t")[0] for l in open(f"{d}/c.agp") if l[0] != "#"}
                self.assertIn(want, names); self.assertFalse(any(x.startswith(("chain8", "chain9")) for x in names))
    def test_utg_names(self):
        # an edit names a chain by a unitig it holds: the name survives renumbering, and an end follows the unitig's own strand
        sys.path.insert(0, f"{R}/scripts"); from chainlib import utg_resolver
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/c.agp", "w").write("chain0001\t1\t10\t1\tW\tuA\t1\t10\t+\nchain8003\t1\t10\t1\tW\tuB\t1\t10\t-\nchain0002\t1\t10\t1\tW\tuC\t1\t10\t+\nchain0002\t11\t20\t2\tW\tuA\t1\t10\t+\n")
            res = utg_resolver(f"{d}/c.agp")
            self.assertEqual([res(x) for x in ("utg:uB", "utg:uB+", "utg:uB.g03-", "utg:uC+", "chain0001+")], ["chain8003", "chain8003-", "chain8003.g03+", "chain0002+", "chain0001+"])
            for bad in ("utg:uA", "utg:uZ"): self.assertRaises(SystemExit, res, bad)                # in two chains; in none
            open(f"{d}/g.tsv", "w").write("#chain\tlength\tgroup\nchain0001\t10\tg01\nchain8003\t10\tunassigned\nchain0002\t20\tg02\n")
            open(f"{d}/g.chroms.tsv", "w").write("g01\tc0\ng02\tc0\n")
            open(f"{d}/e", "w").write("utg:uB\tgroup\twith:utg:uC\tconnector between chain0002's arms\n")
            r = run("apply_group_edits.py", f"{d}/g", f"{d}/e", f"{d}/f", f"{d}/c.agp"); self.assertEqual(r.returncode, 0, r.stderr)
            self.assertIn("chain8003\t10\tg02", open(f"{d}/f.tsv").read())
            open(f"{d}/l", "w").write("#a\tb\tkind\treason\nutg:uB+\tchain0002-\thand\tread evidence\n")
            r = run("resolve_utg.py", f"{d}/c.agp", f"{d}/l", "--keep", "3"); self.assertEqual(r.stdout, "chain8003-\tchain0002-\thand\n")
    def test_decay_norm(self):
        sys.path.insert(0, f"{R}/scripts"); import hic_order as ho
        # one chain ending just inside the lower half of a populated bin (the review's case): the old midpoint opportunity there was 0,
        # clamped to 1, and the extrapolated slope turned positive. Contact must not rise with distance past the data.
        edges = np.logspace(2, 9, 141); k = int(np.searchsorted(edges, 2e6)); L = edges[k] * 1.02
        d = []                                                       # contact ~ 1/x times the positions x apart, (L - x): its counts per bin
        for a, b in zip(edges[:-1], edges[1:]):
            if a >= L: break
            b = min(b, L); n = round(0.15 * (L * np.log(b / a) - (b - a)))
            d += list(np.linspace(a, b, n, endpoint=False))
        self.assertGreaterEqual(sum(x >= edges[k] for x in d), 50)
        lp = ho.decay_model(d, lengths=[L])
        self.assertLessEqual(lp(np.array([1e8]))[0], lp(np.array([1e6]))[0])
        self.assertTrue(np.all(np.isfinite(ho.decay_model([500, 600, 700])(np.array([1e3, 1e6])))))   # sparse: no bin with 50 pairs
        self.assertEqual(ho.decay_model([])(np.array([1e3]))[0], 0.0)
        x = np.array([1e3, 1e5, 1e7]); self.assertTrue(np.array_equal(ho.decay_model(d)(x), ho.decay_model([], hist=ho.decay_hist(d))(x)))
    def test_hic_order_decay_file(self):
        # the decay model written once (--decayOut) and read per chromosome (--decayIn) gives the same order as building it in the job
        import random; random.seed(3)
        with tempfile.TemporaryDirectory() as d:
            L = {"a": 400000, "b": 300000, "c": 200000}
            open(f"{d}/fai", "w").write("".join(f"{c}\t{n}\n" for c, n in L.items()))
            open(f"{d}/g", "w").write("".join(f"{c}\t{n}\tg1\n" for c, n in L.items()))
            off = {"a": 0, "b": 400000, "c": 700000}; rows = []
            for _ in range(20000):
                p = random.randrange(900000); q = min(899999, max(0, p + int(random.expovariate(1 / 30000)) * random.choice((-1, 1))))
                cx = max(c for c in L if off[c] <= p); cy = max(c for c in L if off[c] <= q)
                rows.append(f"{cx}\t{p - off[cx] + 1}\t{cy}\t{q - off[cy] + 1}\n")
            with gzip.open(f"{d}/p.gz", "wt") as o: o.writelines(rows)
            r = run("hic_order.py", f"{d}/g", f"{d}/p.gz", f"{d}/fai", "-", "--decayOut", f"{d}/decay"); self.assertEqual(r.returncode, 0, r.stderr)
            for args, out in (([], "o1"), (["--decayIn", f"{d}/decay"], "o2")):
                r = run("hic_order.py", f"{d}/g", f"{d}/p.gz", f"{d}/fai", f"{d}/{out}", "--minLen", "1000", *args); self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(open(f"{d}/o1").read(), open(f"{d}/o2").read())
            self.assertEqual([l.split("\t")[2] for l in open(f"{d}/o1") if l[0] != "#"], ["a", "b", "c"])
    def test_unhonoured_filter(self):
        # a same-group link below make_order_links.py's bar (< 5 kb overlap, < 3 reads) is not counted as unhonoured
        with tempfile.TemporaryDirectory() as d:
            p = lambda n, x: open(f"{d}/{n}", "w").write(x)
            p("o", "#g\trank\tchain\to\tlen\ng1\t1\tx\t+\t10\ng1\t2\ty\t+\t10\ng1\t3\tz\t+\t10\n")
            p("b", "#a\tb\tsame_group\tpath\n"); p("t", ""); p("ch", "g1\tc1\n")
            for ov, n in (("1000", "count 0"), ("6000", "count 2")):
                p("e", f"#a\tb\treads\tq\tov\nx+\tz-\t0\tNA\t{ov}\n")
                r = run("review/unhonoured.py", f"{d}/o", f"{d}/e", f"{d}/b", f"{d}/t", f"{d}/ch"); self.assertEqual(r.returncode, 0, r.stderr)
                self.assertIn(n, r.stdout)
    def test_first_undecided(self):
        sys.path.insert(0, f"{R}/scripts"); from chainlib import first_undecided
        dec = lambda o0, o1: o0 >= 10 + o1 and o0 >= 2 * max(o1, 1)
        # the review's case: equal reads (200), Hi-C 5 against 30. Reads + Hi-C (205 against 230) would not decide and keep the first.
        self.assertEqual(first_undecided([((200, 5), "p1"), ((200, 30), "p2")], dec), "p2")
        self.assertEqual(first_undecided([((200, 5), "p1"), ((200, 9), "p2")], dec), "p1")          # undecided: the first
        self.assertEqual(first_undecided([((100, 5), "p1"), ((180, 20), "p2")], dec), "p1")         # reads differ: 200 does not decide 105
        self.assertEqual(first_undecided([((10, 5), "p1"), ((180, 30), "p2")], dec), "p2")
    def test_hub_cores(self):
        # Each chain: 300 kb, 1x arms (unique 20x) around a 2x core of unique reads at 100-200 kb, except C (an allele bubble: 1x, MAPQ 0)
        # and D (classed shared). The core's unitigs are those at least half inside the run; breaks go at their outer boundaries.
        #   A: unitigs 0-150, 150-300 kb: neither is half inside the core, so it lies inside unitigs and is not breakable
        #   E: unitigs 0-100, 100-200, 200-300 kb: the middle one is the core; breakable at 100000 and 200000
        #   F: unitigs 0-60, 60-200, 200-300 kb: the core unitig carries 40 kb of 1x arm (core_2x_frac 0.71): not breakable
        with tempfile.TemporaryDirectory() as d:
            U = {"A": [150, 300], "C": [300], "D": [300], "E": [100, 200, 300], "F": [60, 200, 300]}
            with open(f"{d}/w.tsv", "w") as o:
                o.write("#unitig\tstart\tend\ttotal\tunique\tsecondary\n")
                for c in U:
                    for i in range(30):
                        core = c != "C" and 10 <= i < 20
                        o.write(f"{c}\t{i * 10000}\t{(i + 1) * 10000}\t{160 if core or c == 'D' else 80}\t{160 if core else 40 if c == 'C' else 20}\t0\n")
            with open(f"{d}/c.agp", "w") as o:
                for c, ends in U.items():
                    x = 0
                    for k, e in enumerate(ends):
                        o.write(f"{c}\t{x + 1}\t{e * 1000}\t{k + 1}\tW\t{c}u{k}\t1\t{e * 1000 - x}\t+\n"); x = e * 1000
            open(f"{d}/k.tsv", "w").write("#h\n" + "".join(f"{c}\t300000\t{160 if c == 'D' else 80}\t0.3\t0.7\t0.00\t{'shared' if c == 'D' else 'specific'}\n" for c in U))
            o = run("review/hub_cores.py", f"{d}/w.tsv", f"{d}/c.agp", f"{d}/k.tsv", f"{d}/o.tsv", "--sharedLo", "129")
            self.assertEqual(o.returncode, 0, o.stderr)
            F = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(f"{d}/o.tsv") if l[0] != "#"}
            self.assertEqual(sorted(F), ["A", "E", "F"])
            self.assertEqual(F["A"][4:7] + F["A"][10:], ["100000", "200000", "100", ".", ".", "inside/inside", "1.00", "no"])
            self.assertEqual(F["E"][10:], ["100000", "200000", "unitig/unitig", "1.00", "yes"])
            self.assertEqual(F["F"][10:], ["60000", "200000", "unitig/unitig", "0.71", "no"])
            # before D04: class from the windows ("-"), positions lifted onto chains that D03 split (E broken at 50 kb into E and Eb)
            open(f"{d}/l.agp", "w").write("E\t1\t50000\t1\tW\tEu0\t1\t50000\t+\nEb\t1\t50000\t1\tW\tEu0\t50001\t100000\t+\n"
                                       "Eb\t50001\t150000\t2\tW\tEu1\t1\t100000\t+\nEb\t150001\t250000\t3\tW\tEu2\t1\t100000\t+\n")
            o = run("review/hub_cores.py", f"{d}/w.tsv", f"{d}/c.agp", "-", f"{d}/o.tsv", "--sharedLo", "129", "--liftAgp", f"{d}/l.agp")
            self.assertEqual(o.returncode, 0, o.stderr)
            F = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(f"{d}/o.tsv") if l[0] != "#"}
            self.assertEqual(F["Eb"][4:6] + F["Eb"][10:], ["50000", "150000", "50000", "150000", "unitig/unitig", "1.00", "yes"])
            self.assertEqual(F["Eb"][2:4], ["specific", "80.0"])      # median depth from the windows, in its own column
            # E lifted onto chains stored the other way round (R: a 50 kb unitig, then E reversed), split right of the core, and split on both sides
            W = lambda c, k, x0, x1, u, st: f"{c}\t{x0 + 1}\t{x1}\t{k}\tW\t{u}\t1\t{x1 - x0}\t{st}\n"
            for agp, want in (
                    (W("R", 1, 0, 50000, "X", "+") + W("R", 2, 50000, 150000, "Eu2", "-") + W("R", 3, 150000, 250000, "Eu1", "-") + W("R", 4, 250000, 350000, "Eu0", "-"),
                     ("R", ["150000", "250000", "150000", "250000", "unitig/unitig", "1.00", "yes"])),
                    (W("E1", 1, 0, 100000, "Eu0", "+") + W("E1", 2, 100000, 200000, "Eu1", "+") + W("E2", 1, 0, 100000, "Eu2", "+"),
                     ("E1", ["100000", "200000", "100000", ".", "unitig/end", "1.00", "yes"])),
                    (W("E1", 1, 0, 100000, "Eu0", "+") + W("E2", 1, 0, 100000, "Eu1", "+") + W("E3", 1, 0, 100000, "Eu2", "+"),
                     ("E2", ["0", "100000", ".", ".", "end/end", "1.00", "no"]))):
                open(f"{d}/l.agp", "w").write(agp)
                o = run("review/hub_cores.py", f"{d}/w.tsv", f"{d}/c.agp", "-", f"{d}/o.tsv", "--sharedLo", "129", "--liftAgp", f"{d}/l.agp")
                self.assertEqual(o.returncode, 0, o.stderr)
                F = {l.split("\t")[0]: l.rstrip("\n").split("\t") for l in open(f"{d}/o.tsv") if l[0] != "#"}
                self.assertEqual(F[want[0]][4:6] + F[want[0]][10:], want[1])
    def test_rescue_collapsed_redundant(self):
        # a1 (g1) and a2 (g2) anchor chromosome c1. x has a breakable 2x unique-read core: collapsed. y is covered in both groups by short
        # allele blocks (each under --minConflict): redundant. z is placed. w links only to x: it still finds c1 through x.
        with tempfile.TemporaryDirectory() as d:
            L = {"a1": 2000000, "a2": 2000000, "x": 300000, "y": 120000, "z": 200000, "w": 150000}
            open(f"{d}/g.tsv", "w").write("#chain\tlength\tgroup\n" + "".join(f"{c}\t{n}\t{'g1' if c == 'a1' else 'g2' if c == 'a2' else 'unassigned'}\n" for c, n in L.items()))
            open(f"{d}/g.chroms.tsv", "w").write("#group\tchrom\ng1\tc1\ng2\tc1\n")
            with gzip.open(f"{d}/r.tsv.gz", "wt") as o:
                for c, n in (("x", 40), ("y", 30), ("z", 30), ("w", 30)):
                    for i in range(n): o.write(f"{c}\t1\t{'x' if c == 'w' else 'a1' if i % 2 else 'a2'}\t1\n")
            with gzip.open(f"{d}/u.tsv.gz", "wt") as o:
                for i in range(30): o.write("z\t1\ta1\t1\n")
            open(f"{d}/al.tsv", "w").write("#h\n")
            open(f"{d}/rep.tsv", "w").write("#h\n" + "".join(f"{c}\t{n}\t0.0\n" for c, n in L.items()))
            open(f"{d}/b.tsv", "w").write("#h\ny\t0\t40000\ta1\t5000\t45000\t0.0\t+\ny\t40000\t70000\ta1\t9000\t39000\t0.0\t+\ny\t0\t70000\ta2\t1000\t71000\t0.0\t+\n")
            open(f"{d}/h.tsv", "w").write("#chain\n" + "x\t300000\tspecific\t80\t100000\t200000\t100\t160\t160\t0.25\t100000\t200000\tunitig/unitig\t1.00\tyes\n")
            o = run("rescue/rescue_groups.py", f"{d}/g", f"{d}/r.tsv.gz", f"{d}/u.tsv.gz", f"{d}/al.tsv", f"{d}/rep.tsv", "--hubCores", f"{d}/h.tsv", "--blocks", f"{d}/b.tsv")
            self.assertEqual(o.returncode, 0, o.stderr)
            C = {l.split("\t")[0]: (l.split("\t")[3], l.split("\t")[4]) for l in open(f"{d}/g.rescued.tsv") if l[0] != "#"}
            self.assertEqual(C["x"], ("unassigned", "collapsed")); self.assertEqual(C["y"], ("unassigned", "redundant"))
            self.assertEqual(C["z"], ("g1", "hic")); self.assertIn(C["w"][0], ("g1", "g2"))
            # without the options, the old calls
            shutil.copy(f"{d}/g.hic.tsv", f"{d}/g.tsv")
            o = run("rescue/rescue_groups.py", f"{d}/g", f"{d}/r.tsv.gz", f"{d}/u.tsv.gz", f"{d}/al.tsv", f"{d}/rep.tsv")
            self.assertNotIn("collapsed", open(f"{d}/g.rescued.tsv").read())
    def test_rescue_graph_tiebreak(self):
        # x has rescued links to chromosome c1 and no unique Hi-C to either haplotype: a graph link to a2 alone decides it (call graph);
        # without --endLinks, or with links to both haplotypes, it is arbitrary
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/g.chroms.tsv", "w").write("#group\tchrom\ng1\tc1\ng2\tc1\n"); open(f"{d}/al.tsv", "w").write("#h\n")
            open(f"{d}/rep.tsv", "w").write("#h\na1\t2000000\t0.0\na2\t2000000\t0.0\nx\t300000\t0.0\n")
            with gzip.open(f"{d}/r.tsv.gz", "wt") as o: o.write("".join(f"x\t1\t{'a1' if i % 2 else 'a2'}\t1\n" for i in range(40)))
            with gzip.open(f"{d}/u.tsv.gz", "wt") as o: o.write("")
            def call(links):
                open(f"{d}/g.tsv", "w").write("#chain\tlength\tgroup\na1\t2000000\tg1\na2\t2000000\tg2\nx\t300000\tunassigned\n")
                open(f"{d}/el.tsv", "w").write("#end_a\tend_b\treads\tq_gap_median\tgraph_overlap\n" + links)
                r = run("rescue/rescue_groups.py", f"{d}/g", f"{d}/r.tsv.gz", f"{d}/u.tsv.gz", f"{d}/al.tsv", f"{d}/rep.tsv", "--endLinks", f"{d}/el.tsv"); self.assertEqual(r.returncode, 0, r.stderr)
                return {l.split("\t")[0]: tuple(l.split("\t")[3:5]) for l in open(f"{d}/g.rescued.tsv") if l[0] != "#"}["x"]
            self.assertEqual(call("x+\ta2-\t0\tNA\t9000\n"), ("g2", "graph"))
            self.assertEqual(call("x+\ta2-\t0\tNA\t9000\nx-\ta1+\t0\tNA\t9000\n")[1], "arbitrary")
            self.assertEqual(call("x+\ta2-\t1\tNA\t900\n")[1], "arbitrary")                     # too weak a link
    def test_allele_table_dovetails(self):
        # chains A, B, C, D of one unitig each (100 kb). A's end = B's start, identical: a dovetail, set aside. A and C share an identical
        # start on the same strand: parallel alleles (a homozygous start), kept. A's start = D's start on opposite strands: a dovetail.
        # A and B also align at 1% in the middle: an allele block, kept. A near-full copy (C and D, 90 kb) stays an allele too.
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/c.agp", "w").write("".join(f"{c}\t1\t100000\t1\tW\tu{c}\t1\t100000\t+\n" for c in "ABCD"))
            def paf(q, qs, qe, sd, t, ts, te, de): return f"u{q}\t100000\t{qs}\t{qe}\t{sd}\tu{t}\t100000\t{ts}\t{te}\t{qe-qs}\t{qe-qs}\t60\tde:f:{de}\n"
            with gzip.open(f"{d}/s.paf.gz", "wt") as o:
                o.write(paf("A", 70000, 100000, "+", "B", 0, 30000, 0.0) + paf("A", 0, 30000, "+", "C", 0, 30000, 0.0)
                        + paf("A", 0, 25000, "-", "D", 0, 25000, 0.0) + paf("A", 35000, 60000, "+", "B", 35000, 60000, 0.01)
                        + paf("C", 5000, 95000, "+", "D", 5000, 95000, 0.0))
            r = run("allele_table.py", f"{d}/c.agp", f"{d}/s.paf.gz", f"{d}/al"); self.assertEqual(r.returncode, 0, r.stderr)
            rows = lambda f: [l.split("\t") for l in open(f) if l[0] != "#"]
            self.assertEqual(sorted((s[0], s[3], s[1]) for s in rows(f"{d}/al.dovetails.tsv")), [("A", "B", "70000"), ("A", "D", "0")])
            self.assertEqual(sorted((s[0], s[3], s[1]) for s in rows(f"{d}/al.blocks.tsv")), [("A", "B", "35000"), ("A", "C", "0"), ("C", "D", "5000")])
            P = {(s[0], s[1]): int(s[2]) for s in rows(f"{d}/al.pairs.tsv")}
            self.assertEqual(P[("A", "B")], 25000); self.assertNotIn(("A", "D"), P)
    def test_rescue_dovetail_link(self):
        # x's end overlaps a2's end identically (a dovetail, not an allele): it is a graph link to g2, so x goes with a2 (call graph)
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/g.chroms.tsv", "w").write("#group\tchrom\ng1\tc1\ng2\tc1\n"); open(f"{d}/al.tsv", "w").write("#h\n")
            open(f"{d}/rep.tsv", "w").write("#h\na1\t2000000\t0.0\na2\t2000000\t0.0\nx\t300000\t0.0\n")
            with gzip.open(f"{d}/r.tsv.gz", "wt") as o: o.write("".join(f"x\t1\t{'a1' if i % 2 else 'a2'}\t1\n" for i in range(40)))
            with gzip.open(f"{d}/u.tsv.gz", "wt") as o: o.write("")
            open(f"{d}/g.tsv", "w").write("#chain\tlength\tgroup\na1\t2000000\tg1\na2\t2000000\tg2\nx\t300000\tunassigned\n")
            open(f"{d}/dv.tsv", "w").write("#h\na2\t1940000\t2000000\tx\t0\t60000\t0.0\t+\n")
            r = run("rescue/rescue_groups.py", f"{d}/g", f"{d}/r.tsv.gz", f"{d}/u.tsv.gz", f"{d}/al.tsv", f"{d}/rep.tsv", "--dovetails", f"{d}/dv.tsv")
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual({l.split("\t")[0]: tuple(l.split("\t")[3:5]) for l in open(f"{d}/g.rescued.tsv") if l[0] != "#"}["x"], ("g2", "graph"))
            # x's other end overlaps a1 as well: it sits where the haplotypes meet, so it stays unassigned (call fork)
            open(f"{d}/g.tsv", "w").write("#chain\tlength\tgroup\na1\t2000000\tg1\na2\t2000000\tg2\nx\t300000\tunassigned\n")
            open(f"{d}/dv.tsv", "a").write("x\t250000\t300000\ta1\t0\t50000\t0.0\t+\n")
            r = run("rescue/rescue_groups.py", f"{d}/g", f"{d}/r.tsv.gz", f"{d}/u.tsv.gz", f"{d}/al.tsv", f"{d}/rep.tsv", "--dovetails", f"{d}/dv.tsv")
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual({l.split("\t")[0]: tuple(l.split("\t")[3:5]) for l in open(f"{d}/g.rescued.tsv") if l[0] != "#"}["x"], ("unassigned", "fork"))
    def test_rescue_dovetail_fork_later(self):
        # x (longer, so called first) runs into a2's end: graph, g2. z, shorter, is placed in g1 by unique Hi-C later in the round, and its end
        # runs into x's start. x now links to both groups: it sits where the haplotypes meet, so it ends unassigned (fork)
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/g.chroms.tsv", "w").write("#group\tchrom\ng1\tc1\ng2\tc1\n"); open(f"{d}/al.tsv", "w").write("#h\n")
            open(f"{d}/rep.tsv", "w").write("#h\na1\t2000000\t0.0\na2\t2000000\t0.0\nx\t300000\t0.0\nz\t200000\t0.0\n")
            with gzip.open(f"{d}/r.tsv.gz", "wt") as o:
                o.write("".join(f"{c}\t1\t{'a1' if i % 2 else 'a2'}\t1\n" for c in ("x", "z") for i in range(40)))
            with gzip.open(f"{d}/u.tsv.gz", "wt") as o: o.write("".join("z\t1\ta1\t1\n" for i in range(30)))
            open(f"{d}/g.tsv", "w").write("#chain\tlength\tgroup\na1\t2000000\tg1\na2\t2000000\tg2\nx\t300000\tunassigned\nz\t200000\tunassigned\n")
            open(f"{d}/dv.tsv", "w").write("#h\na2\t1940000\t2000000\tx\t250000\t300000\t0.0\t-\nz\t150000\t200000\tx\t0\t50000\t0.0\t+\n")
            r = run("rescue/rescue_groups.py", f"{d}/g", f"{d}/r.tsv.gz", f"{d}/u.tsv.gz", f"{d}/al.tsv", f"{d}/rep.tsv", "--dovetails", f"{d}/dv.tsv")
            self.assertEqual(r.returncode, 0, r.stderr)
            C = {l.split("\t")[0]: tuple(l.split("\t")[3:5]) for l in open(f"{d}/g.rescued.tsv") if l[0] != "#"}
            self.assertEqual(C["z"], ("g1", "hic")); self.assertEqual(C["x"], ("unassigned", "fork"))
    def test_rescue_bubble_arm(self):
        # y (g1) runs into x's start (y's end dovetails x, g1). b's start dovetails into the same start of x, and b's end is identical to
        # y's interior: b is y's other-haplotype arm, so g2 (call bubble_arm), not g1 by the graph link
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/g.chroms.tsv", "w").write("#group\tchrom\ng1\tc1\ng2\tc1\n"); open(f"{d}/al.tsv", "w").write("#h\n")
            open(f"{d}/rep.tsv", "w").write("#h\ny\t2000000\t0.0\nx\t500000\t0.0\nz\t2000000\t0.0\nb\t180000\t0.0\n")
            with gzip.open(f"{d}/r.tsv.gz", "wt") as o: o.write("".join(f"b\t1\t{'y' if i % 2 else 'z'}\t1\n" for i in range(40)))
            with gzip.open(f"{d}/u.tsv.gz", "wt") as o: o.write("")
            open(f"{d}/g.tsv", "w").write("#chain\tlength\tgroup\ny\t2000000\tg1\nx\t500000\tg1\nz\t2000000\tg2\nb\t180000\tunassigned\n")
            open(f"{d}/dv.tsv", "w").write("#h\ny\t1960000\t2000000\tx\t0\t40000\t0.0\t+\nb\t0\t50000\tx\t0\t50000\t0.0\t-\n")
            open(f"{d}/bl.tsv", "w").write("#h\ny\t1830000\t1862000\tb\t148000\t180000\t0.0\t+\n")
            r = run("rescue/rescue_groups.py", f"{d}/g", f"{d}/r.tsv.gz", f"{d}/u.tsv.gz", f"{d}/al.tsv", f"{d}/rep.tsv", "--dovetails", f"{d}/dv.tsv", "--blocks", f"{d}/bl.tsv")
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual({l.split("\t")[0]: tuple(l.split("\t")[3:5]) for l in open(f"{d}/g.rescued.tsv") if l[0] != "#"}["b"], ("g2", "bubble_arm"))
    def test_rescue_copy(self):
        # c is identical (de 0) over its whole length to part of a1, an allele by the table. Together they read 0.75x: one haplotype's
        # sequence twice, so c is left unassigned (call copy). At 1.8x together it is a homozygous stretch, and the allele call stands
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/g.chroms.tsv", "w").write("#group\tchrom\ng1\tc1\ng2\tc1\n")
            open(f"{d}/al.tsv", "w").write("#h\na1\tc\t150000\t150000\t0.0\t1\n")
            open(f"{d}/rep.tsv", "w").write("#h\na1\t2000000\t0.0\na2\t2000000\t0.0\nc\t150000\t0.0\n")
            with gzip.open(f"{d}/r.tsv.gz", "wt") as o: o.write("".join(f"c\t1\t{'a1' if i % 2 else 'a2'}\t1\n" for i in range(40)))
            with gzip.open(f"{d}/u.tsv.gz", "wt") as o: o.write("")
            open(f"{d}/bl.tsv", "w").write("#h\na1\t500000\t650000\tc\t0\t150000\t0.0\t+\n"); open(f"{d}/dv.tsv", "w").write("#h\n")
            for depth, want in ((20, ("unassigned", "copy")), (50, ("g2", "allele"))):
                open(f"{d}/g.tsv", "w").write("#chain\tlength\tgroup\na1\t2000000\tg1\na2\t2000000\tg2\nc\t150000\tunassigned\n")
                open(f"{d}/w.tsv", "w").write("#h\n" + "".join(f"a1\t{x}\t{x + 10000}\t{depth + 10}\nc\t{x - 500000}\t{x - 490000}\t{depth - 5}\n" for x in range(500000, 650000, 10000)))
                r = run("rescue/rescue_groups.py", f"{d}/g", f"{d}/r.tsv.gz", f"{d}/u.tsv.gz", f"{d}/al.tsv", f"{d}/rep.tsv", "--dovetails", f"{d}/dv.tsv",
                        "--blocks", f"{d}/bl.tsv", "--depthWindows", f"{d}/w.tsv", "--hap", "60")
                self.assertEqual(r.returncode, 0, r.stderr)
                self.assertEqual({l.split("\t")[0]: tuple(l.split("\t")[3:5]) for l in open(f"{d}/g.rescued.tsv") if l[0] != "#"}["c"], want)
    def test_rescue_dovetails_redundant(self):
        # y is a2's end plus b2's start (dovetails into both, g2) and 60% an allele of a1 (g1): all its sequence is placed, so redundant
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/g.chroms.tsv", "w").write("#group\tchrom\ng1\tc1\ng2\tc1\n"); open(f"{d}/al.tsv", "w").write("#h\n")
            open(f"{d}/rep.tsv", "w").write("#h\na1\t2000000\t0.0\na2\t2000000\t0.0\nb2\t2000000\t0.0\ny\t100000\t0.0\n")
            with gzip.open(f"{d}/r.tsv.gz", "wt") as o: o.write("".join(f"y\t1\t{'a1' if i % 2 else 'a2'}\t1\n" for i in range(40)))
            with gzip.open(f"{d}/u.tsv.gz", "wt") as o: o.write("")
            open(f"{d}/g.tsv", "w").write("#chain\tlength\tgroup\na1\t2000000\tg1\na2\t2000000\tg2\nb2\t2000000\tg2\ny\t100000\tunassigned\n")
            open(f"{d}/dv.tsv", "w").write("#h\na2\t1940000\t2000000\ty\t0\t60000\t0.0\t+\ny\t40000\t100000\tb2\t0\t60000\t0.0\t+\n")
            open(f"{d}/bl.tsv", "w").write("#h\na1\t500000\t560000\ty\t20000\t80000\t0.02\t+\n")
            r = run("rescue/rescue_groups.py", f"{d}/g", f"{d}/r.tsv.gz", f"{d}/u.tsv.gz", f"{d}/al.tsv", f"{d}/rep.tsv", "--dovetails", f"{d}/dv.tsv", "--blocks", f"{d}/bl.tsv")
            self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual({l.split("\t")[0]: tuple(l.split("\t")[3:5]) for l in open(f"{d}/g.rescued.tsv") if l[0] != "#"}["y"], ("unassigned", "redundant"))
    def test_rescue_swaps(self):
        # c_h1: Z - A | gap | Y - B; c_h2: Z2 - P - X - Q. A links only X; P, B and Q link both X and Y: swapping X and Y closes the gap.
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/j.tsv", "w").write("#h\nc_h1\t50\tZ+\tA+\tdirect\nc_h1\t100\tA+\tY+\thic\nc_h1\t200\tY+\tB+\tdirect\n"
                                       "c_h2\t50\tZ2+\tP+\tdirect\nc_h2\t100\tP+\tX+\tdirect\nc_h2\t200\tX+\tQ+\tdirect\n")
            links = "#h\nZ+\tA-\t0\tNA\t50000\nZ2+\tP-\t0\tNA\t50000\nA+\tX-\t0\tNA\t50000\nX+\tB-\t0\tNA\t50000\nP+\tY-\t0\tNA\t50000\nY+\tQ-\t0\tNA\t50000\nP+\tX-\t0\tNA\t50000\nX+\tQ-\t0\tNA\t50000\nY+\tB-\t0\tNA\t50000\n"
            open(f"{d}/l.tsv", "w").write(links); open(f"{d}/b.tsv", "w").write("#h\n")
            o = run("review/rescue_swaps.py", f"{d}/j.tsv", f"{d}/l.tsv", f"{d}/b.tsv", f"{d}/o.tsv")
            self.assertEqual(o.returncode, 0, o.stderr)
            F = [l.split("\t") for l in open(f"{d}/o.tsv") if l[0] != "#"]
            self.assertEqual([(r[4], r[6], r[7]) for r in F], [("Y", "X", "c_h2")])
            open(f"{d}/l.tsv", "w").write(links.replace("P+\tY-\t0\tNA\t50000\n", ""))   # Y cannot take X's place: no swap
            run("review/rescue_swaps.py", f"{d}/j.tsv", f"{d}/l.tsv", f"{d}/b.tsv", f"{d}/o.tsv")
            self.assertEqual([l for l in open(f"{d}/o.tsv") if l[0] != "#"], [])
    def test_build_curated_overlap(self):
        # A ends with the 30 bp that B starts with, and the graph has no link: overlap:30 joins them with no gap, a wrong length does not
        import random; random.seed(1); q = lambda n: "".join(random.choice("ACGT") for _ in range(n))
        x, sh, y = q(70), q(30), q(70)
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/u.fa", "w").write(f">uA\n{x + sh}\n>uB\n{sh + y}\n")
            open(f"{d}/c.agp", "w").write("A\t1\t100\t1\tW\tuA\t1\t100\t+\nB\t1\t100\t1\tW\tuB\t1\t100\t+\n")
            open(f"{d}/o.tsv", "w").write("#h\ng1\t1\tA\t+\t.\t.\t.\t0\ng1\t2\tB\t+\t.\t.\t.\t0\n")
            open(f"{d}/l.tsv", "w").write("#h\n"); open(f"{d}/ch.tsv", "w").write("#h\n")
            open(f"{d}/s.tsv", "w").write("A\t100\t1000\nB\t100\t1000\n"); open(f"{d}/n.tsv", "w").write("g1\tchr1\tkeep\n")
            args = ["--order", f"{d}/o.tsv", "--agp", f"{d}/c.agp", "--links", f"{d}/l.tsv", "--choices", f"{d}/ch.tsv", "--stats", f"{d}/s.tsv",
                    "--names", f"{d}/n.tsv", "--unitigs", f"{d}/u.fa", "--out", f"{d}/b", "--curated", f"{d}/cur.tsv"]
            for n, want in ((30, x + sh + y), (31, x + sh + "N" * 100 + sh + y)):
                open(f"{d}/cur.tsv", "w").write(f"#h\nA+\tB-\toverlap:{n}\treason\n")
                o = run("build_scaffolds.py", *args)
                self.assertEqual(o.returncode, 0, o.stderr)
                self.assertEqual(open(f"{d}/b.fasta").read().split("\n")[1], want)
    def test_homolog_direction(self):
        # u1 and u2 align forward (unitig strand +), but u2 lies reversed in s_h2: s_h2 is written in reverse against s_h1
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/a.agp", "w").write("s_h1\t1\t50000\t1\tW\tu1\t1\t50000\t+\ns_h2\t1\t50000\t1\tW\tu2\t1\t50000\t-\n"
                                       "t_h1\t1\t50000\t1\tW\tu3\t1\t50000\t+\nt_h2\t1\t50000\t1\tW\tu4\t1\t50000\t+\n")
            with gzip.open(f"{d}/p.paf.gz", "wt") as o:
                o.write("u1\t50000\t0\t40000\t+\tu2\t50000\t0\t40000\t39000\t40000\t60\n")
                o.write("u3\t50000\t0\t40000\t+\tu4\t50000\t0\t40000\t39000\t40000\t60\n")
            o = run("review/homolog_direction.py", f"{d}/a.agp", f"{d}/p.paf.gz", f"{d}/o.tsv")
            self.assertEqual(o.returncode, 0, o.stderr)
            R = {l.split("\t")[0]: l.rstrip("\n").split("\t")[4] for l in open(f"{d}/o.tsv") if l[0] != "#"}
            self.assertEqual(R, {"s_h2": "inverted", "t_h2": "forward"})
    def test_stable_key(self):
        # a renumbered chain keeps its key; copies keep their suffix; a chain missing from the AGP falls back to its name
        sys.path.insert(0, f"{R}/scripts"); import chainlib
        P = {"chain0007": [("utg9", "+"), ("utg2", "-")], "chain8001": [("utg5", "+")]}
        Q = {"chain0012": [("utg9", "+"), ("utg2", "-")], "chain8003": [("utg5", "+")]}
        self.assertEqual(chainlib.stable_key(P, "chain0007.g05"), chainlib.stable_key(Q, "chain0012.g05"))
        self.assertEqual(chainlib.stable_path(P, ["chain8001+", "chain0007-"]), chainlib.stable_path(Q, ["chain8003+", "chain0012-"]))
        self.assertEqual(chainlib.stable_key(P, "chainX"), (("chainX",), ""))
if __name__ == "__main__": unittest.main(verbosity=1)
