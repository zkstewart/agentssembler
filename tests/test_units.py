#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
"""Fast checks with no data: every script parses and states its purpose, no hard-coded paths or samples, and the edit-applying scripts behave. Run: python3 tests/test_units.py"""
import ast, glob, gzip, os, re, shutil, subprocess, sys, tempfile, textwrap, unittest
R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = glob.glob(f"{R}/scripts/*.py") + glob.glob(f"{R}/scripts/*/*.py")
def run(script, *a): return subprocess.run([sys.executable, f"{R}/scripts/{script}", *a], capture_output=True, text=True)
class T(unittest.TestCase):
    def test_parse_and_doc(self):
        for f in SCRIPTS:
            t = ast.parse(open(f).read())
            if not f.endswith("chainlib.py"): self.assertTrue(ast.get_docstring(t), f)
    def test_no_hard_coding(self):
        pat = re.compile(r"/home/|research/genomes|13Q010|WatsonsChoice|00_inputs|16_validate/|\bBusk\b")
        for f in SCRIPTS + glob.glob(f"{R}/stages/*.sh"):
            if f.endswith("lib.sh"): continue
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
    def test_agp_to_fasta(self):
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/u.fa", "w").write(">a\nAACC\n>b\nGGTT\n")
            open(f"{d}/a.agp", "w").write("s\t1\t4\t1\tW\ta\t1\t4\t+\ns\t5\t7\t2\tU\t3\tscaffold\tyes\tx\ns\t8\t11\t3\tW\tb\t1\t4\t-\n")
            r = run("agp_to_fasta.py", f"{d}/a.agp", f"{d}/u.fa", f"{d}/o.fa"); self.assertEqual(r.returncode, 0, r.stderr)
            self.assertEqual(open(f"{d}/o.fa").read(), ">s\nAACCNNNAACC\n")
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
            open(f"{d}/e2", "w").write("s\tshared_in\tg01,g05\ttwo chromosomes\n"); self.assertNotEqual(run("apply_group_edits.py", f"{d}/g", f"{d}/e2", f"{d}/x").returncode, 0)
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
    def test_homolog_and_tail_alarms(self):
        sys.path.insert(0, R)
        from agentssembler import report as r
        with tempfile.TemporaryDirectory() as d:
            open(f"{d}/f.fai", "w").write("c1_h1\t30000000\t0\t60\t61\nc1_h2\t25000000\t0\t60\t61\nc2_h1\t30000000\t0\t60\t61\nc2_h2\t29500000\t0\t60\t61\n")
            open(f"{d}/a.agp", "w").write("c1_h2\t1\t100\t1\tW\tu\t1\t100\t+\nc1_h2\t101\t200\t2\tN\t100\tscaffold\tyes\tna\n")
            A = r.homolog_lengths(f"{d}/f.fai", f"{d}/a.agp")
            self.assertEqual(len(A), 1); self.assertIn("c1_h2", A[0]); self.assertIn("1 gap:", A[0])
            open(f"{d}/t.tsv", "w").write("#h\nc1_h2\t9\tright\t14\tunplaced_a\t14\t1\t1\t2\t+\t60\nc1_h2\t9\tleft\t10\tunplaced_b\t2\t1\t1\t2\t+\t60\nc1_h2\t9\tleft\t4\tc1_h2\t4\t1\t1\t2\t+\t60\nc1_h2\t7\tleft\t4\tc1_h1\t4\t1\t1\t2\t+\t60\nc2_h1\t5\tleft\t3\t.\t0\t0\t.\t.\t.\t.\n")
            T = r.tail_alarms(f"{d}/t.tsv")
            self.assertEqual(len(T), 2); self.assertIn("14 of 14", T[0]); self.assertIn("1 gap edges", T[1])
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
            open(f"{d}/all.fa", "w").write(f">s\n{L}{'N' * 100}{Rt}\n>un\n{U}\n"); open(f"{d}/r.fa", "w").write(f">r1\n{read}\n")
            with gzip.open(f"{d}/m.paf.gz", "wt") as o: o.write("r1\t13000\t7000\t13000\t+\ts\t40100\t20100\t26100\t6000\t6000\t60\n")
            x = run("review/gap_tails.py", f"{d}/a.agp", f"{d}/r.fa", f"{d}/all.fa", f"{d}/o.tsv", "--paf", f"{d}/m.paf.gz", "--minimap2", os.environ.get("MINIMAP2") or shutil.which("minimap2"), "--threads", "2")
            self.assertEqual(x.returncode, 0, x.stderr)
            F = [l.split("\t") for l in open(f"{d}/o.tsv") if l[0] != "#"]
            self.assertEqual(len(F), 1); self.assertEqual(F[0][:6], ["s", "20000", "right", "1", "un", "1"])
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
if __name__ == "__main__": unittest.main(verbosity=1)
