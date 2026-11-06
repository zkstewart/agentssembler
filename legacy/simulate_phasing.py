#!/usr/bin/env python3
'''
Test hic_phase.py on a planted phasing with this assembly's real reference placement.

1. A random valid phasing (no allelic conflicts) is made by running hic_phase.py on random
   links; it is the "truth".
2. Hi-C links are simulated from the truth: expected links between two unitigs on the same
   chromosome decay with the distance d between their reference positions as
   lam = 50 * (1e6 / (d + 1e6))^1.5, scaled by --ratio when they are on different
   haplotypes (e.g. 0.33 = 3:1 within vs between homologs).
3. hic_phase.py is run on the simulated links, and the result is scored as the share of
   unitig pairs within 3 Mb whose same/different-haplotype relationship matches the truth
   (label permutations do not matter). Random labelling scores about 0.56.
'''
import os, re, sys, random, argparse, subprocess

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--paf", required=True); p.add_argument("--phasing", required=True)
    p.add_argument("--ploidy", type=int, default=3)
    p.add_argument("--ratio", type=float, default=0.33)
    p.add_argument("--seed", type=int, default=7)
    p.add_argument("--workdir", required=True)
    args = p.parse_args()
    os.makedirs(args.workdir, exist_ok=True)
    phaser = os.path.join(os.path.dirname(os.path.abspath(__file__)), "hic_phase.py")
    W = lambda f: os.path.join(args.workdir, f)
    run = lambda links, out, beam: subprocess.run([sys.executable, phaser, "--paf", args.paf, "--phasing", args.phasing,
        "--links", links, "--out", out, "--beam", str(beam), "--ploidy", str(args.ploidy)], check=True, stdout=subprocess.DEVNULL)

    random.seed(args.seed)
    ids = [ l.split("\t")[0] for l in open(args.phasing) if "\tassigned\t" in l ]
    with open(W("random_links.tsv"), "w") as f:
        for a in ids: f.write(f"{a}\t.\t1\n")
        for _ in range(200000):
            a, b = random.sample(ids, 2); f.write(f"{a}\t{b}\t{random.randint(1, 5)}\n")
    run(W("random_links.tsv"), W("truth"), 50)

    truth, pos = {}, {}
    for l in open(W("truth.phasing.tsv")):
        r = l.rstrip("\n").split("\t")
        if r[3] == "assigned":
            truth[r[0]] = r[2]
            m = re.match(r"([^:]+):(\d+)-(\d+)", r[6]); pos[r[0]] = (m.group(1), (int(m.group(2)) + int(m.group(3))) // 2)
    random.seed(args.seed + 1)
    tids = list(truth)
    with open(W("simulated_links.tsv"), "w") as out:
        for i, a in enumerate(tids):
            out.write(f"{a}\t.\t{random.randint(500, 2000)}\n")
            for b in tids[i+1:]:
                if pos[a][0] != pos[b][0]: continue
                d = abs(pos[a][1] - pos[b][1])
                lam = 50 * (1e6 / (d + 1e6)) ** 1.5 * (1.0 if truth[a] == truth[b] else args.ratio)
                n = sum(1 for _ in range(3) if random.random() < 0.5) * int(random.expovariate(1 / lam))
                if n: out.write(f"{a}\t{b}\t{n}\n")
    run(W("simulated_links.tsv"), W("recovered"), 500)

    for name, f in (("input phasing", args.phasing), ("recovered", W("recovered.phasing.tsv"))):
        lab = {}
        for l in open(f):
            r = l.rstrip("\n").split("\t")
            if not l.startswith("#") and r[0] in truth and r[2] != ".": lab[r[0]] = r[2]
        q = list(lab); ok = tot = 0
        for i, a in enumerate(q):
            for b in q[i+1:]:
                if pos[a][0] == pos[b][0] and abs(pos[a][1] - pos[b][1]) < 3e6:
                    tot += 1; ok += ((lab[a] == lab[b]) == (truth[a] == truth[b]))
        print(f"{name}: pairwise same/different-haplotype agreement with the planted truth (pairs < 3 Mb apart) = {ok/tot:.3f} of {tot}")

if __name__ == "__main__":
    main()
