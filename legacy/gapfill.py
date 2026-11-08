#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Haplotype-aware gap filling of a pseudochromosome AGP from the hifiasm unitig graph,
with Hi-C deciding between alternative alleles.

Every component stays a hifiasm unitig (or a trimmed sub-range of one), so the result is
still an AGP over the unitigs and the rest of the pipeline (agp_to_fasta.py, contact maps,
Juicebox) works on it unchanged. No sequence is invented and no reads are needed.

For each gap between two unitigs L and R (in their scaffold orientations):

1. Paths in the p_utg GFA from the end of L to the start of R are enumerated: the direct
   link L -> R, and paths through up to --maxNodes interior unitigs adding <= --maxBp.
   Paths through a unitig placed on the same scaffold or on a scaffold of another locus are
   not used (they point to a misorder and are reported).
2. Each interior unitig X is scored by its allele-specific Hi-C links:
     own   = links to the unitigs of this scaffold within +-W bp of the gap
     other = links to the unitigs of the best competing scaffold at the same reference
             locus (unitigs whose reference footprint overlaps the gap's reference interval)
   both normalised as sum(links / sqrt(depth_X * depth_w)) like hic_phase.py.
3. The gap is filled with:
     direct     L and R overlap in the graph: trim the overlap. No new sequence.
     graph      exactly one path, and Hi-C does not contradict it (other < 2 x own, or
                too few links to say)
     hic        several paths; the chosen one's distinguishing unitigs have >= --minLinks
                own links, prefer this scaffold >= --ratio x over any competitor, and beat
                the alternative's distinguishing unitigs >= --ratio x
   otherwise the gap stays ('ambiguous', 'hic_rejects', 'no_path', ...).
   Gaps with no graph path but an exact overlap of the flank sequences (>= --minAlnOverlap
   bp, identity >= --minAlnIdentity, found with minimap2) are closed as 'alignment'.
4. An unused unitig chosen by gaps on two scaffolds is given to the one whose Hi-C prefers
   it (>= --ratio x); if neither does, both keep it and it is reported as 'shared'
   (homozygous sequence present in both haplotypes).

Trimming follows the graph: each node keeps the sequence beyond what the previous nodes
already spell, so link overlaps are never duplicated.

Outputs (--out prefix): .agp, .report.tsv (one row per gap), .membership.tsv (new unitigs
and how they were placed), .unphased.bed (fills without independent Hi-C phase support).
'''
import sys, argparse, collections, subprocess, os, tempfile, math, threading, concurrent.futures
from common import natural_key

def flip(s): return "-" if s == "+" else "+"

def homolog_key(scaf):
    """chr5_h2 -> chr5: scaffolds sharing this key are homologs"""
    return scaf.rsplit("_", 1)[0]

def read_fai(f): return { l.split("\t")[0]: int(l.split("\t")[1]) for l in open(f) }

def read_gfa_links(f):
    adj = collections.defaultdict(dict) # (u, us) -> {(v, vs): overlap}
    for l in open(f):
        if l[0] != "L": continue
        s = l.split("\t"); u, us, v, vs = s[1], s[2], s[3], s[4]
        ov = int(s[5].rstrip("M")) if s[5].rstrip("M").isdigit() else 0
        adj[(u, us)][(v, vs)] = ov
        adj[(v, flip(vs))][(u, flip(us))] = ov
    return adj

def read_links(f):
    links = collections.defaultdict(dict); depth = {}
    for l in open(f):
        a, b, n = l.rstrip("\n").split("\t")
        if b == ".": depth[a] = int(n)
        else: links[a][b] = int(n); links[b][a] = int(n)
    return links, depth

def read_footprints(f):
    '''annotarium phasing.tsv: main reference region of each unitig'''
    fp = {}
    for l in open(f):
        if l.startswith("#"): continue
        s = l.rstrip("\n").split("\t")
        if len(s) < 8 or s[6] in ("", "."): continue
        regs, bps = s[6].split(","), [ int(x) for x in s[7].split(",") ]
        c, rng = regs[bps.index(max(bps))].rsplit(":", 1)
        a, b = rng.split("-"); fp[s[0]] = (c, int(a), int(b))
    return fp

def read_agp(f):
    rows = []
    for l in open(f):
        if l.startswith("#"): continue
        rows.append(l.rstrip("\n").split("\t"))
    return rows

def shared_units(note):
    """unitigs named in 'shared:utgA@chr1_h1,utgB@...' or 'shared:utgA' notes"""
    out = set()
    for w in note.split():
        if w.startswith("shared:"):
            out.update(x.split("@")[0] for x in w[7:].split(","))
    return out

def rc(s): return s[::-1].translate(str.maketrans("ACGTNacgtn", "TGCANtgcan"))

class Seqs:
    '''random access to the unitig FASTA through its .fai'''
    def __init__(self, fasta):
        self.f = open(fasta, "rb"); self.idx = {}; self.lock = threading.Lock()
        for l in open(fasta + ".fai"):
            s = l.split("\t"); self.idx[s[0]] = (int(s[1]), int(s[2]), int(s[3]), int(s[4]))
    def get(self, name, start, end): # 0-based half-open
        ln, off, lb, by = self.idx[name]
        start = max(0, start); end = min(ln, end)
        a = off + (start // lb) * by + start % lb
        b = off + (end // lb) * by + end % lb
        with self.lock:
            self.f.seek(a); raw = self.f.read(b - a)
        return raw.decode().replace("\n", "").replace("\r", "")

def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--agp", required=True)
    p.add_argument("--gfa", required=True, help="hifiasm *.p_utg.gfa (the graph of the unitigs in the AGP)")
    p.add_argument("--fasta", required=True, help="unitig FASTA (with .fai)")
    p.add_argument("--links", required=True, help="unitig-unitig Hi-C links + depths (03_links/links.tsv)")
    p.add_argument("--footprints", required=True, help="annotarium phasing.tsv (reference footprint of every unitig)")
    p.add_argument("--out", required=True)
    p.add_argument("--maxNodes", type=int, default=3)
    p.add_argument("--maxBp", type=int, default=1500000)
    p.add_argument("--maxPaths", type=int, default=200)
    p.add_argument("--window", type=int, default=1000000, help="Hi-C window around the gap (bp)")
    p.add_argument("--minLinks", type=int, default=10)
    p.add_argument("--ratio", type=float, default=2.0)
    p.add_argument("--minimap2", default=None, help="minimap2 binary: also close exact flank overlaps with no graph link")
    p.add_argument("--flank", type=int, default=200000)
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--minAlnOverlap", type=int, default=10000)
    p.add_argument("--minAlnIdentity", type=float, default=0.9999)
    p.add_argument("--minEquivIdentity", type=float, default=0.999, help="alternative fills this similar are interchangeable")
    a = p.parse_args()

    ln = read_fai(a.fasta + ".fai")
    adj = read_gfa_links(a.gfa)
    links, depth = read_links(a.links)
    fp = read_footprints(a.footprints)
    rows = read_agp(a.agp)

    def spell(lu, ls, path, ru, rs):
        """Pieces spelled by a graph path between the end of (lu, ls) and the start of (ru, rs):
        each node keeps only the sequence beyond what earlier nodes already spell.
        Returns the interior pieces [unitig, strand, start, end] (1-based) and the number of
        bp to trim from the start of the right unitig."""
        nodes = [(lu, ls)] + list(path) + [(ru, rs)]
        end_sofar = ln[lu]; start = 0; pieces = []; skipR = 0
        for (u0, s0), (u1, s1) in zip(nodes[:-1], nodes[1:]):
            start = start + ln[u0] - adj[(u0, s0)][(u1, s1)]
            skip = max(0, end_sofar - start)
            if (u1, s1) == nodes[-1]:
                skipR = skip
            elif skip < ln[u1]:
                pieces.append([u1, s1, skip + 1, ln[u1]] if s1 == "+" else [u1, s1, 1, ln[u1] - skip])
            end_sofar = max(end_sofar, start + ln[u1])
        return pieces, skipR

    # Current placement of every unitig, and component positions per scaffold
    used = collections.defaultdict(set)
    comps = collections.defaultdict(list) # scaffold -> [(start, end, unitig)]
    for r in rows:
        if r[4] == "W":
            used[r[5]].add(r[0]); comps[r[0]].append((int(r[1]), int(r[2]), r[5]))

    def locus(scaf, pos):
        '''reference interval of the components of scaf within +-window of pos'''
        hits = [ fp[u] for s, e, u in comps[scaf] if e >= pos - a.window and s <= pos + a.window and u in fp ]
        if not hits: return None
        c = collections.Counter(h[0] for h in hits).most_common(1)[0][0]
        xs = [ h for h in hits if h[0] == c ]
        return c, min(h[1] for h in xs), max(h[2] for h in xs)

    def window_units(scaf, pos):
        return { u for s, e, u in comps[scaf] if e >= pos - a.window and s <= pos + a.window }

    def competitor_units(scaf, loc):
        '''unitigs of every other scaffold whose footprint overlaps the locus'''
        out = collections.defaultdict(set)
        if loc is None: return out
        c, lo, hi = loc
        for s2, cs in comps.items():
            if s2 == scaf: continue
            for _, _, u in cs:
                f = fp.get(u)
                if f and f[0] == c and f[2] >= lo and f[1] <= hi: out[s2].add(u)
        return out

    def score(x, units):
        raw = 0; norm = 0.0; dx = depth.get(x, 0)
        for w in units:
            n = links[x].get(w, 0)
            if n and dx and depth.get(w, 0):
                raw += n; norm += n / math.sqrt(dx * depth[w])
        return raw, norm

    gaps = [ i for i, r in enumerate(rows) if r[4] in "UN" ]
    decisions = {} # row index -> dict
    claims = collections.defaultdict(list) # unused unitig -> [gap row idx]

    for i in gaps:
        r, L, R = rows[i], rows[i-1], rows[i+1]
        gid = f"{r[0]}:{r[1]}"; scaf = r[0]; pos = int(r[1])
        d = dict(gid=gid, scaf=scaf, left=L[5]+L[8], right=R[5]+R[8], method=".", path=[], added=0, own_links=0, other_links=0,
                 own=0.0, other=0.0, competitor=".", n_paths=0, note="")
        decisions[i] = d
        if L[4] != "W" or R[4] != "W":
            d["method"] = "not_between_components"; continue
        lu, ls, ru, rs = L[5], L[8], R[5], R[8]
        # Graph links apply only to unitig ends; skip if the AGP already trimmed the facing ends
        lfull = (int(L[7]) == ln[lu]) if ls == "+" else (int(L[6]) == 1)
        rfull = (int(R[6]) == 1) if rs == "+" else (int(R[7]) == ln[ru])
        if not (lfull and rfull):
            d["method"] = "flank_not_unitig_end"; continue

        # enumerate paths
        paths = []
        if (ru, rs) in adj[(lu, ls)]: paths.append([])
        stack = [ ((lu, ls), [], 0) ]
        conflict = set()
        while stack and len(paths) < a.maxPaths:
            node, path, added = stack.pop()
            for (v, vs), ov in adj[node].items():
                if (v, vs) == (ru, rs):
                    if path: paths.append(path)
                    continue
                if v in (lu, ru) or any(v == x for x, _ in path) or len(path) >= a.maxNodes: continue
                na = added + ln[v] - ov
                if na > a.maxBp: continue
                stack.append(((v, vs), path + [(v, vs)], na))
        loc = locus(scaf, pos)
        comp = competitor_units(scaf, loc)
        own_units = window_units(scaf, pos)
        valid = []
        for path in paths:
            # a placed interior unitig is allowed only on a homolog (same chromosome prefix): shared sequence
            bad = [ v for v, _ in path if v in used and (scaf in used[v] or any(homolog_key(s2) != homolog_key(scaf) for s2 in used[v])) ]
            if bad: conflict.update(f"{v}@{','.join(sorted(used[v]))}" for v in bad); continue
            valid.append(path)
        d["n_paths"] = len(valid)
        if conflict: d["note"] = "path_via_placed:" + ";".join(sorted(conflict))[:300]

        # Hi-C score of each node
        node_score = {}
        def ns(v):
            if v not in node_score:
                o = score(v, own_units)
                best = (0, 0.0, ".")
                for s2, us in comp.items():
                    sc = score(v, us)
                    if sc[1] > best[1]: best = (sc[0], sc[1], s2)
                node_score[v] = (o[0], o[1], best[0], best[1], best[2])
            return node_score[v]
        def agg(nodes):
            o = sum(ns(v)[0] for v in nodes); on = sum(ns(v)[1] for v in nodes)
            ot = sum(ns(v)[2] for v in nodes); otn = sum(ns(v)[3] for v in nodes)
            comps_ = { ns(v)[4] for v in nodes if ns(v)[4] != "." }
            return o, on, ot, otn, ",".join(sorted(comps_)) or "."

        if not valid:
            d["method"] = "no_path"; continue
        uniq = { tuple(v for v, _ in pth): pth for pth in valid }
        valid = list(uniq.values())
        valid.sort(key=lambda pth: (sum(v in used for v, _ in pth), -sum(ns(v)[1] for v, _ in pth), len(pth)))
        best = valid[0]
        bn = [ v for v, _ in best ]
        o, on, ot, otn, cp = agg(bn)
        d.update(path=best, own_links=o, other_links=ot, own=on, other=otn, competitor=cp)
        if not best:
            # direct link: if longer alternatives exist, they must not be clearly better
            d["method"] = "direct"; continue
        if len(valid) == 1:
            if ot >= a.minLinks and otn >= a.ratio * on:
                d["method"] = "hic_rejects"; continue
            d["method"] = "graph"
        else:
            alt = [ v for v, _ in valid[1] ]
            d1 = [ v for v in bn if v not in alt ]; d2 = [ v for v in alt if v not in bn ]
            o1 = agg(d1); o2 = agg(d2)
            ok = (o1[0] >= a.minLinks and o1[1] >= a.ratio * o1[3] and o1[1] >= a.ratio * o2[1])
            if not ok:
                d["method"] = "ambiguous"; d["alts"] = valid
                d["note"] = (d["note"] + " " if d["note"] else "") + f"alt={';'.join(alt) or 'direct'}:own_links={o2[0]}"
                continue
            d["method"] = "hic"

    # Ambiguous gaps whose alternative fills are the same sequence (hifiasm often keeps
    # near-identical copies of homozygous sequence as separate unitigs; they carry no
    # allele-specific Hi-C, and choosing either cannot create a mosaic)
    if a.minimap2:
        seqs = Seqs(a.fasta)
        def oriented(u, s, start, end): # 0-based on the oriented sequence
            if s == "+": return seqs.get(u, start, end)
            n = ln[u]; return rc(seqs.get(u, n - end, n - start))
        def piece_seq(pc):
            u, s, st, en = pc
            return oriented(u, s, st - 1, en) if s == "+" else oriented(u, s, ln[u] - en, ln[u] - st + 1)
        def fill_seqs(i):
            d = decisions[i]; L, R = rows[i-1], rows[i+1]
            sp = [ spell(L[5], L[8], pth, R[5], R[8]) for pth in d["alts"] ]
            q = max(s for _, s in sp) + 1000 # compare up to a common point in the right unitig
            return [ "".join(piece_seq(pc) for pc in pcs) + oriented(R[5], R[8], skip, q) for pcs, skip in sp ]
        def identity(x, y, tmp):
            if x == y: return 1.0
            if min(len(x), len(y)) < 1000: return 0.0
            open(f"{tmp}/x.fa", "w").write(f">x\n{x}\n"); open(f"{tmp}/y.fa", "w").write(f">y\n{y}\n")
            out = subprocess.run([a.minimap2, "-x", "asm5", "-c", f"{tmp}/y.fa", f"{tmp}/x.fa"], capture_output=True, text=True).stdout
            m = sum(int(l.split("\t")[9]) for l in out.splitlines() if "tp:A:P" in l)
            return m / max(len(x), len(y))
        def equiv(i):
            fs = fill_seqs(i)
            with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(a.out))) as tmp:
                return i, min(identity(fs[0], f, tmp) for f in fs[1:]), len(fs)
        amb = [ i for i in gaps if decisions[i]["method"] == "ambiguous" ]
        with concurrent.futures.ThreadPoolExecutor(a.threads) as ex:
            for i, idn, n in ex.map(equiv, amb):
                d = decisions[i]
                d["note"] = (d["note"] + " " if d["note"] else "") + f"alt_identity={idn:.5f}"
                if idn >= a.minEquivIdentity:
                    d["method"] = "equivalent"

    # Each unused unitig fills at most one gap unless that is impossible. Gaps are taken in
    # order of evidence (hic, graph, equivalent). A later gap that wants a taken unitig
    # switches to an equivalent alternative path if it has one, loses it to the holder if the
    # unitig's Hi-C clearly prefers the holder's scaffold, and otherwise shares it (a
    # homozygous copy, reported).
    rank = { "hic": 0, "graph": 1, "equivalent": 2 }
    taken = {} # unitig -> gap row index
    for i in sorted((i for i in gaps if decisions[i]["method"] in rank), key=lambda i: (rank[decisions[i]["method"]], i)):
        d = decisions[i]
        def clash(path): return [ v for v, _ in path if v not in used and v in taken ]
        if clash(d["path"]) and d["method"] == "equivalent":
            free = [ pth for pth in d["alts"] if not clash(pth) and not any(v in used for v, _ in pth) ]
            if free:
                d["path"] = free[0]; d["note"] = (d["note"] + " " if d["note"] else "") + "switched_to_free_twin"
        for v in clash(d["path"]):
            h = taken[v]
            s_h = score(v, window_units(decisions[h]["scaf"], int(rows[h][1])))
            s_i = score(v, window_units(d["scaf"], int(rows[i][1])))
            if s_h[0] >= a.minLinks and s_h[1] >= a.ratio * s_i[1]:
                d["method"] = "lost_to_" + decisions[h]["scaf"]; d["path"] = []; break
            for j in (h, i):
                decisions[j]["note"] = (decisions[j]["note"] + " " if decisions[j]["note"] else "") + f"shared:{v}"
        if d["method"] in rank:
            for v, _ in d["path"]:
                if v not in used: taken.setdefault(v, i)
            homs = [ v for v, _ in d["path"] if v in used ]
            if homs and "shared:" + homs[0] not in d["note"]:
                d["note"] = (d["note"] + " " if d["note"] else "") + "shared:" + ",".join(f"{v}@{','.join(sorted(used[v]))}" for v in homs)

    # Exact flank overlaps without a graph link
    if a.minimap2:
        todo = [ i for i in gaps if decisions[i]["method"] == "no_path" ]
        def align(i):
            L, R = rows[i-1], rows[i+1]
            lu, ls, ru, rs = L[5], L[8], R[5], R[8]
            fl = min(a.flank, ln[lu]); fr = min(a.flank, ln[ru])
            with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(a.out))) as tmp:
                open(f"{tmp}/L.fa", "w").write(f">L\n{oriented(lu, ls, ln[lu]-fl, ln[lu])}\n")
                open(f"{tmp}/R.fa", "w").write(f">R\n{oriented(ru, rs, 0, fr)}\n")
                out = subprocess.run([a.minimap2, "-x", "asm5", "-c", "-N", "20", f"{tmp}/R.fa", f"{tmp}/L.fa"],
                                     capture_output=True, text=True).stdout
            best = None
            for l in out.splitlines():
                s = l.split("\t")
                qlen, qs, qe, st, ts, te, m, bl = int(s[1]), int(s[2]), int(s[3]), s[4], int(s[7]), int(s[8]), int(s[9]), int(s[10])
                if st == "+" and qlen - qe <= 100 and ts <= 100 and m / bl >= a.minAlnIdentity and bl >= a.minAlnOverlap:
                    if best is None or bl > best[0]: best = (bl, qlen - qe, te, m / bl)
            return i, best
        with concurrent.futures.ThreadPoolExecutor(a.threads) as ex:
            for i, best in ex.map(align, todo):
                if best:
                    d = decisions[i]
                    d["method"] = "alignment"; d["aln"] = best
                    d["note"] = (d["note"] + " " if d["note"] else "") + f"aln_overlap={best[0]}:identity={best[3]:.5f}"

    # Build the new AGP. Each scaffold is a list of oriented pieces [unitig, strand, start, end]
    # (1-based on the unitig) and gaps; a closed gap splices the path in, trimming overlaps.
    def take_suffix(piece, n):
        '''drop the first n bp of a piece in its scaffold orientation'''
        u, s, st, en = piece
        if s == "+": return [u, s, st + n, en]
        return [u, s, st, en - n]
    def drop_end(piece, n):
        u, s, st, en = piece
        if s == "+": return [u, s, st, en - n]
        return [u, s, st + n, en]

    out_rows = collections.defaultdict(list)
    order = []
    i = 0
    filled_methods = { "direct", "graph", "hic", "equivalent", "alignment" }
    bed = []; members = []
    by_scaf = collections.defaultdict(list)
    for j, r in enumerate(rows):
        by_scaf[r[0]].append(j)
    for scaf in sorted(by_scaf, key=natural_key):
        items = []
        for j in by_scaf[scaf]:
            r = rows[j]
            if r[4] == "W": items.append(["W", [r[5], r[8], int(r[6]), int(r[7])], j])
            else: items.append(["G", r, j])
        merged = []
        k = 0
        while k < len(items):
            it = items[k]
            if it[0] == "G" and decisions[it[2]]["method"] in filled_methods and merged and merged[-1][0] == "W":
                d = decisions[it[2]]
                left = merged[-1][1]; right = items[k+1][1]
                if d["method"] == "alignment":
                    _, dropL, te, _ = d["aln"]
                    left_new = drop_end(left, dropL); right_new = take_suffix(right, te)
                    merged[-1][1] = left_new
                    new = []
                else:
                    new, skipR = spell(left[0], left[1], d["path"], right[0], right[1])
                    right_new = take_suffix(right, skipR)
                    d["fill_bp"] = sum(p[3] - p[2] + 1 for p in new)
                L_ok = merged[-1][1][3] >= merged[-1][1][2]
                R_ok = right_new[3] >= right_new[2]
                if not (L_ok and R_ok):
                    d["method"] = "trim_exceeds_component"; merged.append(it); k += 1; continue
                for pc in new:
                    merged.append(["W", pc, None, d])
                items[k+1][1] = right_new
                d["filled"] = True
                k += 1
                continue
            merged.append(it); k += 1
        # write
        posn = 1; part = 1
        for it in merged:
            if it[0] == "W":
                u, s, st, en = it[1]; L = en - st + 1
                out_rows[scaf].append([scaf, posn, posn + L - 1, part, "W", u, st, en, s])
                if len(it) > 3 and it[3] is not None:
                    d = it[3]
                    sh = u in shared_units(d["note"])
                    if d["method"] != "hic" or sh:
                        bed.append((scaf, posn - 1, posn + L - 1, f"{u}:{d['method']}{':shared' if sh else ''}"))
                    members.append((u, scaf, "gapfill_" + d["method"] + ("_shared" if sh else "")))
                posn += L
            else:
                r = it[1]; L = int(r[2]) - int(r[1]) + 1
                out_rows[scaf].append([scaf, posn, posn + L - 1, part] + r[4:])
                posn += L
            part += 1

    with open(a.out + ".agp", "w") as o:
        o.write("##agp-version 2.1\n")
        for scaf in sorted(out_rows, key=natural_key):
            for r in out_rows[scaf]: o.write("\t".join(map(str, r)) + "\n")
    with open(a.out + ".report.tsv", "w") as o:
        o.write("#gap_id\tleft\tright\tmethod\tfilled\tn_paths\tpath\tfill_bp\town_links\tcompetitor_links\town_score\tcompetitor_score\tcompetitor\tnote\n")
        for i in gaps:
            d = decisions[i]
            o.write(f"{d['gid']}\t{d['left']}\t{d['right']}\t{d['method']}\t{'yes' if d.get('filled') else 'no'}\t{d['n_paths']}\t"
                    f"{';'.join(v + s for v, s in d['path']) or '.'}\t{d.get('fill_bp', 0)}\t{d['own_links']}\t{d['other_links']}\t"
                    f"{d['own']:.4f}\t{d['other']:.4f}\t{d['competitor']}\t{d['note'] or '.'}\n")
    with open(a.out + ".membership.tsv", "w") as o:
        o.write("#contig\tscaffold\thow\n")
        for m in members: o.write("\t".join(m) + "\n")
    with open(a.out + ".unphased.bed", "w") as o:
        for b in bed: o.write("\t".join(map(str, b)) + "\n")
    c = collections.Counter((decisions[i]["method"], bool(decisions[i].get("filled"))) for i in gaps)
    print("#method\tfilled\tgaps")
    for (m, f), n in sorted(c.items(), key=lambda x: -x[1]): print(f"{m}\t{'yes' if f else 'no'}\t{n}")
    print(f"# gaps filled: {sum(n for (m, f), n in c.items() if f)} of {len(gaps)}")

if __name__ == "__main__":
    main()
