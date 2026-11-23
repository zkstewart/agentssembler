#!/usr/bin/env python3
# This Source Code Form is subject to the terms of the Mozilla Public License, v. 2.0.
# If a copy of the MPL was not distributed with this file, You can obtain one at https://mozilla.org/MPL/2.0/.
'''
Reference-free grouping of chains into haplotypes from Hi-C and allelic conflicts.

Score between chains = Hi-C links / sqrt(depth_a * depth_b). Two chains conflict (cannot share a
haplotype) when alleles.pairs.tsv shows >= --minConflict bp and >= --minConflictFrac of the shorter.
  1. Chains >= --seedLen are clustered into --chroms chromosomes by average linkage on Hi-C.
  2. Within each, average-linkage clustering into --ploidy haplotypes, never across a conflict, then
     single-chain moves while the same-haplotype link density improves.
  3. Shorter chains join the best non-conflicting group if they have >= --minLinks links there and
     >= --ratio times the next best; otherwise they stay unassigned.
Chains whose within-chain Hi-C density is below --minWithin (see chain_stats.py; --stats) do not seed
steps 1-2: their reads rarely place uniquely, so they carry little linkage and their "allelic" hits are
often repeats. They are only placed in step 3.

Shared chains (--shared, one chain per line): sequence collapsed over all haplotypes of a chromosome, found from read
depth. Their Hi-C comes from every haplotype, so they help step 1 but would merge haplotypes in steps 2-3. They are kept
out of those steps, then assigned to a chromosome: >= --minLinks links, and >= --ratio times the next chromosome by
links per group depth. Written to out_prefix.shared.tsv, with out_prefix.chroms.tsv (group to chromosome). A shared chain
is never in out_prefix.tsv as a member of a group: see make_shared_copies.py.

Dosage (--classDepth with --hap): above ploidy 2 a shared chain need not be in every haplotype. Its dosage is its median depth
over haplotype depth, from 2 to the number of haplotype groups of its chromosome. A chain below that number goes to
out_prefix.dosage.tsv with the groups chosen to carry it, ranked by: a graph link to a member (--endLinks), then no allele of
the chain among the members, then Hi-C links per group depth. The call column says what separated the chosen groups from
the rest: graph, alleles, hic (>= --ratio times the best group left out), or arbitrary. No file is written when every shared
chain is in every haplotype, which is always so at ploidy 2.

--hapFirst replaces steps 1-2: all seeds are clustered straight to chroms x ploidy haplotype groups under the conflict
constraint, and groups are then gathered into chromosomes, --ploidy at a time, by allele alignment. Each aligned chain pair
counts its allelic bp times (1 - divergence)^--deWeight, so that homologs outrank homeologs of similar aligned length.
Use it when homologs are divergent: Hi-C between homologs can then be as weak as between chromosomes, so forcing the
chromosome count first splits one chromosome by haplotype and merges two others. The merge scores around the target are
printed: a cliff right after the target is independent support for it. Tested at ploidy 2 only.
Usage: hic_group.py links.tsv alleles.pairs.tsv chains.fai out_prefix --chroms N --ploidy N [--stats chain_stats.tsv] [--shared shared.txt]
         [--hapFirst] [--classDepth chain_class.tsv --hap H [--endLinks end_links.tsv]]
'''
import argparse, math, collections
p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("links"); p.add_argument("alleles"); p.add_argument("fai"); p.add_argument("out")
p.add_argument("--chroms",type=int,required=True); p.add_argument("--ploidy",type=int,required=True); p.add_argument("--shared"); p.add_argument("--hapFirst",action="store_true"); p.add_argument("--seedLen",type=int,default=1000000)
p.add_argument("--maxConflictDe",type=float,help="ignore allele pairs whose median divergence (alleles.pairs.tsv median_de) is above this: in a genome with old duplications, paralogs align like alleles and would forbid true merges. Off by default; the citrus genomes did not need it"); p.add_argument("--minConflict",type=int,default=50000); p.add_argument("--minConflictFrac",type=float,default=0.1)
p.add_argument("--stats"); p.add_argument("--minWithin",type=float,default=700)
p.add_argument("--repeat",help="chain_repeat.tsv: conflicts involving a chain above --maxRepeat are ignored, because satellite arrays align to each other without being alleles")
p.add_argument("--maxRepeat",type=float,default=0.5)
p.add_argument("--minLinks",type=int,default=10); p.add_argument("--ratio",type=float,default=2.0)
p.add_argument("--deWeight",type=float,default=10,help="--hapFirst: allelic bp of a chain pair count times (1 - median divergence) to this power")
p.add_argument("--classDepth",help="chain_class.tsv: median depth per chain, for the dosage of shared chains"); p.add_argument("--hap",type=float,help="haplotype depth, with --classDepth")
p.add_argument("--endLinks",help="end_links.tsv: graph links from a shared chain to a group's members rank that group first for its dosage"); p.add_argument("--minOverlap",type=int,default=5000)
a=p.parse_args()
if a.classDepth and not a.hap: p.error("--classDepth needs --hap")
L={l.split("\t")[0]:int(l.split("\t")[1]) for l in open(a.fai)}
lk=collections.defaultdict(dict); dep={}
for l in open(a.links):
    x,y,n=l.rstrip().split("\t")
    if y==".": dep[x]=int(n)
    else: lk[x][y]=int(n); lk[y][x]=int(n)
w=lambda x,y: lk[x].get(y,0)/math.sqrt(max(dep.get(x,1),1)*max(dep.get(y,1),1))
conf=collections.defaultdict(set); sat=set()
if a.repeat:
    sat={l.split("\t")[0] for l in open(a.repeat) if l[0]!="#" and float(l.split("\t")[2])>a.maxRepeat}
    print(f"# {len(sat)} chains above repeat fraction {a.maxRepeat}: {sum(L[c] for c in sat)/1e6:.1f} Mb; their conflicts are ignored")
for l in open(a.alleles):
    if l[0]=="#": continue
    s=l.split("\t"); x,y,bx,by=s[0],s[1],int(s[2]),int(s[3])
    if x in sat or y in sat: continue
    if a.maxConflictDe is not None and len(s)>4 and s[4].strip()!="nan" and float(s[4])>a.maxConflictDe: continue
    if min(bx,by)>=a.minConflict and min(bx,by)>=a.minConflictFrac*min(L[x],L[y]): conf[x].add(y); conf[y].add(x)   # fraction of the shorter
lowinfo=set()
if a.stats:
    for l in open(a.stats):
        if l[0]=="#": continue
        r=l.split("\t")
        if float(r[2])<a.minWithin: lowinfo.add(r[0])
shared=set(l.split()[0] for l in open(a.shared) if l.strip() and l[0]!="#") if a.shared else set()
seeds=[c for c in L if L[c]>=a.seedLen and dep.get(c,0)>0 and c not in lowinfo]
print(f"# {len(seeds)} seeds; low-information chains >= seedLen held back: {sum(1 for c in lowinfo if L[c]>=a.seedLen)}, {sum(L[c] for c in lowinfo if L[c]>=a.seedLen)/1e6:.1f} Mb")
def dsum(S): return max(sum(dep.get(x,0) for x in S),1)
def avg(A,B): return sum(lk[x].get(y,0) for x in A for y in B)/(dsum(A)*dsum(B))  # average linkage: unbiased by size
def clash(A,B): return any(conf[x]&B for x in A)
def agglomerate(items,k,constrained):
    C={i:{c} for i,c in enumerate(items)}
    while len(C)>k:
        best=None; ks=list(C)
        for i,g in enumerate(ks):
            for h in ks[i+1:]:
                s=avg(C[g],C[h])
                if (best is None or s>best[0]) and not (constrained and clash(C[g],C[h])): best=(s,g,h)
        if best is None or best[0]==0: print(f"# stopped at {len(C)} clusters (wanted {k})"); break
        C[best[1]]|=C.pop(best[2])
    return C
G={}
if a.hapFirst:
    target=a.chroms*a.ploidy; hs=[c for c in sorted(seeds,key=lambda c:-L[c]) if c not in shared]
    print(f"# haplotype-first: {len(hs)} seeds, {sum(L[c] for c in hs)/1e6:.1f} Mb, target {target} groups")
    C={i:{c} for i,c in enumerate(hs)}; scores=[]; snap=None
    while len(C)>target-4:
        if len(C)==target: snap={k:set(v) for k,v in C.items()}
        best=None; ks=list(C)
        for i,g in enumerate(ks):
            for h in ks[i+1:]:
                x=avg(C[g],C[h])
                if (best is None or x>best[0]) and not clash(C[g],C[h]): best=(x,g,h)
        if best is None or best[0]==0: break
        scores.append((len(C),best[0]*1e9))
        C[best[1]]|=C.pop(best[2])
    print("# merge scores (clusters before merge, score):"," ".join(f"{k}:{v:.1f}" for k,v in scores if k<=target+4))
    if snap is None: snap={k:set(v) for k,v in C.items()}; print(f"# WARNING: clustering stopped at {len(C)} groups and never passed through {target}")
    # Gather groups into chromosomes, --ploidy at a time, by allelic bp. The sets are chosen on identity-weighted bp; plain bp is
    # what is reported, and what numbers the chromosomes.
    gl=list(snap); owner={c:g for g,S in snap.items() for c in S}; bp=collections.Counter(); wbp=collections.Counter()
    for l in open(a.alleles):
        if l[0]=="#": continue
        r=l.split("\t"); x,y=r[0],r[1]
        if x in sat or y in sat or x not in owner or y not in owner or owner[x]==owner[y]: continue
        k=frozenset((owner[x],owner[y])); n=min(int(r[2]),int(r[3])); bp[k]+=n; wbp[k]+=n*max(0.0,1-float(r[4]))**a.deWeight
    S={g:[g] for g in gl}; made=[]
    while True:
        best=None; ks=list(S)
        for i,g in enumerate(ks):
            for h in ks[i+1:]:
                if len(S[g])+len(S[h])>a.ploidy: continue
                x=sum(wbp[frozenset((u,v))] for u in S[g] for v in S[h])/(len(S[g])*len(S[h]))
                if x>0 and (best is None or x>best[0]): best=(x,g,h)
        if best is None: break
        S[best[1]]+=S.pop(best[2])
        if best[1] not in made: made.append(best[1])
        if best[2] in made: made.remove(best[2])
    inside=lambda M: [bp[frozenset((u,v))] for i,u in enumerate(M) for v in M[i+1:]]
    keys={k:k for k in bp}; cid=0; word="pair" if a.ploidy==2 else "set"
    for g in sorted(made,key=lambda g:-sum(inside(S[g]))):
        M=S[g]; M=list(keys.get(frozenset(M),M)) if len(M)==2 else M
        for j,x in enumerate(M): G[(cid,j)]=set(snap[x])
        out=max((w for q,w in bp.items() if len(q&set(M))==1),default=0)
        print(f"# {word} c{cid}: {' + '.join(f'{sum(L[c] for c in snap[x])/1e6:.1f}' for x in M)} Mb, allelic {min(inside(M))/1e6:.1f} Mb{'' if len(M)==2 else ' in its weakest pair'}; next best for {'either' if len(M)==2 else 'any member'} {out/1e6:.1f} Mb"+("" if len(M)==a.ploidy else f"; WARNING: {len(M)} of {a.ploidy} haplotypes"))
        cid+=1
    for g in gl:
        if g in S and g not in made: G[(cid,0)]=set(snap[g]); cid+=1; print(f"# WARNING: unpaired group of {sum(L[c] for c in snap[g])/1e6:.1f} Mb")
else:
    # Stage 1: chromosomes. Stage 2: haplotypes within each, never across an allelic conflict.
    chroms=agglomerate(seeds,a.chroms,False)
    print("# chromosome cluster sizes Mb:",sorted((round(sum(L[c] for c in S)/1e6,1) for S in chroms.values()),reverse=True))
    for cid,S in chroms.items():
        hs=sorted((c for c in S if c not in shared),key=lambda c:-L[c])
        print(f"# chromosome {cid}: {len(S)} seeds, {len(hs)} haplotype-specific, {sum(L[c] for c in S if c in shared)/1e6:.1f} Mb shared, {sum(L[c] for c in hs)/1e6:.1f} Mb specific")
        for j,H in enumerate(agglomerate(hs,a.ploidy,True).values()): G[(cid,j)]=set(H)
# Single-chain moves while same-haplotype links per depth improve. A chain may change chromosome, because
# whole-chain linkage in step 1 can be diluted for a chain that only touches its neighbours at one end
of={c:g for g,S in G.items() for c in S}
changed=True; rounds=0
while changed and rounds<50:
    changed=False; rounds+=1
    for c in sorted(of,key=lambda c:L[c]):
        cur=of[c]
        val=lambda g:sum(lk[c].get(x,0) for x in G[g] if x!=c)/max(sum(dep.get(x,0) for x in G[g] if x!=c),1)
        bg,bs=cur,val(cur)
        for g in G:
            if g==cur or conf[c]&G[g]: continue
            if val(g)>bs*1.05: bg,bs=g,val(g)
        if bg!=cur and len(G[cur])>1: G[cur].discard(c); G[bg].add(c); of[c]=bg; changed=True
print("# move rounds",rounds)
# shorter chains
unas=[]
for c in sorted((c for c in L if c not in of and c not in shared),key=lambda c:-L[c]):
    sc=[]
    for g,S in G.items():
        n=sum(lk[c].get(x,0) for x in S)
        sc.append((n/math.sqrt(max(sum(dep.get(x,0) for x in S),1)),n,g,bool(conf[c]&S)))
    sc.sort(reverse=True)
    allowed=[x for x in sc if not x[3]]
    if allowed and allowed[0][1]>=a.minLinks and (len(allowed)<2 or allowed[0][0]>=a.ratio*max(allowed[1][0],1e-12)):
        G[allowed[0][2]].add(c); of[c]=allowed[0][2]
    else: unas.append(c)
order=sorted(G,key=lambda g:-sum(L[c] for c in G[g]))
with open(a.out+".tsv","w") as o:
    o.write("#chain\tlength\tgroup\n")
    for i,g in enumerate(order,1):
        for c in sorted(G[g],key=lambda c:-L[c]): o.write(f"{c}\t{L[c]}\tg{i:02d}\n")
    for c in unas: o.write(f"{c}\t{L[c]}\tunassigned\n")
    for c in sorted(shared,key=lambda c:-L[c]): o.write(f"{c}\t{L[c]}\tshared\n")
rank={g:f"g{i:02d}" for i,g in enumerate(order,1)}
with open(a.out+".chroms.tsv","w") as o:
    for g in order: o.write(f"{rank[g]}\tc{g[0]}\n")
if shared:
    # Links per unit of chromosome depth. Assigned shared chains join their chromosome for the next round, because a
    # mostly collapsed chromosome has little haplotype-specific sequence to link to.
    cids=sorted({g[0] for g in G}); mem={k:set().union(*(G[g] for g in G if g[0]==k)) for k in cids}; call={}
    for rnd in range(5):
        new={}
        for c in shared:
            sc=[]
            for k in cids:
                X=mem[k]|{x for x,y in call.items() if y==k and x!=c}
                n=sum(lk[c].get(x,0) for x in X); sc.append((n/max(sum(dep.get(x,0) for x in X),1),n,k))
            sc.sort(reverse=True); b=sc[0]; r=b[0]/max(sc[1][0],1e-12)
            if b[1]>=a.minLinks and r>=a.ratio: new[c]=b[2]
        if new==call: break
        call=new
    print(f"# shared assignment rounds {rnd+1}")
    # Dosage: which haplotypes of its chromosome carry each assigned shared chain
    med={}; glink=collections.defaultdict(collections.Counter)
    if a.classDepth: med={l.split("\t")[0]:float(l.split("\t")[2]) for l in open(a.classDepth) if l[0]!="#"}
    if a.endLinks:
        for l in open(a.endLinks):
            if l[0]=="#": continue
            r=l.rstrip("\n").split("\t")
            if int(r[4])<a.minOverlap: continue
            u,v=r[0][:-1],r[1][:-1]
            for s,t in ((u,v),(v,u)):
                if s in shared and t in of: glink[s][of[t]]+=1
    carry={}; dos=[]
    for c in sorted(call,key=lambda c:-L[c]):
        gs=[g for g in G if g[0]==call[c]]; d=min(max(int(round(med[c]/a.hap)),2),len(gs)) if c in med else len(gs)
        if d>=len(gs): carry[c]=gs; continue
        R=sorted(((glink[c][g]>0,not conf[c]&G[g],sum(lk[c].get(x,0) for x in G[g])/dsum(G[g]),g) for g in gs),key=lambda r:r[:3],reverse=True); In,Out=R[:d],R[d:]
        ratio=min(r[2] for r in In)/max(max(r[2] for r in Out),1e-12)
        how="graph" if all(r[0] for r in In) and not any(r[0] for r in Out) else "alleles" if all(r[1] for r in In) and not any(r[1] for r in Out) else "hic" if ratio>=a.ratio else "arbitrary"
        carry[c]=[g for g in gs if g in {r[3] for r in In}]
        dos.append((c,L[c],f"c{call[c]}",f"{med[c]:.1f}",d,",".join(rank[g] for g in carry[c]),",".join(rank[r[3]] for r in Out),f"{min(ratio,999):.2f}",how))
    if dos:
        with open(a.out+".dosage.tsv","w") as o:
            o.write("#chain\tlength\tchrom\tmedian_depth\tdosage\tgroups\tleft_out\thic_ratio\tcall\n")
            for r in dos: o.write("\t".join(map(str,r))+"\n")
        print(f"# dosage: {len(dos)} shared chains in fewer haplotypes than their chromosome has, {sum(r[1] for r in dos)/1e6:.1f} Mb; calls {dict(collections.Counter(r[-1] for r in dos))}")
    with open(a.out+".shared.tsv","w") as o:
        o.write("#chain\tlength\tchrom\tlinks\tgroups\tscore_ratio\n")
        for c in sorted(shared,key=lambda c:-L[c]):
            sc=[]
            for k in cids:
                X=mem[k]|{x for x,y in call.items() if y==k and x!=c}
                n=sum(lk[c].get(x,0) for x in X); sc.append((n/max(sum(dep.get(x,0) for x in X),1),n,k))
            sc.sort(reverse=True); b=sc[0]; r=b[0]/max(sc[1][0],1e-12); ok=c in call
            gs=",".join(rank[g] for g in carry[c]) if ok else "."
            o.write(f"{c}\t{L[c]}\t{'c'+str(call[c]) if ok else 'unassigned'}\t{b[1]}\t{gs}\t{r:.2f}\n")
