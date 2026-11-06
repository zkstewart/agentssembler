#!/usr/bin/env python3
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
--hapFirst replaces steps 1-2: all seeds are clustered straight to chroms x ploidy haplotype groups under the conflict
constraint, and groups are then paired into chromosomes by allele alignment (ploidy 2 only). Use it when homologs are
divergent: Hi-C between homologs can then be as weak as between chromosomes, so forcing the chromosome count first
splits one chromosome by haplotype and merges two others. The merge scores around the target are printed: a cliff
right after the target is independent support for it.
Usage: hic_group.py links.tsv alleles.pairs.tsv chains.fai out_prefix --chroms N --ploidy N [--stats chain_stats.tsv] [--shared shared.txt]
'''
import argparse, math, collections
p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("links"); p.add_argument("alleles"); p.add_argument("fai"); p.add_argument("out")
p.add_argument("--chroms",type=int,required=True); p.add_argument("--ploidy",type=int,required=True); p.add_argument("--shared"); p.add_argument("--hapFirst",action="store_true"); p.add_argument("--seedLen",type=int,default=1000000)
p.add_argument("--minConflict",type=int,default=50000); p.add_argument("--minConflictFrac",type=float,default=0.1)
p.add_argument("--stats"); p.add_argument("--minWithin",type=float,default=700)
p.add_argument("--repeat",help="chain_repeat.tsv: conflicts involving a chain above --maxRepeat are ignored, because satellite arrays align to each other without being alleles")
p.add_argument("--maxRepeat",type=float,default=0.5)
p.add_argument("--minLinks",type=int,default=10); p.add_argument("--ratio",type=float,default=2.0); a=p.parse_args()
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
    if min(bx,by)>=a.minConflict and min(bx/L[x],by/L[y])>=a.minConflictFrac: conf[x].add(y); conf[y].add(x)
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
    if a.ploidy!=2: raise SystemExit("--hapFirst pairs groups by alleles for ploidy 2 only")
    target=a.chroms*a.ploidy; hs=[c for c in sorted(seeds,key=lambda c:-L[c]) if c not in shared]
    print(f"# haplotype-first: {len(hs)} seeds, {sum(L[c] for c in hs)/1e6:.1f} Mb, target {target} groups")
    C={i:{c} for i,c in enumerate(hs)}; scores=[]
    while len(C)>target-4:
        best=None; ks=list(C)
        for i,g in enumerate(ks):
            for h in ks[i+1:]:
                x=avg(C[g],C[h])
                if (best is None or x>best[0]) and not clash(C[g],C[h]): best=(x,g,h)
        if best is None or best[0]==0: break
        scores.append((len(C),best[0]*1e9))
        if len(C)==target: snap={k:set(v) for k,v in C.items()}
        C[best[1]]|=C.pop(best[2])
    print("# merge scores (clusters before merge, score):"," ".join(f"{k}:{v:.1f}" for k,v in scores if k<=target+4))
    # pair groups into chromosomes by allelic bp
    gl=list(snap); owner={c:g for g,S in snap.items() for c in S}; bp=collections.Counter()
    for l in open(a.alleles):
        if l[0]=="#": continue
        r=l.split("\t"); x,y=r[0],r[1]
        if x in sat or y in sat or x not in owner or y not in owner or owner[x]==owner[y]: continue
        bp[frozenset((owner[x],owner[y]))]+=min(int(r[2]),int(r[3]))
    paired=set(); cid=0
    for k,v in bp.most_common():
        x,y=tuple(k)
        if x in paired or y in paired: continue
        paired|={x,y}; G[(cid,0)]=set(snap[x]); G[(cid,1)]=set(snap[y]); cid+=1
        print(f"# pair c{cid-1}: {sum(L[c] for c in snap[x])/1e6:.1f} + {sum(L[c] for c in snap[y])/1e6:.1f} Mb, allelic {v/1e6:.1f} Mb; next best for either {max((w for q,w in bp.items() if q!=k and (x in q or y in q)),default=0)/1e6:.1f} Mb")
    for g in gl:
        if g not in paired: G[(cid,0)]=set(snap[g]); cid+=1; print(f"# WARNING: unpaired group of {sum(L[c] for c in snap[g])/1e6:.1f} Mb")
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
    with open(a.out+".shared.tsv","w") as o:
        o.write("#chain\tlength\tchrom\tlinks\tgroups\tscore_ratio\n")
        for c in sorted(shared,key=lambda c:-L[c]):
            sc=[]
            for k in cids:
                X=mem[k]|{x for x,y in call.items() if y==k and x!=c}
                n=sum(lk[c].get(x,0) for x in X); sc.append((n/max(sum(dep.get(x,0) for x in X),1),n,k))
            sc.sort(reverse=True); b=sc[0]; r=b[0]/max(sc[1][0],1e-12); ok=c in call
            gs=",".join(rank[g] for g in G if g[0]==call[c]) if ok else "."
            o.write(f"{c}\t{L[c]}\t{'c'+str(call[c]) if ok else 'unassigned'}\t{b[1]}\t{gs}\t{r:.2f}\n")
