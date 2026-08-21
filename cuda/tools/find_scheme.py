import itertools
from math import comb
G = 6.2e9

def contiguous(pi):
    seen=set()
    for p in pi:
        seen.add(p)
        if max(seen)-min(seen)+1 != len(seen): return False
    return True

def graded(P,K):
    return [v for v in itertools.product(range(K+1), repeat=P) if all(v[i]<=v[i+1] for i in range(P-1))]

def admits(pi,L,U,d):
    c=0
    for j,p in enumerate(pi):
        c += d[p]
        if c < L[j] or c > U[j]: return False
    return True

def scost(pi,L,U,P,rlen):
    t=0.0
    for j in range(P):
        ln = int(rlen*(j+1)/P)
        nb = sum(comb(ln,e)*(3**e) for e in range(U[j]+1))
        t += min(G, G*nb/(4.0**ln)) + nb
    return t

def solve(K,P,rlen,verbose=True):
    dists=[d for d in itertools.product(range(K+1),repeat=P) if sum(d)<=K]
    pis=[pi for pi in itertools.permutations(range(P)) if contiguous(pi)]
    gr=graded(P,K)
    cands=[]
    for pi in pis:
        for U in gr:
            if U[-1]!=K: continue
            for L in gr:
                if any(L[i]>U[i] for i in range(P)): continue
                if L[0]!=0: continue
                cov=frozenset(i for i,d in enumerate(dists) if admits(pi,L,U,d))
                if cov: cands.append((scost(pi,L,U,P,rlen), pi,L,U, cov))
    if verbose: print("  candidate searches: %d, distributions: %d" % (len(cands), len(dists)))
    # weighted greedy set cover
    need=set(range(len(dists))); chosen=[]
    while need:
        best=None
        for c,pi,L,U,cov in cands:
            g=len(cov & need)
            if g==0: continue
            r=c/g
            if best is None or r<best[0]: best=(r,c,pi,L,U,cov)
        if best is None: return None
        _,c,pi,L,U,cov = best
        chosen.append((c,pi,L,U)); need -= cov
    # local improvement: drop redundant searches
    changed=True
    while changed:
        changed=False
        for i in range(len(chosen)):
            rest=[chosen[j] for j in range(len(chosen)) if j!=i]
            covd=set()
            for _,pi,L,U in rest:
                covd |= {k for k,d in enumerate(dists) if admits(pi,L,U,d)}
            if len(covd)==len(dists):
                chosen=rest; changed=True; break
    return chosen

for (K,P) in ((3,4),(3,5),(4,5),(4,6)):
    for rlen in (32,):
        r=solve(K,P,rlen,verbose=(P==4))
        if r is None: print("K=%d P=%d: no cover"%(K,P)); continue
        tot=sum(c for c,_,_,_ in r)
        print("K=%d P=%d rlen=%d: %d searches, model cost %.3g" % (K,P,rlen,len(r),tot))
        for c,pi,L,U in r:
            print("     pi=%s L=%s U=%s   (%.2g)" % ("".join(map(str,pi)),"".join(map(str,L)),"".join(map(str,U)),c))
