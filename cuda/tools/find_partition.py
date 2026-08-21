import itertools, time
from math import comb
G=6.2e9
def contiguous(pi):
    s=set()
    for p in pi:
        s.add(p)
        if max(s)-min(s)+1!=len(s): return False
    return True
def graded(P,K): return [v for v in itertools.product(range(K+1),repeat=P) if all(v[i]<=v[i+1] for i in range(P-1))]

def run(K,P,rlen=32,maxs=12):
    dists=[d for d in itertools.product(range(K+1),repeat=P) if sum(d)<=K]
    N=len(dists); full=(1<<N)-1
    Us=[u for u in graded(P,K) if u[-1]==K and u[0]<=1]
    Ls=[l for l in graded(P,K) if l[0]==0]
    cands=[]
    for pi in [p for p in itertools.permutations(range(P)) if contiguous(p)]:
        cum=[list(itertools.accumulate(d[p] for p in pi)) for d in dists]
        LE=[[0]*(K+1) for _ in range(P)]; GE=[[0]*(K+1) for _ in range(P)]
        for j in range(P):
            for v in range(K+1):
                a=b=0
                for i,c in enumerate(cum):
                    if c[j]<=v: a|=1<<i
                    if c[j]>=v: b|=1<<i
                LE[j][v]=a; GE[j][v]=b
        for U in Us:
            mu=full
            for j in range(P): mu&=LE[j][U[j]]
            if not mu: continue
            for L in Ls:
                if any(L[j]>U[j] for j in range(P)): continue
                m=mu
                for j in range(P):
                    if L[j]: m&=GE[j][L[j]]
                    if not m: break
                if m: cands.append((m,bin(m).count("1"),pi,L,U))
    # greedy DISJOINT cover (exact partition)
    cands.sort(key=lambda c:-c[1])
    chosen=[]; cov=0
    while cov!=full and len(chosen)<=maxs:
        best=None
        for c in cands:
            if c[0]&cov: continue
            if best is None or c[1]>best[1]: best=c
        if best is None: break
        chosen.append(best); cov|=best[0]
    return N,len(cands),(chosen if cov==full else None)

for K,P in ((5,6),(5,7),(6,7)):
    t0=time.time(); N,nc,r=run(K,P)
    print("K=%d P=%d: %d distributions, %d candidates, %.1fs" % (K,P,N,nc,time.time()-t0))
    if r is None: print("   no exact partition by greedy\n"); continue
    print("   PARTITION, %d searches:" % len(r))
    for m,n,pi,L,U in r:
        print('      {"%s"} L="%s" U="%s"  covers %d' % ("".join(map(str,pi)),"".join(map(str,L)),"".join(map(str,U)),n))
    print()
