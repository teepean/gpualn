#!/usr/bin/env python3
"""Paired bootstrap over het SITES for reference-bias differences between strategies.

All strategies see the same simulated reads, so the comparison is paired and the site is the
natural resampling unit (observations within a site are not independent). Resampling sites with
replacement and recomputing each strategy's bias gives a CI on the DIFFERENCE that accounts for
that pairing -- an unpaired SE (0.5/sqrt(n) ~ 0.0024 here) would be far too conservative.
"""
import random, sys
from collections import defaultdict

def load(path):
    d = {}
    for line in open(path):
        p, r, a = line.split()
        d[int(p)] = (int(r), int(a))
    return d

names = sys.argv[1::2]
paths = sys.argv[2::2]
tabs = [load(p) for p in paths]
sites = sorted(set().union(*[set(t) for t in tabs]))
B = 2000
rnd = random.Random(1)

def bias(tab, idx):
    r = a = 0
    for i in idx:
        s = sites[i]
        if s in tab:
            r += tab[s][0]; a += tab[s][1]
    return 0.5 - a / (r + a) if (r + a) else float('nan')

n = len(sites)
point = [bias(t, range(n)) for t in tabs]
boots = [[] for _ in tabs]
for _ in range(B):
    idx = [rnd.randrange(n) for _ in range(n)]
    for k, t in enumerate(tabs):
        boots[k].append(bias(t, idx))

print('paired bootstrap over %d het sites, B=%d' % (n, B))
for k, nm in enumerate(names):
    lo, hi = sorted(boots[k])[int(.025*B)], sorted(boots[k])[int(.975*B)]
    print('  %-10s bias %+0.5f   95%% CI [%+0.5f, %+0.5f]' % (nm, point[k], lo, hi))
print('  --- paired differences (negative = first has LESS bias) ---')
for i in range(len(tabs)):
    for j in range(i+1, len(tabs)):
        d = [boots[i][b] - boots[j][b] for b in range(B)]
        ds = sorted(d)
        lo, hi = ds[int(.025*B)], ds[int(.975*B)]
        pt = point[i] - point[j]
        sig = 'SIGNIFICANT' if (lo > 0 or hi < 0) else 'not significant'
        print('  %-10s - %-10s  %+0.5f   95%% CI [%+0.5f, %+0.5f]   %s'
              % (names[i], names[j], pt, lo, hi, sig))
