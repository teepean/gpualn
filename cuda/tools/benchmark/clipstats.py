"""Fraction of mapped primary records with soft clipping, per length bin group.   clipstats.py SAM..."""
import sys, re
for path in sys.argv[1:]:
    n = c5 = c3 = any_ = 0; short = [0, 0]
    for line in open(path):
        if line[0] == '@': continue
        f = line.split('\t', 6)
        fl = int(f[1])
        if fl & 0x904: continue
        cig = f[5]; n += 1
        lead = re.match(r'^\d+S', cig) is not None; trail = re.search(r'\d+S$', cig) is not None
        c5 += lead; c3 += trail; any_ += (lead or trail)
        if int(f[0].split(':')[4]) < 64:
            short[0] += 1; short[1] += (lead or trail)
    print('%-40s mapped %8d  clipped %6.3f%%  (leading %6.3f%%, trailing %6.3f%%)  <64bp clipped %6.3f%%'
          % (path.rsplit('/', 1)[-1], n, 100 * any_ / n, 100 * c5 / n, 100 * c3 / n, 100 * short[1] / max(1, short[0])))
