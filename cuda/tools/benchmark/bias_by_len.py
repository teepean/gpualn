#!/usr/bin/env python3
"""Transversion-only reference bias per read-length bin, several SAMs, one pass each."""
import sys
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from collections import defaultdict
from ref_bias import load_het, cigar_pairs, TRANSITION

variants, mapq = sys.argv[1], int(sys.argv[2])
het = load_het(variants)
cols = []
for spec in sys.argv[3:]:
    path, lab = spec.split('=', 1)
    c = defaultdict(lambda: [0, 0])
    for line in open(path):
        if line[0] == '@':
            continue
        f = line.split('\t', 11)
        flag = int(f[1])
        if flag & 0x904 or int(f[4]) < mapq:
            continue
        refpos, seq = int(f[3]), f[9]
        lbin = min(14, (int(f[0].split(':')[4]) - 20) // 10)
        roff = 0
        for n, op in cigar_pairs(f[5]):
            if op in 'M=X':
                for k in range(n):
                    Q = refpos + k
                    if Q in het:
                        ref, alt = het[Q]
                        if (ref, alt) not in TRANSITION:
                            b = seq[roff + k]
                            if b == ref: c[lbin][0] += 1
                            elif b == alt: c[lbin][1] += 1
                refpos += n; roff += n
            elif op in 'DN': refpos += n
            elif op in 'IS': roff += n
    cols.append((lab, c))

print('transversion bias (0.5 - alt_frac), MAPQ>=%d; negative = anti-reference' % mapq)
print('  %-9s' % 'bin' + ''.join('%16s' % lab for lab, _ in cols))
for b in range(15):
    if all(sum(c[b]) < 150 for _, c in cols):
        continue
    row = '  %3d-%-5d' % (20 + 10 * b, 29 + 10 * b)
    for _, c in cols:
        r, a = c[b]
        row += '%16s' % ('%+0.5f (n=%d)' % (0.5 - a / (r + a), r + a) if r + a else '-')
    print(row)
