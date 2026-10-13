#!/usr/bin/env python3
"""Score a SAM against the truth encoded in the read name by simulate_adna.py.

QNAME = r<idx>:<chrom>:<true 1-based leftmost ref pos>:<strand>:<len>:<hap>

A read counts as CORRECT when it is mapped, lands on the true chromosome, and |POS - true| <= tol.
POS is SAM's leftmost reference coordinate, which is what the simulator encodes, so no
strand-dependent or clipping-dependent correction is needed.

Reports, overall and per length bin: mapped fraction (sensitivity), precision (correct/mapped) and
accuracy (correct/total), optionally after a MAPQ filter. Secondary alignments and supplementary
records are ignored; only the primary record of each read is scored.
"""
import argparse, sys
from collections import defaultdict

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('sam')
    ap.add_argument('--tol', type=int, default=0, help='max |POS-true| still called correct')
    ap.add_argument('--mapq', type=int, default=0, help='discard alignments below this MAPQ')
    ap.add_argument('--label', default='')
    ap.add_argument('--tsv', help='also append a machine-readable row per length bin here')
    a = ap.parse_args()

    tot = defaultdict(int); mapped = defaultdict(int); corr = defaultdict(int)
    fh = sys.stdin if a.sam == '-' else open(a.sam)
    for line in fh:
        if line[0] == '@':
            continue
        f = line.split('\t', 11)
        flag = int(f[1])
        if flag & 0x900:                      # secondary / supplementary
            continue
        q = f[0].split(':')
        true_chr, true_pos, L = q[1], int(q[2]), int(q[3 + 1])
        b = min(14, (L - 20) // 10)
        tot[b] += 1
        if flag & 4:
            continue
        if int(f[4]) < a.mapq:
            continue
        mapped[b] += 1
        if f[2] == true_chr and abs(int(f[3]) - true_pos) <= a.tol:
            corr[b] += 1

    T = sum(tot.values()); M = sum(mapped.values()); C = sum(corr.values())
    print('%s  tol=%d mapq>=%d' % (a.label or a.sam, a.tol, a.mapq))
    print('  bin        reads   mapped%    precision%   accuracy%')
    for b in range(15):
        if not tot[b]:
            continue
        print('  %3d-%3d  %7d   %7.3f   %9.3f   %9.3f'
              % (20 + 10 * b, 29 + 10 * b, tot[b],
                 100.0 * mapped[b] / tot[b],
                 100.0 * corr[b] / mapped[b] if mapped[b] else 0.0,
                 100.0 * corr[b] / tot[b]))
    print('  ALL      %7d   %7.3f   %9.3f   %9.3f'
          % (T, 100.0 * M / T, 100.0 * C / M if M else 0.0, 100.0 * C / T))
    if a.tsv:
        with open(a.tsv, 'a') as out:
            for b in range(15):
                if not tot[b]:
                    continue
                out.write('%s\t%d\t%d\t%d\t%d\t%d\t%d\n'
                          % (a.label, a.tol, a.mapq, 20 + 10 * b, tot[b], mapped[b], corr[b]))

if __name__ == '__main__':
    main()
