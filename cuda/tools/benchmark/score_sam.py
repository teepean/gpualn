#!/usr/bin/env python3
"""Score a SAM against the truth encoded in the read name by simulate_adna.py.

QNAME = r<idx>:<chrom>:<true 1-based leftmost ref pos>:<strand>:<len>:<hap>

A read counts as CORRECT when it is mapped, lands on the true chromosome, and |start - true| <= tol,
where start is POS minus any leading soft/hard clip (the unclipped start). The simulator encodes the
fragment's leftmost reference coordinate, which equals POS for an end-to-end alignment on either
strand; local aligners (strobealign, mem without -L) can clip a damaged or mismatching end, which
moves POS without the read being misplaced. --raw-pos compares POS itself (the behaviour before
2026-11-15, under which a clipped read was always wrong).

Reports, overall and per length bin: mapped fraction (sensitivity), precision (correct/mapped) and
accuracy (correct/total), optionally after a MAPQ filter. Secondary alignments and supplementary
records are ignored; only the primary record of each read is scored.
"""
import argparse, sys
from collections import defaultdict

def lead_clip(cig):
    """Length of the leading S/H operations of a CIGAR string."""
    n = c = 0
    for ch in cig:
        if ch.isdigit():
            n = n * 10 + int(ch)
        elif ch in 'SH':
            c += n; n = 0
        else:
            break
    return c

def start_of(f, raw):
    """Unclipped 1-based leftmost reference coordinate of a split SAM record."""
    return int(f[3]) if raw else int(f[3]) - lead_clip(f[5])

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('sam')
    ap.add_argument('--tol', type=int, default=0, help='max |POS-true| still called correct')
    ap.add_argument('--mapq', type=int, default=0, help='discard alignments below this MAPQ')
    ap.add_argument('--label', default='')
    ap.add_argument('--tsv', help='also append a machine-readable row per length bin here')
    ap.add_argument('--raw-pos', action='store_true', help='compare POS, not the unclipped start')
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
        if f[2] == true_chr and abs(start_of(f, a.raw_pos) - true_pos) <= a.tol:
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
