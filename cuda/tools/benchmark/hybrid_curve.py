#!/usr/bin/env python3
"""Exact accuracy of the length-split hybrid at every cutoff, from one aln SAM + one mem SAM.

The pipeline sends each read to exactly one mapper based on its length, and single-end bwa aln and
bwa mem are both per-read independent, so the hybrid's output for cutoff c is exactly:
    read.len <  c  -> the aln record
    read.len >= c  -> the mem record
No separate run per cutoff is needed, and the curve is exact at 1 bp resolution.
"""
import argparse, sys

def load(path, tol, mapq):
    """-> {idx: (length, correct)} for primary records"""
    d = {}
    for line in open(path):
        if line[0] == '@':
            continue
        f = line.split('\t', 11)
        flag = int(f[1])
        if flag & 0x900:
            continue
        q = f[0].split(':')
        idx, true_chr, true_pos, L = int(q[0][1:]), q[1], int(q[2]), int(q[4])
        ok = 0
        if not (flag & 4) and int(f[4]) >= mapq:
            if f[2] == true_chr and abs(int(f[3]) - true_pos) <= tol:
                ok = 1
        d[idx] = (L, ok)
    return d

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--aln', required=True)
    ap.add_argument('--mem', required=True)
    ap.add_argument('--tol', type=int, default=0)
    ap.add_argument('--mapq', type=int, default=0)
    ap.add_argument('--cutoffs', default='0,30,40,50,55,60,61,64,66,70,75,80,90,100,170')
    a = ap.parse_args()

    A = load(a.aln, a.tol, a.mapq)
    M = load(a.mem, a.tol, a.mapq)
    common = A.keys() & M.keys()
    if len(common) != len(A) or len(common) != len(M):
        print('WARNING: aln has %d reads, mem has %d, %d in common'
              % (len(A), len(M), len(common)), file=sys.stderr)
    T = len(common)
    print('hybrid accuracy (correct/total), tol=%d mapq>=%d, %d reads' % (a.tol, a.mapq, T))
    print('  cutoff   accuracy%   vs all-aln   vs all-mem      note')
    all_aln = sum(A[i][1] for i in common) / T * 100
    all_mem = sum(M[i][1] for i in common) / T * 100
    for c in [int(x) for x in a.cutoffs.split(',')]:
        corr = sum((A[i][1] if A[i][0] < c else M[i][1]) for i in common)
        acc = corr / T * 100
        note = ''
        if c == 0:   note = '= all mem'
        if c >= 170: note = '= all aln'
        print('  %6d   %9.4f   %+10.4f   %+10.4f   %s'
              % (c, acc, acc - all_aln, acc - all_mem, note))
    print('  (all-aln %.4f   all-mem %.4f)' % (all_aln, all_mem))

if __name__ == '__main__':
    main()
