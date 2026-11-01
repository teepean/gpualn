#!/usr/bin/env python3
"""Split a FASTQ into the reads that have >=1 alignment in a .sai and the reads that have none.

    split_by_sai.py IN.sai IN.fq HIT.fq NOHIT.fq

.sai layout (bwa aln): 4-byte magic + 64-byte gap_opt_t header, then per read an int32 n_aln
followed by n_aln 24-byte bwt_aln1_t records, in input order.
"""
import struct
import sys


def main():
    sai, fq, out_hit, out_nohit = sys.argv[1:5]
    with open(sai, 'rb') as f:
        data = f.read()
    off, n_aln = 68, []
    while off < len(data):
        (n,) = struct.unpack_from('<i', data, off)
        n_aln.append(n)
        off += 4 + 24 * n
    if off != len(data):
        sys.exit('truncated .sai')
    n_hit = n_no = 0
    with open(fq) as fin, open(out_hit, 'w') as fh, open(out_nohit, 'w') as fn:
        for i, n in enumerate(n_aln):
            rec = [fin.readline() for _ in range(4)]
            if not rec[3]:
                sys.exit('FASTQ shorter than .sai at read %d' % i)
            if n > 0:
                fh.writelines(rec)
                n_hit += 1
            else:
                fn.writelines(rec)
                n_no += 1
        if fin.readline():
            sys.exit('FASTQ longer than .sai')
    print('%d reads: %d with hits (%.2f%%), %d without' % (len(n_aln), n_hit, 100.0 * n_hit / len(n_aln), n_no))


if __name__ == '__main__':
    main()
