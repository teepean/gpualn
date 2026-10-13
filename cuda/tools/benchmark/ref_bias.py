#!/usr/bin/env python3
"""Reference bias at heterozygous SNPs, measured against known truth.

At a het site the simulator drew each read's haplotype 50/50, so the EXPECTED alt-allele fraction
is exactly 0.500. Anything below that is reference bias: reads carrying the non-reference allele
mismatch the reference and are mapped less readily (or mapped elsewhere).

Two modes, and the comparison between them is the point:

  --truth <fq>   reads the simulated FASTQ and counts alleles at the sites each read TRULY spans,
                 using the encoded true position/strand. This is the perfect-mapper control: it
                 isolates bias contributed by post-mortem damage and sequencing error alone.
  --sam <sam>    counts alleles from the alignments, walking the CIGAR. The difference from the
                 control is the bias the MAPPER adds.

Transitions (A<->G, C<->T) are reported separately from transversions, because C>T / G>A damage
masquerades as a real allele at transition sites and cannot at transversion sites. Transversions
are therefore the damage-immune read on mapper bias.
"""
import argparse, sys
from collections import defaultdict

COMP = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
TRANSITION = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}

def load_het(path):
    het = {}
    for line in open(path):
        p, ref, alt, gt = line.rstrip('\n').split('\t')
        if gt in ('0|1', '1|0') and len(ref) == 1 and len(alt) == 1:
            het[int(p)] = (ref, alt)
    return het

def cigar_pairs(cig):
    n = 0
    for ch in cig:
        if ch.isdigit():
            n = n * 10 + int(ch)
        else:
            yield n, ch
            n = 0

def tally(counts, site_class, lbin, ref, alt, base):
    if base == ref:
        counts[(site_class, lbin)][0] += 1
    elif base == alt:
        counts[(site_class, lbin)][1] += 1
    else:
        counts[(site_class, lbin)][2] += 1

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--variants', required=True)
    ap.add_argument('--sam')
    ap.add_argument('--truth')
    ap.add_argument('--mapq', type=int, default=0)
    ap.add_argument('--trim', type=int, default=0, help='ignore this many bases at each read end')
    ap.add_argument('--label', default='')
    ap.add_argument('--dump', help='write per-site transversion ref/alt counts here (for bootstrap)')
    a = ap.parse_args()

    het = load_het(a.variants)
    counts = defaultdict(lambda: [0, 0, 0])      # (class, lbin) -> [ref, alt, other]
    persite = defaultdict(lambda: [0, 0])        # pos -> [ref, alt], transversions only

    if a.truth:
        fh = open(a.truth)
        for i, line in enumerate(fh):
            if i % 4 == 0:
                q = line[1:].rstrip('\n').split(':')
                pos, strand, L = int(q[2]), int(q[3]), int(q[4])
            elif i % 4 == 1:
                read = line.rstrip('\n')
                lbin = min(14, (L - 20) // 10)
                for off in range(L):
                    Q = pos + off
                    if Q not in het:
                        continue
                    roff = off if strand == 0 else L - 1 - off
                    if roff < a.trim or roff >= L - a.trim:
                        continue
                    b = read[roff] if strand == 0 else COMP.get(read[roff], 'N')
                    ref, alt = het[Q]
                    cls = 'ti' if (ref, alt) in TRANSITION else 'tv'
                    tally(counts, cls, lbin, ref, alt, b)
                    if cls == 'tv':
                        if b == ref: persite[Q][0] += 1
                        elif b == alt: persite[Q][1] += 1
    else:
        for line in open(a.sam):
            if line[0] == '@':
                continue
            f = line.split('\t', 11)
            flag = int(f[1])
            if flag & 0x904:                      # unmapped / secondary / supplementary
                continue
            if int(f[4]) < a.mapq:
                continue
            refpos, seq = int(f[3]), f[9]
            L = len(seq)
            lbin = min(14, (int(f[0].split(':')[4]) - 20) // 10)
            roff = 0
            for n, op in cigar_pairs(f[5]):
                if op in 'M=X':
                    for k in range(n):
                        Q = refpos + k
                        if Q in het:
                            ro = roff + k
                            if ro >= a.trim and ro < L - a.trim:
                                ref, alt = het[Q]
                                cls = 'ti' if (ref, alt) in TRANSITION else 'tv'
                                tally(counts, cls, lbin, ref, alt, seq[ro])
                                if cls == 'tv':
                                    if seq[ro] == ref: persite[Q][0] += 1
                                    elif seq[ro] == alt: persite[Q][1] += 1
                    refpos += n; roff += n
                elif op in 'D N':
                    refpos += n
                elif op in 'IS':
                    roff += n
                # H consumes neither

    def summarise(cls):
        r = sum(v[0] for (c, _), v in counts.items() if c == cls)
        al = sum(v[1] for (c, _), v in counts.items() if c == cls)
        o = sum(v[2] for (c, _), v in counts.items() if c == cls)
        return r, al, o

    print('%s   mapq>=%d trim=%d' % (a.label or (a.sam or a.truth), a.mapq, a.trim))
    print('  class          ref      alt    other   alt_frac   bias(0.5-alt_frac)')
    for cls, name in (('tv', 'transversion'), ('ti', 'transition'), (None, 'ALL')):
        if cls is None:
            r = sum(v[0] for v in counts.values()); al = sum(v[1] for v in counts.values())
            o = sum(v[2] for v in counts.values())
        else:
            r, al, o = summarise(cls)
        if r + al == 0:
            continue
        frac = al / (r + al)
        print('  %-12s %8d %8d %8d   %8.5f   %+8.5f' % (name, r, al, o, frac, 0.5 - frac))
    print('  --- transversion alt_frac by read length ---')
    for b in range(15):
        v = counts.get(('tv', b))
        if not v or v[0] + v[1] < 200:
            continue
        print('    %3d-%3d  n=%7d  alt_frac %8.5f  bias %+8.5f'
              % (20 + 10 * b, 29 + 10 * b, v[0] + v[1], v[1] / (v[0] + v[1]),
                 0.5 - v[1] / (v[0] + v[1])))
    if a.dump:
        with open(a.dump, 'w') as fh:
            for q in sorted(persite):
                fh.write('%d\t%d\t%d\n' % (q, persite[q][0], persite[q][1]))

def _unused():
    pass

if __name__ == '__main__':
    main()
