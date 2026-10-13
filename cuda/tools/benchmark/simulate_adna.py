#!/usr/bin/env python3
"""Simulate aDNA reads from chr22 with a known true position per read.

Follows the EXPERIMENTAL DESIGN of Oliva et al. 2021 (Briefings in Bioinformatics) rather than
their code, which does not run (the published pipeline has a tab/space SyntaxError and the scoring
script has an R parse error). Deliberate design choices, all flagged in the README:

  * read lengths: 15 equal-sized bins spanning 20-169 bp, as in their cutting_modern_reads_to_aDNA
  * damage: La Brana mapDamage profile (the same table they passed to gargammel deamSim),
    applied as C>T from the 5' end and G>A from the 3' end
  * sample variation: NA19471 (LWK), biallelic SNVs only. Indels are skipped ON PURPOSE -- they
    shift coordinates between haplotype and reference and so blur the very ground truth the whole
    benchmark rests on. Reference bias at SNPs, the metric of interest, is unaffected.
  * truth: the 1-based leftmost REFERENCE coordinate of the fragment, encoded in the read name.
    Directly comparable to SAM POS for an end-to-end alignment on either strand, which removes
    the hardcoded `-(250 - Read.Length)` correction their scoring script needs.
"""
import argparse, bisect, os, random, sys
# run under `python3 -I` (no CWD / user-site on sys.path, so a planted module cannot hijack an
# import) and add only this script's own directory back, for the sibling damage_model module.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from damage_model import load_rates

COMP = str.maketrans('ACGTN', 'TGCAN')

def load_fasta(path):
    return ''.join(l.strip() for l in open(path) if not l.startswith('>')).upper()

def load_variants(path):
    pos, alleles = [], []
    for line in open(path):
        p, ref, alt, gt = line.rstrip('\n').split('\t')
        if gt == '0|0':
            continue
        h = gt.split('|')
        pos.append(int(p))
        alleles.append((ref, alt, int(h[0]), int(h[1])))
    return pos, alleles

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--ref', required=True)
    ap.add_argument('--variants', required=True)
    ap.add_argument('--damage', required=True)
    ap.add_argument('--n', type=int, default=300000)
    ap.add_argument('--out', required=True)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--qual', type=int, default=37, help='baseline base quality')
    ap.add_argument('--no-damage', action='store_true')
    a = ap.parse_args()

    rnd = random.Random(a.seed)
    seq = load_fasta(a.ref)
    vpos, vall = load_variants(a.variants)
    rates = load_rates(a.damage)
    maxd = max(p for (_, p) in rates) if rates else 0
    err = 10.0 ** (-a.qual / 10.0)
    qstr_cache = {}

    bins = [(20 + 10 * i, 29 + 10 * i) for i in range(15)]      # 20-29 ... 160-169
    N = len(seq)
    out = open(a.out, 'w')
    made = 0; tries = 0
    while made < a.n:
        tries += 1
        if tries > a.n * 50:
            sys.exit('too many rejections; reference may be mostly N')
        lo, hi = bins[made % 15]
        L = rnd.randint(lo, hi)
        start = rnd.randrange(0, N - L)                          # 0-based
        frag = seq[start:start + L]
        if 'N' in frag:
            continue
        hap = rnd.randrange(2)
        # apply this haplotype's SNVs inside the fragment
        i = bisect.bisect_left(vpos, start + 1)
        if i < len(vpos) and vpos[i] <= start + L:
            frag = list(frag)
            while i < len(vpos) and vpos[i] <= start + L:
                ref, alt, h0, h1 = vall[i]
                if (h0 if hap == 0 else h1) == 1:
                    off = vpos[i] - (start + 1)
                    if frag[off] == ref:                         # guard against ref mismatch
                        frag[off] = alt
                i += 1
            frag = ''.join(frag)
        strand = rnd.randrange(2)                                # 0 = forward, 1 = reverse
        read = frag.translate(COMP)[::-1] if strand else frag
        # post-mortem damage, in sequencing orientation
        if not a.no_damage:
            b = list(read)
            for k in range(min(maxd, L)):
                if b[k] == 'C' and rnd.random() < rates.get(('5p', k + 1), 0.0):
                    b[k] = 'T'
            for k in range(min(maxd, L)):
                j = L - 1 - k
                if b[j] == 'G' and rnd.random() < rates.get(('3p', k + 1), 0.0):
                    b[j] = 'A'
            read = ''.join(b)
        # sequencing error
        b = list(read)
        for k in range(L):
            if rnd.random() < err:
                b[k] = rnd.choice([c for c in 'ACGT' if c != b[k]])
        read = ''.join(b)
        if L not in qstr_cache:
            qstr_cache[L] = chr(33 + a.qual) * L
        out.write('@r%d:22:%d:%d:%d:%d\n%s\n+\n%s\n'
                  % (made, start + 1, strand, L, hap, read, qstr_cache[L]))
        made += 1
    out.close()
    print('wrote %d reads to %s (%d position draws)' % (made, a.out, tries), file=sys.stderr)

if __name__ == '__main__':
    main()
