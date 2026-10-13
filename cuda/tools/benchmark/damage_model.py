"""Position-dependent post-mortem damage rates from a mapDamage misincorporation table.

Source table: gargammel examplesMapDamage/results_LaBrana/misincorporation.txt, i.e. the same
damage profile Oliva et al. fed to `deamSim -mapdamage`. We read C>T per position from the 5' end
and G>A per position from the 3' end (the double-stranded-library signature that this table
actually measures).
"""
import sys
from collections import defaultdict

def load_rates(path, maxpos=25):
    num = defaultdict(float); den = defaultdict(float)
    with open(path) as fh:
        hdr = None
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if hdr is None:
                hdr = {name: i for i, name in enumerate(f)}
                continue
            end, pos = f[hdr['End']], int(f[hdr['Pos']])
            if pos > maxpos:
                continue
            if end == '5p':
                num[('5p', pos)] += float(f[hdr['C>T']]); den[('5p', pos)] += float(f[hdr['C']])
            elif end == '3p':
                num[('3p', pos)] += float(f[hdr['G>A']]); den[('3p', pos)] += float(f[hdr['G']])
    return {k: (num[k] / den[k] if den[k] else 0.0) for k in den}

if __name__ == '__main__':
    r = load_rates(sys.argv[1])
    print("pos    5' C>T     3' G>A")
    for p in range(1, 16):
        print("%3d   %7.4f   %7.4f" % (p, r.get(('5p', p), 0), r.get(('3p', p), 0)))
