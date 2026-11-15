"""Precision and accuracy at every MAPQ threshold (unclipped start, tol 0), and the best accuracy each
mapper reaches at a given precision -- the fair comparison when mappers put MAPQ on different scales.
   mapq_curve.py SAM..."""
import sys, re
from collections import defaultdict
def lead(c):
    m = re.match(r'^(\d+)[SH]', c); return int(m.group(1)) if m else 0
for path in sys.argv[1:]:
    tot = 0; m = defaultdict(int); c = defaultdict(int)
    for line in open(path):
        if line[0] == '@': continue
        f = line.split('\t', 6); fl = int(f[1])
        if fl & 0x900: continue
        tot += 1
        if fl & 4: continue
        q = f[0].split(':'); mq = int(f[4])
        m[mq] += 1
        if f[2] == q[1] and int(f[3]) - lead(f[5]) == int(q[2]): c[mq] += 1
    print(path)
    cm = cc = 0; rows = []
    for mq in sorted(m, reverse=True):
        cm += m[mq]; cc += c[mq]; rows.append((mq, 100 * cc / cm, 100 * cc / tot, 100 * cm / tot))
    for mq, p, a, s in sorted(rows):
        if mq in (0, 1, 5, 10, 15, 20, 25, 30, 37, 40, 50, 60) or mq == rows[0][0]:
            print('  MAPQ>=%-3d precision %7.3f  accuracy %7.3f  kept %7.3f' % (mq, p, a, s))
    for target in (99.9, 99.95, 99.99):
        best = max((a for mq, p, a, s in rows if p >= target), default=0)
        print('  best accuracy with precision >= %.2f%%: %.3f' % (target, best))
