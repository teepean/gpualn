"""Length-split hybrid SAM: primary records of reads shorter than CUT from SHORT_SAM, the rest from LONG_SAM.
   merge_hybrid.py SHORT_SAM LONG_SAM CUT OUT_SAM"""
import sys
short, long_, cut, out = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
with open(out, 'w') as o:
    for path, keep in ((short, lambda L: L < cut), (long_, lambda L: L >= cut)):
        for line in open(path):
            if line[0] == '@':
                if path == short: o.write(line)
                continue
            f = line.split('\t', 2)
            if int(f[1]) & 0x900: continue
            if keep(int(f[0].split(':')[4])): o.write(line)
