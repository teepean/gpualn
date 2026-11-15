#!/usr/bin/env python3
"""Write the strobealign-HRY benchmark results as a flat OpenDocument text file (.fodt) with tables.

   hry_report_odt.py WORKDIR SPEEDDIR OUT.fodt
   soffice --headless --convert-to odt OUT.fodt          # -> OUT.odt

WORKDIR holds the simulated-set outputs of hry_run.sh and the scoring steps (BENCHMARK_ACCURACY.md);
SPEEDDIR the outputs of hry_speed.sh. Everything is read from those files, except the MAPQ curves,
which are recomputed from the SAMs so that every threshold is present for every mapper.
"""
import os, re, sys
from xml.sax.saxutils import escape

W, P, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
STRATS = ['aln', 'mem', 'hry', 'hrydef', 'dk', 'strobe']
NAME = {'aln': 'bwa aln', 'mem': 'bwa-mem3 -L 15', 'hry': 'HRY (paper)', 'hrydef': 'HRY (--adna only)',
        'dk': 'strobealign-DK', 'strobe': 'strobealign default', 'hry50': 'HRY -r 50'}
BINS = list(range(20, 170, 10))

def rd(name, base=W):
    with open(os.path.join(base, name)) as f:
        return f.read()

def f3(x): return '%.3f' % x
def f5s(x): return '%+.5f' % x

# ---------------------------------------------------------------- data
def acc_tsv(n, q):
    d = {}
    for line in rd('acc_%s_%d.tsv' % (n, q)).splitlines():
        lab, tol, mq, b, tot, m, c = line.split('\t')
        d[int(b)] = (int(tot), int(m), int(c))
    return d

def all_line(text):
    m = re.search(r'^\s*ALL\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)', text, re.M)
    return [float(m.group(i)) for i in (2, 3, 4)]

def time_log(path):
    t = open(path).read()
    el = re.search(r'Elapsed \(wall clock\) time.*: (\S+)', t).group(1)
    s = 0.0
    for part in el.split(':'):
        s = s * 60 + float(part)
    rss = int(re.search(r'Maximum resident set size \(kbytes\): (\d+)', t).group(1))
    user = float(re.search(r'User time \(seconds\): ([\d.]+)', t).group(1))
    return s, rss / 1048576, user, t

def sa_times(text):
    out = []
    for key in ('Estimated read length', 'Read profile', 'Time reading reference', 'Total time reading index',
                'Total time indexing', 'Total time mapping', 'Total time creating strobemers',
                'Total time finding hits (non-rescue mode)', 'Total time finding hits (rescue mode)',
                'Total time chaining (non-rescue mode)', 'Total time chaining (rescue mode)',
                'Total time extending and pairing seeds', 'Total time writing alignments'):
        m = re.search(re.escape(key) + r': (.*)', text)
        if m:
            out.append((key, m.group(1).strip()))
    return out

def mapq_curve(path):
    """cumulative (mapped, correct) per MAPQ threshold, using the unclipped start"""
    tot = 0; m = {}; c = {}
    for line in open(path):
        if line[0] == '@':
            continue
        f = line.split('\t', 6)
        fl = int(f[1])
        if fl & 0x900:
            continue
        tot += 1
        if fl & 4:
            continue
        q = f[0].split(':'); mq = int(f[4])
        lead = re.match(r'^(\d+)[SH]', f[5])
        ok = f[2] == q[1] and int(f[3]) - (int(lead.group(1)) if lead else 0) == int(q[2])
        m[mq] = m.get(mq, 0) + 1
        c[mq] = c.get(mq, 0) + ok
    return tot, m, c

def cum_at(tot, m, c, thr):
    mm = sum(v for k, v in m.items() if k >= thr); cc = sum(v for k, v in c.items() if k >= thr)
    return (100 * cc / mm if mm else None), 100 * cc / tot, 100 * mm / tot

def best_at(tot, m, c, target):
    best = None
    for thr in sorted(m):
        p, a, k = cum_at(tot, m, c, thr)
        if p is not None and p >= target and (best is None or a > best):
            best = a
    return best

def rb_classes(text):
    d = {}
    for cls in ('transversion', 'transition', 'ALL'):
        m = re.search(r'^\s+%s\s+(\d+)\s+(\d+)\s+(\d+)\s+([\d.]+)\s+([+-][\d.]+)' % cls, text, re.M)
        if m:
            d[cls] = (int(m.group(1)), int(m.group(2)), int(m.group(3)), float(m.group(4)), float(m.group(5)))
    return d

# ---------------------------------------------------------------- ODF writer
body = []
ntab = [0]

def p(text, style='Body'):
    body.append('<text:p text:style-name="%s">%s</text:p>' % (style, escape(text)))

def h(text, level=1):
    if level == 0:                       # document title: a paragraph, not an outline heading
        return p(text, 'H0')
    body.append('<text:h text:style-name="H%d" text:outline-level="%d">%s</text:h>' % (level, level, escape(text)))

def bullets(items):
    for it in items:
        p('• ' + it, 'Bullet')

def table(caption, header, rows, widths=None, note=None, bold=None, numeric_from=1):
    """rows: list of lists of str; bold: set of (row, col) to emphasise"""
    ntab[0] += 1
    tn = 'Table%d' % ntab[0]
    p('Table %d. %s' % (ntab[0], caption), 'Caption')
    ncol = len(header)
    widths = widths or ([2] + [1] * (ncol - 1))
    tot = float(sum(widths))
    cols = []
    for i, w in enumerate(widths):
        cs = '%s.C%d' % (tn, i)
        autostyles.append('<style:style style:name="%s" style:family="table-column">'
                          '<style:table-column-properties style:rel-column-width="%d*"/></style:style>'
                          % (cs, int(1000 * w / tot)))
        cols.append('<table:table-column table:style-name="%s"/>' % cs)
    out = ['<table:table table:name="%s" table:style-name="Tbl">' % tn] + cols
    out.append('<table:table-header-rows><table:table-row>')
    for x in header:
        out.append('<table:table-cell table:style-name="HdrCell" office:value-type="string">'
                   '<text:p text:style-name="TH">%s</text:p></table:table-cell>' % escape(x))
    out.append('</table:table-row></table:table-header-rows>')
    for ri, r in enumerate(rows):
        out.append('<table:table-row>')
        for ci, x in enumerate(r):
            st = 'TDnum' if ci >= numeric_from else 'TD'
            if bold and (ri, ci) in bold:
                st += 'B'
            cell = 'Cell' if ri % 2 == 0 else 'CellAlt'
            out.append('<table:table-cell table:style-name="%s" office:value-type="string">'
                       '<text:p text:style-name="%s">%s</text:p></table:table-cell>' % (cell, st, escape(str(x))))
        out.append('</table:table-row>')
    out.append('</table:table>')
    body.append(''.join(out))
    if note:
        p(note, 'Note')

autostyles = []

# ---------------------------------------------------------------- content
h('strobealign-HRY vs bwa aln vs bwa-mem3: benchmark results', 0)
p('2026-11-15, gpualn project. Summary and discussion: cuda/BENCHMARK_ACCURACY.md, section '
  '"strobealign-HRY: tested 2026-11-15, not adopted". Plain-text version: cuda/HRY_BENCHMARK_RESULTS.txt.',
  'Subtitle')

h('Verdict')
p('Not adopted for the human aDNA pipeline. On ground truth (1.5 M simulated reads, hs37d5):')
bullets([
    'accuracy 92.744% vs bwa aln 93.334%; HRY is below aln in 14 of 15 length bins;',
    'at MAPQ≥25: 80.075% vs 86.731%; at equal precision (99.9%): 81.769% vs 89.836%;',
    'transversion reference bias (MAPQ≥25) +0.01819 vs aln +0.00952 (significant); the paper does not measure bias;',
    'either arm of the length-split hybrid is worse than production aln + bwa-mem3 (93.347%, bias +0.00913).',
    'Speed on real libraries (32 threads, saved index): 6-10x faster than CPU bwa aln, slower than gpualn on the '
    'two low-endogenous libraries, faster on the 94%-endogenous one; 21 GB RAM, 21 GB index.',
    'On real libraries HRY reports nearly every exogenous read as mapped (soft-clipped local hits) unless '
    'filtered on alignment score ≥ 60.',
])

h('The paper')
p('Tolstoganov I, Sahlin K. Hybrid RY-mers increase aDNA alignment sensitivity. bioRxiv '
  '10.64898/2026.11.05.764514, v2 posted 2026-11-10. Code: github.com/Itolstoganov/strobealign (tag hry); '
  'evaluation: github.com/Itolstoganov/adna_bench.')
table('How the paper\'s benchmark differs from this one', ['', 'the paper', 'here'], [
    ['reads', '40 bp only (gargammel), 80% endogenous + 20% bacterial', '20-169 bp in 15 bins, all endogenous; plus 3 real libraries'],
    ['references / divergence', 'fruit fly, ZymoBIOMICS, human CHM13 with 0.5 / 2 / 5% synthetic substitutions + indels (Mason)',
     'hs37d5 chr22 with one person\'s real SNVs (NA19471), ~0.1%'],
    ['"correct"', 'primary overlaps the true interval AND alignment score ≥ 60 (strobealign scoring; bwa rescored from MD/CIGAR)',
     'unclipped start equals the true start (tolerance 0)'],
    ['MAPQ-filtered results', 'none', 'MAPQ ≥ 0 and ≥ 25, plus precision-matched'],
    ['reference bias', 'not measured (no heterozygous sites in the design)', 'transversion bias at ~10,800 covered het sites, paired bootstrap'],
], widths=[1.3, 3, 3], numeric_from=99)

h('1. Setup')
table('Machine and tools', ['item', 'value'], [
    ['CPU', 'AMD Ryzen 9 9950X, 16 cores / 32 threads'],
    ['RAM / OS', '123 GB / Linux 7.3.1-1-cachyos'],
    ['GPU', 'not used (reserved for another process); gpualn / G-ALN rows are from 2026-11-01'],
    ['bwa', '0.7.18-r1243 (/usr/bin/bwa)'],
    ['bwa-mem3', 'v0.14.1, the pipeline\'s long arm'],
    ['strobealign', '0.18.0-alpha, tag hry = commit 4be2dcc, cargo build --release (rustc 1.99.1)'],
    ['reference', 'hs37d5'],
], widths=[1, 4], numeric_from=99)
table('Strategies and exact options', ['label', 'command / options', 'note'], [
    ['bwa aln', 'bwa aln -l 1024 -n 0.01 -o 2, then bwa samse', 'BWA1; also the paper\'s bwa aln options'],
    ['bwa-mem3 -L 15', 'bwa-mem3 mem -L 15 --supp-rep-hard-cap 5', 'pipeline MEM_OPTS since 2026-10-14'],
    ['HRY (paper)', 'strobealign --adna --ry-len 8 -M 500 --ssw -k 17 -s 15', 'adna_bench config: versions[0], main_k 17, s = k-3+1'],
    ['HRY -r 50', 'as HRY (paper) plus -r 50', 'byte-identical records to HRY (paper)'],
    ['HRY (--adna only)', 'strobealign --adna', 'release default'],
    ['strobealign-DK', 'strobealign -l 1 -u 1 -M 500 -k 17 -s 15', 'the paper\'s non-aDNA baseline'],
    ['strobealign default', 'strobealign', 'plain defaults, same binary'],
    ['hybrid X + Y', 'reads < 64 bp to X, ≥ 64 bp to Y', 'hybrid_curve.py / merge_hybrid.py'],
], widths=[1.3, 3, 2.2], numeric_from=99)
table('Datasets', ['dataset', 'reads', 'lengths', 'mapped by bwa aln', 'description'], [
    ['simulated (sim15M)', '1,500,000', '20-169 bp, 15 equal bins', '99.96%',
     'chr22 of hs37d5 with NA19471 SNVs (38,599 het sites), La Brana damage, Q37, seed 7; all endogenous'],
    ['sub100k', '100,000', '30-91 bp', '0.51%', 'head of AVA1B (soil)'],
    ['p17_shuf1m', '999,552', '30-49 bp', '94.42%', 'pipeline short branch, shuffled 1 M subsample'],
    ['p1_shuf1m', '1,000,528', '30-49 bp', '0.20%', 'pipeline short branch, shuffled 1 M subsample'],
], widths=[1.4, 1, 1.3, 1.1, 3.5], numeric_from=99)
table('Scoring definitions', ['term', 'definition'], [
    ['correct', 'primary record on the true chromosome with unclipped start (POS minus leading S/H) = true start'],
    ['mapped %', 'mapped / total'], ['precision %', 'correct / mapped'], ['accuracy %', 'correct / total'],
    ['bias', '0.5 - alt/(ref+alt) at heterozygous transversion sites, reads mapped to chr22 only'],
    ['control', 'alleles counted at the true positions (damage and sequencing error only)'],
    ['bootstrap', 'paired, het site resampled, B = 2000'],
], widths=[1, 5], numeric_from=99)

h('2. Simulated set: headline')
rows = []; accs = {}
for n in STRATS:
    a0 = all_line(rd('score_%s_0.txt' % n)); a25 = all_line(rd('score_%s_25.txt' % n)); raw = all_line(rd('scoreraw_%s_0.txt' % n))
    accs[n] = (a0[2], a25[2])
    rows.append([NAME[n]] + [f3(x) for x in a0] + [f3(x) for x in a25] + [f3(raw[2])])
bold = {(STRATS.index(max(STRATS, key=lambda n: accs[n][0])), 3), (STRATS.index(max(STRATS, key=lambda n: accs[n][1])), 6)}
table('All 1.5 M reads (bold: best accuracy)', ['strategy', 'mapped %', 'precision %', 'accuracy %', 'MAPQ≥25 mapped %',
      'MAPQ≥25 precision %', 'MAPQ≥25 accuracy %', 'accuracy, old raw-POS rule'], rows,
      widths=[1.8, 1, 1, 1, 1, 1, 1, 1.1], bold=bold,
      note='HRY -r 50 gave byte-identical records to HRY (paper). The raw-POS column is the scoring rule before '
           '2026-11-15, under which every read clipped at its left end counted as wrong.')

h('3. Accuracy by read length')
for q in (0, 25):
    d = {n: acc_tsv(n, q) for n in STRATS}
    rows = []; bold = set()
    for ri, b in enumerate(BINS):
        vals = [100 * d[n][b][2] / d[n][b][0] for n in STRATS]
        best = max(range(len(vals)), key=lambda i: vals[i]); bold.add((ri, best + 1))
        rows.append(['%d-%d' % (b, b + 9)] + [f3(v) for v in vals])
    table('Accuracy %% by length bin, MAPQ ≥ %d (bold: best in the row; 100,000 reads per bin)' % q,
          ['bin (bp)'] + [NAME[n] for n in STRATS], rows, widths=[1] + [1] * 6, bold=bold)
for q in (0, 25):
    d = {n: acc_tsv(n, q) for n in STRATS}
    rows = []
    for b in BINS:
        base = 100 * d['aln'][b][2] / d['aln'][b][0]
        rows.append(['%d-%d' % (b, b + 9), f3(base)] + ['%+.3f' % (100 * d[n][b][2] / d[n][b][0] - base) for n in STRATS[1:]])
    table('Accuracy difference to bwa aln (percentage points), MAPQ ≥ %d' % q,
          ['bin (bp)', 'bwa aln accuracy %'] + [NAME[n] for n in STRATS[1:]], rows, widths=[1, 1.1] + [1] * 5)

h('4. Per-strategy detail by length bin')
for n in STRATS:
    d0, d25 = acc_tsv(n, 0), acc_tsv(n, 25)
    rows = []
    for b in BINS + ['ALL']:
        if b == 'ALL':
            t0 = [sum(d0[x][i] for x in BINS) for i in range(3)]; t25 = [sum(d25[x][i] for x in BINS) for i in range(3)]
            lab = 'ALL'
        else:
            t0, t25 = d0[b], d25[b]; lab = '%d-%d' % (b, b + 9)
        def trio(t):
            tot, m, c = t
            return [f3(100 * m / tot), f3(100 * c / m) if m else '-', f3(100 * c / tot)]
        rows.append([lab, '{:,}'.format(t0[0])] + trio(t0) + trio(t25))
    table('%s: mapped / precision / accuracy by length bin' % NAME[n],
          ['bin (bp)', 'reads', 'mapped %', 'precision %', 'accuracy %', 'MAPQ≥25 mapped %', 'MAPQ≥25 precision %',
           'MAPQ≥25 accuracy %'], rows, widths=[1, 1.1, 1, 1, 1, 1, 1, 1])

h('5. Precision and accuracy at every MAPQ threshold')
p('MAPQ is on different scales in bwa and strobealign, so a fixed MAPQ cut is not a like-for-like comparison. '
  'These tables give precision and accuracy at each threshold, and the best accuracy each mapper reaches at a '
  'given precision. Recomputed from the SAMs (unclipped start, tolerance 0).')
curves = {n: mapq_curve(os.path.join(W, 'big.%s.sam' % n)) for n in STRATS}
THR = [0, 1, 5, 10, 15, 20, 25, 30, 35, 37, 40, 50, 60]
for what, idx in (('precision %', 0), ('accuracy %', 1), ('reads kept %', 2)):
    rows = []
    for t in THR:
        r = ['≥ %d' % t]
        for n in STRATS:
            tot, m, c = curves[n]
            if t > max(m):
                r.append('-')
            else:
                v = cum_at(tot, m, c, t)[idx]
                r.append(f3(v) if v is not None else '-')
        rows.append(r)
    table('%s by MAPQ threshold' % what.capitalize(), ['MAPQ'] + [NAME[n] for n in STRATS], rows, widths=[0.8] + [1] * 6,
          note='"-": above the mapper\'s highest MAPQ (bwa aln max 37, the others 60).' if idx == 0 else None)
rows = []
for target in (99.0, 99.5, 99.9, 99.95, 99.99):
    r = ['≥ %.2f%%' % target]
    for n in STRATS:
        v = best_at(*curves[n], target)
        r.append(f3(v) if v is not None else 'never')
    rows.append(r)
table('Best accuracy % reached at a given precision (any MAPQ threshold)', ['precision'] + [NAME[n] for n in STRATS], rows,
      widths=[0.9] + [1] * 6)

h('6. Checks')
rows = []
for n in ['aln', 'mem', 'hry', 'dk', 'strobe']:
    rows.append([NAME[n], f3(all_line(rd('score_%s_0.txt' % n))[2]), f3(all_line(rd('tol100_%s.txt' % n))[2])])
table('Position tolerance: is HRY\'s deficit near misses? (accuracy %, MAPQ ≥ 0)',
      ['strategy', 'tolerance 0 bp', 'tolerance 100 bp'], rows, widths=[2, 1, 1],
      note='A 100 bp tolerance (close to the paper\'s overlap rule) moves every mapper by 0.03-0.04 points.')
rows = []
fn = {'big.aln.sam': 'aln', 'big.mem.sam': 'mem', 'big.hry.sam': 'hry', 'big.hrydef.sam': 'hrydef', 'big.dk.sam': 'dk', 'big.strobe.sam': 'strobe'}
for line in rd('clipstats.txt').splitlines():
    m = re.match(r'(\S+)\s+mapped\s+(\d+)\s+clipped\s+([\d.]+)%\s+\(leading\s+([\d.]+)%, trailing\s+([\d.]+)%\)\s+<64bp clipped\s+([\d.]+)%', line)
    if m:
        rows.append([NAME[fn[m.group(1)]], '{:,}'.format(int(m.group(2))), m.group(3), m.group(4), m.group(5), m.group(6)])
table('Soft clipping among mapped primary records (%)', ['strategy', 'mapped records', 'clipped', 'leading clip',
      'trailing clip', 'clipped, reads < 64 bp'], rows, widths=[1.8, 1.2, 1, 1, 1, 1.2])
rows = []
for line in rd('hry_as60.txt').splitlines():
    m = re.match(r'\s*(\d+)-\s*(\d+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s+\|\s+([\d.]+)\s+([\d.]+)', line)
    if m:
        b = int(m.group(1)); a = acc_tsv('aln', 0)[b]
        rows.append(['%s-%s' % (m.group(1), m.group(2)), m.group(4), m.group(5), m.group(6), m.group(7), f3(100 * a[2] / a[0])])
table('HRY with the paper\'s alignment-score filter (AS ≥ 60), MAPQ ≥ 0', ['bin (bp)', 'mapped %', 'accuracy %',
      'AS≥60 mapped %', 'AS≥60 accuracy %', 'bwa aln accuracy %'], rows, widths=[1, 1, 1, 1, 1, 1],
      note='The maximum score for an L bp read is 2L + 20, so a 20 bp read passes AS ≥ 60 only as a perfect end-to-end match.')

h('7. Length-split hybrids')
hc = {}
for f in ('aln_mem', 'hry_mem', 'aln_hry'):
    for q in (0, 25):
        for line in rd('hc_%s_%d.txt' % (f, q)).splitlines():
            m = re.match(r'\s+(\d+)\s+([\d.]+)\s+[+-]', line)
            if m:
                hc[(f, q, int(m.group(1)))] = float(m.group(2))
rows = []
for c in (0, 40, 64, 100, 170):
    lab = {0: '0 (all long arm)', 170: '170 (all short arm)'}.get(c, str(c))
    rows.append([lab] + ['%.3f' % hc[(f, q, c)] for f in ('aln_mem', 'hry_mem', 'aln_hry') for q in (0, 25)])
bold = {(2, 1), (2, 2)}
table('Hybrid accuracy % by cutoff (short arm below the cutoff, long arm at or above)',
      ['cutoff (bp)', 'aln + mem3', 'aln + mem3, MAPQ≥25', 'HRY + mem3', 'HRY + mem3, MAPQ≥25', 'aln + HRY',
       'aln + HRY, MAPQ≥25'], rows, widths=[1.3, 1, 1, 1, 1, 1, 1], bold=bold,
      note='Production is aln + bwa-mem3 at cutoff 64 (bold).')

h('8. Reference bias at heterozygous sites')
ctrl = rb_classes(rd('rb_control.txt'))
p('Expected alt-allele fraction 0.500. Control (perfect mapper, alleles at the true positions): transversions '
  'ref %d / alt %d, alt fraction %.5f, bias %+.5f; transitions bias %+.5f.' % (
      ctrl['transversion'][0], ctrl['transversion'][1], ctrl['transversion'][3], ctrl['transversion'][4], ctrl['transition'][4]))
SB = ['aln', 'mem', 'hry', 'dk', 'strobe', 'h_aln_mem', 'h_hry_mem', 'h_aln_hry']
SBN = dict(NAME, h_aln_mem='aln + mem3 (production)', h_hry_mem='HRY + mem3', h_aln_hry='aln + HRY')
rows = []
for n in SB:
    r0 = rb_classes(rd('rb0_%s.txt' % n)); r25 = rb_classes(rd('rb25_%s.txt' % n))
    rows.append([SBN[n], f5s(r0['transversion'][4]), f5s(r25['transversion'][4]),
                 '{:,}'.format(r25['transversion'][0] + r25['transversion'][1]), f5s(r25['transition'][4]), f5s(r25['ALL'][4])])
table('Reference bias by strategy (0.5 - alt fraction; positive = toward the reference)',
      ['strategy', 'transversions, MAPQ≥0', 'transversions, MAPQ≥25', 'tv alleles counted, MAPQ≥25',
       'transitions, MAPQ≥25', 'all sites, MAPQ≥25'], rows, widths=[1.8, 1, 1, 1.1, 1, 1],
      note='Transversions are the damage-immune measure: C>T/G>A damage can fake an allele at a transition site, not at a transversion.')
rows = []
for n in SB:
    r25 = rb_classes(rd('rb25_%s.txt' % n))
    for cls in ('transversion', 'transition'):
        r, a, o, fr, b = r25[cls]
        rows.append([SBN[n], cls, '{:,}'.format(r), '{:,}'.format(a), '{:,}'.format(o), '%.5f' % fr, f5s(b)])
table('Allele counts at het sites, MAPQ ≥ 25', ['strategy', 'site class', 'ref', 'alt', 'other', 'alt fraction', 'bias'],
      rows, widths=[1.8, 1.1, 1, 1, 0.8, 1, 1], numeric_from=2)
rows = []
for n in ('aln', 'mem', 'h_aln_mem'):
    rows.append([SBN[n], f5s(rb_classes(rd('rbold25_%s.txt' % n))['transversion'][4]), f5s(rb_classes(rd('rb25_%s.txt' % n))['transversion'][4])])
table('Chromosome fix: old rule (reads on any chromosome counted) vs new, transversions, MAPQ ≥ 25',
      ['strategy', 'old rule', 'new rule (chr22 only)'], rows, widths=[2, 1, 1])
cols = ['aln', 'mem', 'hry', 'dk']; rows = []
for line in rd('bbl25.txt').splitlines():
    m = re.match(r'\s+(\d+)-(\d+)\s*(.*)$', line)
    if not m:
        continue
    toks = re.findall(r'[+-]\d\.\d{5} \(n=\d+\)|-', m.group(3))
    r = ['%s-%s' % (m.group(1), m.group(2))]
    for t in toks:
        if t == '-':
            r += ['unmapped', '-']
        else:
            v, nn = re.match(r'([+-]\d\.\d{5}) \(n=(\d+)\)', t).groups(); r += [v, nn]
    rows.append(r)
table('Transversion bias by read length, MAPQ ≥ 25 (n = alleles counted)',
      ['bin (bp)', 'bwa aln', 'n', 'bwa-mem3', 'n', 'HRY', 'n', 'DK', 'n'], rows, widths=[1, 1, 0.7, 1, 0.7, 1, 0.7, 1, 0.7],
      note='Negative values mean anti-reference; at these counts single bins carry roughly ±0.01-0.03 of noise.')
bs = rd('bootstrap25.txt'); rows = []; prs = []
for line in bs.splitlines():
    m = re.match(r'^\s+(\S+)\s+bias\s+([+-][\d.]+)\s+95% CI \[([+-][\d.]+), ([+-][\d.]+)\]', line)
    if m:
        rows.append([m.group(1), m.group(2), '[%s, %s]' % (m.group(3), m.group(4))])
    m = re.match(r'^\s+(\S+)\s+-\s+(\S+)\s+([+-][\d.]+)\s+95% CI \[([+-][\d.]+), ([+-][\d.]+)\]\s+(.*)$', line)
    if m:
        prs.append([m.group(1) + ' - ' + m.group(2), m.group(3), '[%s, %s]' % (m.group(4), m.group(5)), m.group(6).strip()])
nsite = re.search(r'over (\d+) het sites', bs).group(1)
table('Transversion bias with 95%% bootstrap CI, MAPQ ≥ 25 (%s het sites, B = 2000)' % nsite,
      ['strategy', 'bias', '95% CI'], rows, widths=[2, 1, 1.6])
table('Paired bootstrap differences (negative = the first has less bias)', ['comparison', 'difference', '95% CI', 'result'],
      prs, widths=[2.4, 1, 1.6, 1.2])

h('9. Simulated set: run time and memory (16 threads; not a speed benchmark)')
rows = []
for n in ['aln', 'mem', 'hry', 'hry50', 'hrydef', 'dk', 'strobe']:
    s, gb, user, _ = time_log(os.path.join(W, 'big.%s.log' % n))
    rows.append([NAME[n], '%.1f' % s, '%.1f' % user, '%.1f' % gb])
table('Wall time, CPU time and peak memory', ['strategy', 'wall (s)', 'user CPU (s)', 'peak RSS (GB)'], rows, widths=[2, 1, 1, 1],
      note='strobealign rows build the hs37d5 index inside the run; bwa-mem3 includes a cold index load (warm: about 3 s). '
           'bwa aln excludes samse.')
rows = []
keys = None; per = {}
for n in ('hry', 'dk', 'strobe'):
    per[n] = dict(sa_times(rd('big.%s.log' % n)))
    keys = keys or [k for k, _ in sa_times(rd('big.%s.log' % n))]
for k in keys:
    rows.append([k] + [per[n].get(k, '-') for n in ('hry', 'dk', 'strobe')])
table('strobealign time breakdown on the simulated set (strobealign -v)', ['step', NAME['hry'], NAME['dk'], NAME['strobe']],
      rows, widths=[2.4, 1.2, 1.2, 1.2],
      note='Per-step "Total time" lines are summed over all threads; "Total time indexing" and "Total time mapping" are wall-clock.')

h('10. Real libraries: speed')
p('32 threads, strictly one run at a time with nothing else running, page cache warmed before each tool. '
  'HRY uses an index built once (-r 50, --create-index) and loaded with --use-index; its output is byte-identical to '
  'the run that builds the index itself.')
LIBS = ['sub100k', 'p17_shuf1m', 'p1_shuf1m']
sp = {}
for lib in LIBS:
    sp[(lib, 'aln')] = time_log(os.path.join(P, '%s.aln.log' % lib))
    for i in (1, 2, 3):
        sp[(lib, 'hry%d' % i)] = time_log(os.path.join(P, '%s.hry.%d.log' % (lib, i)))
    sp[(lib, 'fly')] = time_log(os.path.join(P, '%s.hryfly.log' % lib))
def med(lib):
    return sorted(sp[(lib, 'hry%d' % i)][0] for i in (1, 2, 3))[1]
GPU = {'sub100k': (2.22, 1.97), 'p17_shuf1m': (17.4, 9.63), 'p1_shuf1m': (7.77, 6.95)}
OLD = {'sub100k': 15.3, 'p17_shuf1m': 56.1, 'p1_shuf1m': 149.1}
rows = [
    ['bwa aln -t 32 (today)'] + ['%.2f' % sp[(l, 'aln')][0] for l in LIBS],
    ['bwa aln -t 32 (2026-11-01)'] + ['%.1f' % OLD[l] for l in LIBS],
    ['strobealign-HRY, saved index (median of 3)'] + ['%.2f' % med(l) for l in LIBS],
    ['strobealign-HRY, index built in the run'] + ['%.2f' % sp[(l, 'fly')][0] for l in LIBS],
    ['gpualn, GPU (2026-11-01)'] + ['%.2f' % GPU[l][0] for l in LIBS],
    ['G-ALN 0.2.4 + scheme prefilter (2026-11-01)'] + ['%.2f' % GPU[l][1] for l in LIBS],
    ['HRY speed-up over CPU bwa aln'] + ['%.1fx' % (sp[(l, 'aln')][0] / med(l)) for l in LIBS],
    ['HRY vs gpualn'] + [('%.1fx faster' % (GPU[l][0] / med(l))) if med(l) < GPU[l][0] else ('%.1fx slower' % (med(l) / GPU[l][0])) for l in LIBS],
    ['peak RSS bwa aln / HRY (GB)'] + ['%.1f / %.1f' % (sp[(l, 'aln')][1], sp[(l, 'hry3')][1]) for l in LIBS],
]
table('Wall time in seconds', ['', 'sub100k (100 k)', 'p17 (1 M, 94% mapped)', 'p1 (1 M, 0.2% mapped)'], rows, widths=[2.6, 1, 1.1, 1.1])
rows = []
for lib in LIBS:
    for key, lab in (('aln', 'bwa aln'), ('hry1', 'HRY saved index, run 1'), ('hry2', 'HRY saved index, run 2'),
                     ('hry3', 'HRY saved index, run 3'), ('fly', 'HRY index built in run')):
        s, gb, user, _ = sp[(lib, key)]
        rows.append([lib, lab, '%.2f' % s, '%.1f' % user, '%.1f' % gb])
s, gb, user, _ = time_log(os.path.join(P, 'index.log'))
rows.append(['(all)', 'HRY index build, -r 50', '%.2f' % s, '%.1f' % user, '%.1f' % gb])
table('Every timed run', ['library', 'run', 'wall (s)', 'user CPU (s)', 'peak RSS (GB)'], rows, widths=[1.2, 2, 1, 1, 1], numeric_from=2,
      note='HRY index on disk: 21,413,792,756 bytes (hs37d5.fa.sti).')
rows = []
per = {lib: dict(sa_times(rd('%s.hry.3.log' % lib, P))) for lib in LIBS}
keys = [k for k, _ in sa_times(rd('p17_shuf1m.hry.3.log', P))]
for k in keys:
    rows.append([k] + [per[l].get(k, '-') for l in LIBS])
table('HRY time breakdown, saved index (run 3)', ['step', 'sub100k', 'p17', 'p1'], rows, widths=[2.4, 1, 1, 1],
      note='Per-step lines are summed over the 32 threads; "Total time reading index" and "Total time mapping" are wall-clock.')

h('11. Real libraries: what gets reported as mapped')
rows = []
for line in rd('realmap.txt', P).splitlines():
    m = re.match(r'\s*(\S+)\s+(aln|hry)\s+reads\s+(\d+)\s+\|\s+mapped\s+([\d.]+)%\s+MAPQ>=25\s+([\d.]+)%\s+MAPQ>=30\s+([\d.]+)%\s+\|\s+soft-clipped\s+([\d.]+)%[^|]*(?:\|\s+AS>=60\s+([\d.]+)%\s+\+MAPQ>=25\s+([\d.]+)%\s+\+MAPQ>=30\s+([\d.]+)%)?', line)
    if m:
        g = m.groups()
        rows.append([g[0], 'bwa aln' if g[1] == 'aln' else 'HRY', '{:,}'.format(int(g[2])), g[3], g[4], g[5], g[6]] +
                    ([g[7], g[8], g[9]] if g[7] else ['-', '-', '-']))
table('Mapped fractions (% of all reads)', ['library', 'mapper', 'reads', 'mapped', 'MAPQ≥25', 'MAPQ≥30',
      'soft-clipped (% of mapped)', 'AS≥60', 'AS≥60 + MAPQ≥25', 'AS≥60 + MAPQ≥30'], rows,
      widths=[1.2, 0.9, 1, 0.8, 0.8, 0.8, 1, 0.8, 1, 1], numeric_from=2)
table('p1: HRY MAPQ ≥ 25 records on reads bwa aln leaves unmapped (examples)', ['position', 'MAPQ', 'CIGAR', 'length', 'NM', 'AS'], [
    ['22:46372284', '27', '9S20M18S', '47', '0', '40'], ['1:157958409', '25', '17M28S', '45', '0', '44'],
    ['17:51508179', '25', '17M20S', '37', '0', '44'], ['19:19335811', '25', '17S19M', '36', '0', '48'],
    ['4:819975', '25', '30S19M', '49', '0', '48']], widths=[1.5, 0.7, 1.3, 0.8, 0.6, 0.6],
    note='Of HRY\'s 3,587 MAPQ ≥ 25 records on p1, 72.7% are soft-clipped and 70.8% have AS < 60; only 1,022 of these reads '
         'are mapped by bwa aln at all. MAPQ does not remove such local hits; the paper\'s AS ≥ 60 does.')

h('12. Caveats')
bullets([
    'The simulated reads are all endogenous; exogenous behaviour is only seen on the real libraries (no truth there).',
    'SNVs only (no indels) in the simulation, so placement truth is exact; the paper\'s Mason indels are not modelled.',
    'Divergence here is one person\'s real SNVs (~0.1%), the human-aDNA case. The paper tests human CHM13 at 0.5-5% synthetic '
    'divergence; its CHM13 0.5% win does not reproduce here (HRY -0.69 points at 40-49 bp). Not resolved without running its pipeline.',
    'The paper\'s correctness rule (overlap AND AS ≥ 60) differs; the overlap part is ruled out as the cause (section 6), the '
    'AS part penalises HRY\'s shortest reads here and bwa aln\'s 5+-mismatch alignments there.',
    'strobealign and bwa MAPQ are on different scales; section 5 compares at equal precision.',
    'gpualn and G-ALN speeds are from 2026-11-01 (GPU reserved today); the bwa aln rows were re-measured and agree within 1-7%.',
    'adna_bench\'s compile rule uses cargo build (dev profile, opt-level 1); if the paper\'s timings came from it they understate '
    'HRY\'s speed. Ours is a release build.',
])
table('Files', ['what', 'where'], [
    ['simulated set, SAMs, logs', '/home/teemu/sorsa/adna_simbench/work'],
    ['real-library runs', '/home/teemu/sorsa/adna_simbench/speed'],
    ['scripts', 'cuda/tools/benchmark/ (hry_run.sh, hry_speed.sh, score_sam.py, mapq_curve.py, hybrid_curve.py, '
                'merge_hybrid.py, ref_bias.py, bias_by_len.py, bootstrap_bias.py, clipstats.py, hry_report_odt.py)'],
    ['HRY source and binary', '/home/teemu/sorsa/hry/strobealign-hry (tag hry)'],
    ['authors\' benchmark code', '/home/teemu/sorsa/hry/adna_bench'],
], widths=[1.4, 4], numeric_from=99)

# ---------------------------------------------------------------- document
NS = ('xmlns:office="urn:oasis:names:tc:opendocument:xmlns:office:1.0" '
      'xmlns:style="urn:oasis:names:tc:opendocument:xmlns:style:1.0" '
      'xmlns:text="urn:oasis:names:tc:opendocument:xmlns:text:1.0" '
      'xmlns:table="urn:oasis:names:tc:opendocument:xmlns:table:1.0" '
      'xmlns:fo="urn:oasis:names:tc:opendocument:xmlns:xsl-fo-compatible:1.0" '
      'xmlns:svg="urn:oasis:names:tc:opendocument:xmlns:svg-compatible:1.0"')
FONT = 'Liberation Sans'
def pstyle(name, size, bold=False, align=None, mt='0cm', mb='0.1cm', color=None, italic=False, keep=False, indent=None, family='paragraph'):
    tp = 'fo:font-size="%spt" style:font-name="%s"' % (size, FONT)
    if bold: tp += ' fo:font-weight="bold"'
    if italic: tp += ' fo:font-style="italic"'
    if color: tp += ' fo:color="%s"' % color
    pp = 'fo:margin-top="%s" fo:margin-bottom="%s"' % (mt, mb)
    if align: pp += ' fo:text-align="%s"' % align
    if keep: pp += ' fo:keep-with-next="always"'
    if indent: pp += ' fo:margin-left="%s" fo:text-indent="-%s"' % (indent, indent)
    return ('<style:style style:name="%s" style:family="paragraph"><style:paragraph-properties %s/>'
            '<style:text-properties %s/></style:style>' % (name, pp, tp))
styles = ''.join([
    pstyle('Body', 10, mb='0.15cm'),
    pstyle('Subtitle', 9, italic=True, color='#555555', mb='0.3cm'),
    pstyle('Bullet', 10, mb='0.08cm', indent='0.4cm'),
    pstyle('Caption', 9, bold=True, mt='0.3cm', mb='0.1cm', keep=True),
    pstyle('Note', 8, italic=True, color='#444444', mt='0.08cm', mb='0.25cm'),
    pstyle('TH', 8, bold=True, align='center', mb='0cm'),
    pstyle('TD', 8, mb='0cm'), pstyle('TDB', 8, bold=True, mb='0cm'),
    pstyle('TDnum', 8, align='end', mb='0cm'), pstyle('TDnumB', 8, bold=True, align='end', mb='0cm'),
    pstyle('H0', 17, bold=True, mt='0cm', mb='0.15cm', color='#1F3864', keep=True),
    pstyle('H1', 13, bold=True, mt='0.45cm', mb='0.15cm', color='#1F3864', keep=True),
    pstyle('H2', 11, bold=True, mt='0.3cm', mb='0.1cm', keep=True),
])
cell = ('<style:style style:name="%s" style:family="table-cell"><style:table-cell-properties fo:padding="0.05cm" '
        'fo:border="0.5pt solid #9AA5B1"%s/></style:style>')
auto = ''.join([
    '<style:page-layout style:name="PL"><style:page-layout-properties fo:page-width="21cm" fo:page-height="29.7cm" '
    'style:print-orientation="portrait" fo:margin-top="1.6cm" fo:margin-bottom="1.6cm" fo:margin-left="1.7cm" '
    'fo:margin-right="1.7cm"/></style:page-layout>',
    '<style:style style:name="Tbl" style:family="table"><style:table-properties style:width="17.6cm" '
    'table:align="margins" fo:margin-top="0cm" fo:margin-bottom="0cm"/></style:style>',
    cell % ('HdrCell', ' fo:background-color="#DCE6F1"'), cell % ('Cell', ''), cell % ('CellAlt', ' fo:background-color="#F4F6F9"'),
] + autostyles)
doc = ('<?xml version="1.0" encoding="UTF-8"?>\n<office:document %s office:version="1.3" '
       'office:mimetype="application/vnd.oasis.opendocument.text">'
       '<office:font-face-decls><style:font-face style:name="%s" svg:font-family="\'%s\'" '
       'style:font-family-generic="swiss"/></office:font-face-decls>'
       '<office:styles>%s</office:styles><office:automatic-styles>%s</office:automatic-styles>'
       '<office:master-styles><style:master-page style:name="Standard" style:page-layout-name="PL"/></office:master-styles>'
       '<office:body><office:text>%s</office:text></office:body></office:document>\n'
       % (NS, FONT, FONT, styles, auto, ''.join(body)))
with open(OUT, 'w') as f:
    f.write(doc)
print('wrote %s: %d tables' % (OUT, ntab[0]))
