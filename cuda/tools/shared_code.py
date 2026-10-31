"""Lines and identifiers that G-ALN's own code shares with gpualn's code but NOT with upstream bwa.
usage: shared_code.py GALN_SRC OURS_DIR UPSTREAM_DIR [EXTRA_UPSTREAM_DIR...]
UPSTREAM_DIR = upstream bwa, e.g. `git archive origin/master | tar -x -C DIR` (lh3/bwa 0.7.19);
pass G-ALN's own src/*/bwa_upstream dirs as extras. Used for GALN_COMPARISON.md "Provenance"."""
import os, re, sys
EXT = ('.c', '.h', '.cu', '.cuh', '.cpp')
def files(root, skip=()):
    for d, _, fs in os.walk(root):
        if any(s in d for s in skip): continue
        for f in fs:
            if f.endswith(EXT): yield os.path.join(d, f)
def norm(line): return re.sub(r'\s+', ' ', line.strip())
TRIVIAL = re.compile(r'^[\s{}();,]*$|^(#include|#endif|#else|return|break;|continue;|else|\} else)')
def lines(paths):
    out = {}
    for p in paths:
        for i, l in enumerate(open(p, errors='replace'), 1):
            n = norm(l)
            if len(n) >= 30 and not TRIVIAL.match(n): out.setdefault(n, (p, i))
    return out
IDENT = re.compile(r'\b[A-Za-z_][A-Za-z0-9_]{5,}\b')
def idents(paths):
    s = set()
    for p in paths: s.update(IDENT.findall(open(p, errors='replace').read()))
    return s
galn = list(files(sys.argv[1], skip=('bwa_upstream', 'third_party')))
ours = list(files(sys.argv[2]))
up = [p for r in sys.argv[3:] for p in files(r)]
G, O, U = lines(galn), lines(ours), set(lines(up))
shared = sorted(k for k in G if k in O and k not in U)
print(f"G-ALN own files: {len(galn)}, gpualn files: {len(ours)}, upstream files: {len(up)}")
print(f"distinct non-trivial lines: G-ALN {len(G)}, gpualn {len(O)}")
print(f"\n=== {len(shared)} lines shared by G-ALN and gpualn, absent from upstream bwa ===")
for k in shared: print(f"  {k[:110]}\n      galn: {os.path.relpath(G[k][0], sys.argv[1])}:{G[k][1]}   ours: {os.path.relpath(O[k][0], sys.argv[2])}:{O[k][1]}")
gi, oi, ui = idents(galn), idents(ours), idents(up)
common = sorted(x for x in gi & oi - ui if not x.startswith(('cuda', '__', 'std', 'uint', 'int', 'size_t', 'unsigned', 'printf', 'fprintf', 'memcpy', 'memset', 'malloc', 'calloc', 'realloc', 'pthread', 'atomic', 'const', 'static', 'struct', 'return', 'define', 'include', 'ifndef', 'endif', 'typedef', 'extern', 'inline', 'double', 'float', 'sizeof', 'stderr', 'stdout', 'exit', 'assert', 'blockIdx', 'blockDim', 'threadIdx', 'gridDim', 'warpSize')))
print(f"\n=== {len(common)} identifiers (>=6 chars) in both, absent from upstream bwa and CUDA/libc ===")
print('  ' + ' '.join(common))
