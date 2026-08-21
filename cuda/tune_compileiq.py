#!/usr/bin/env python3
"""Tune ptxas Advanced Controls for a gpualn DFS kernel with NVIDIA CompileIQ.

The kernel sits on a jagged codegen landscape (PROGRESS.md: maxrregcount 48 -> 6,445 r/s,
40 -> 5,243, default -> 9,934), which is exactly what CompileIQ's search over the undocumented
ptxas scheduler controls is for. Source-level tuning has plateaued, so this is the remaining
CUDA-13-specific lever (nothing else in 13.x targets sm_86).

Objective: GPU-kernel seconds (lower is better), taken as the best of N runs to suppress thermal
drift. Hard correctness gate: any config whose .sai md5 differs from the golden value scores
INVALID, so a "fast" but wrong schedule can never win.

ONE ACF PER ENGINE. The schedules do not transfer: the ACF tuned for `k_dfs_warp2` gave 1.174x
there and ~1.0x on `k_dfs_scheme`. Pick the engine with --engine; the default is `scheme`, matching
the binary's own default.

TUNE ON REPRESENTATIVE READS. `k_dfs_scheme` is ONE kernel that contains both hot paths: reads with
K>SCH_MAXK or len>SCH_MAXLEN take the `d_dfs_has_hit_warp2` branch inlined into the same __global__
(scheme_engine.cuh:284). One register allocation, one ptxas schedule, both paths. So the fallback
rate of the tuning set is not a contaminant to be eliminated -- it is a property of the workload the
ACF is being fitted to, and it should match production. Tuning at 0% fallback and shipping to a set
with 10% fits the ACF to only half the kernel. The measured rate is printed before the search starts.

Run:   /home/teemu/sorsa/CompileIQ/.venv/bin/python cuda/tune_compileiq.py [--engine scheme|exact]
       [--generations G] [--pool P] [--repeats R] [--reads FQ --golden MD5]
Apply: nvcc ... -Xptxas=--apply-controls=cuda/gpualn_scheme.acf
"""
import argparse, hashlib, os, re, subprocess
from uuid import uuid4

from compileiq.ciq import Search
from compileiq.search_spaces.compilers import PtxasSearchSpace
from compileiq.types import INVALID_SCORE, SearchConfiguration
from compileiq.utils.helpers import save_compiler_config

REPO   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REF    = "/home/dnastorage/aDNApipeline/hs37d5.fa"
READS  = os.path.join(REPO, "test_data/sub100k.fq")
GOLDEN = "eecf35c15db452bfe0500ce0b9dbd723"       # sub100k .sai, == CPU bwa aln
OBJS   = ["bwtgap.o", "bwtaln.o", "bwaseqio.o", "bamlite.o", "bwase.o"]
ACF    = {"scheme": "cuda/gpualn_scheme.acf", "exact": "cuda/gpualn.acf"}

ARGS = None


def build(acf, exe):
    """acf=None builds with the default schedule (ptxas rejects an empty controls file)."""
    cmd = ["nvcc", "-O3", "-std=c++14", "-arch=sm_86", "-DALN_GPU_MAIN", "-I."]
    if acf:
        cmd.append(f"-Xptxas=--apply-controls={acf}")
    cmd += ["cuda/aln_gpu.cu", *OBJS, "-o", exe, "-L.", "-lbwa", "-lm", "-lz", "-lpthread", "-lrt"]
    r = subprocess.run(cmd, cwd=REPO, capture_output=True, text=True, timeout=600)
    return r.returncode == 0


def run_once(exe):
    """returns (kernel_seconds, md5, fallback_percent); any element may be None"""
    r = subprocess.run([exe, "-l", "1024", "-n", "0.01", "-o", "2", "-t", "16",
                        "-x", ARGS.engine, REF, ARGS.reads],
                       capture_output=True, timeout=1200, cwd=REPO,
                       env=dict(os.environ, GPUALN_HISTO="1"))
    if r.returncode != 0:
        return None, None, None
    err = r.stderr.decode(errors="replace")
    t = re.search(r"GPU-kernel ([0-9.]+)s", err)
    fb = re.search(r"scheme-fallback \d+ \(([0-9.]+)%\)", err)
    return (float(t.group(1)) if t else None), hashlib.md5(r.stdout).hexdigest(), \
           (float(fb.group(1)) if fb else 0.0)


def objective(config) -> float:
    uid = uuid4().hex[:8]
    acf, exe = f"/tmp/ciq_{uid}.acf", f"/tmp/ciq_{uid}.bin"
    try:
        if config is None:
            acf = None
        else:
            save_compiler_config(acf, config)
        if not build(acf, exe):
            return INVALID_SCORE
        best = None
        for _ in range(ARGS.repeats):
            t, md5, _ = run_once(exe)
            if t is None or md5 != ARGS.golden:  # correctness gate
                return INVALID_SCORE
            best = t if best is None else min(best, t)
        return best
    except Exception:
        return INVALID_SCORE
    finally:
        for f in (acf, exe):
            if not f: continue
            try: os.remove(f)
            except OSError: pass


def probe_workload():
    """Validate the reads/golden pair and report what mix of the kernel the search will be fitting."""
    exe = f"/tmp/ciq_probe_{uuid4().hex[:8]}.bin"
    try:
        if not build(None, exe):
            print("[tune] probe build FAILED -- aborting"); return False
        t, md5, fb = run_once(exe)
        if t is None:
            print("[tune] probe run FAILED -- aborting"); return False
        if md5 != ARGS.golden:
            print(f"[tune] probe md5 {md5} != golden {ARGS.golden} -- wrong reads/golden pair, aborting")
            return False
        if ARGS.engine == "scheme":
            print(f"[tune] {os.path.basename(ARGS.reads)}: {fb:.2f}% of reads take the fallback branch "
                  f"inside k_dfs_scheme (baseline {t:.3f}s)")
            if fb == 0.0:
                print("[tune] warning: no fallback reads here, so the ACF will be fitted to the scheme "
                      "path only; production aDNA sets run ~10% fallback (K>=5)")
        else:
            print(f"[tune] {os.path.basename(ARGS.reads)}: baseline {t:.3f}s")
        return True
    finally:
        try: os.remove(exe)
        except OSError: pass


def main():
    global ARGS
    ap = argparse.ArgumentParser()
    ap.add_argument("--engine", choices=("scheme", "exact"), default="scheme",
                    help="which kernel to tune; schedules do not transfer between them")
    ap.add_argument("--scheme", action="store_true", help=argparse.SUPPRESS)  # back-compat alias
    ap.add_argument("--generations", type=int, default=4)
    ap.add_argument("--pool", type=int, default=16)
    ap.add_argument("--repeats", type=int, default=2, help="runs per config; score is the best")
    ap.add_argument("--reads", default=READS)
    ap.add_argument("--golden", default=GOLDEN, help="expected .sai md5 for --reads")
    ap.add_argument("--out", default=None, help="ACF path [cuda/gpualn_scheme.acf or cuda/gpualn.acf]")
    ARGS = ap.parse_args()
    if ARGS.scheme:
        ARGS.engine = "scheme"
    if ARGS.out is None:
        ARGS.out = os.path.join(REPO, ACF[ARGS.engine])

    ver = re.search(r"release (\d+\.\d+),", subprocess.run(
        ["ptxas", "--version"], capture_output=True, text=True, check=True).stdout).group(1)
    print(f"[tune] ptxas {ver}, engine={ARGS.engine}, "
          f"{ARGS.generations} generations x pool {ARGS.pool}, {ARGS.repeats} repeats")

    if not probe_workload():
        return

    base = objective(None)                     # default schedule, same measurement path
    if not isinstance(base, float):
        print("[tune] baseline FAILED (build/run/md5) -- aborting"); return
    print(f"[tune] baseline (no controls): {base:.3f}s")

    tuner = Search(objective_function=objective,
                   search_space=PtxasSearchSpace(version=ver),
                   search_config=SearchConfiguration(problem_type="min",
                                                     generations=ARGS.generations,
                                                     pool_size=ARGS.pool))
    results = tuner.start(num_workers=1)       # serialise: one GPU, timings must not contend
    best = results.get_best_result()
    sc = best.get("score_1") if best else None
    if not isinstance(sc, (int, float)):
        print(f"[tune] no valid configuration found (best={sc!r}); keeping the default schedule")
        return
    print(f"[tune] best {sc:.3f}s vs baseline {base:.3f}s = {base/sc:.3f}x")
    if sc < base:
        save_compiler_config(ARGS.out, best["params"])
        print(f"[tune] wrote {ARGS.out}  (apply with -Xptxas=--apply-controls={ARGS.out})")
    else:
        print("[tune] best is not better than the default schedule; nothing written")


if __name__ == "__main__":
    main()
