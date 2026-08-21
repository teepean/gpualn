# gpualn — GPU BWA-backtrack for ancient DNA

`gpualn` runs the **`bwa aln` (BWA-backtrack) bounded-error search on the GPU** and emits a
`.sai` file that is **byte-identical to CPU `bwa aln`**. It is a standalone project derived from
[bwa](https://github.com/lh3/bwa) 0.7.19-r1273 by Heng Li; the rest of bwa (`index`, `mem`, `samse`,
`sampe`, `bwasw`, …) is kept intact and works exactly as upstream.

The target workload is ancient DNA: millions of short (30–90 bp) reads searched with
`-l 1024 -n 0.01 -o 2`, where the overwhelming majority of reads have **no hit** and the unbounded
backtracking over the FM-index is what costs all the time.

> gpualn is not endorsed by or affiliated with the bwa project. Report gpualn issues here, not to
> upstream bwa. See `NOTICE` for attribution and `COPYING` for the licence (GPLv3).

## Why it is fast

Two engines, both bit-exact at the output level:

- **Exact engine** — a per-read iterative DFS over the same node set bwa's best-first queue visits
  when the `max_entries` cap is never reached (measured: it never is, on aDNA data). It answers only
  the question "does this read have any hit?", which is the cheap half of the work.
- **Search-scheme engine** (default) — a **bidirectional** FM-index search using optimal
  search schemes (Kianfar for K≤2, SeqAn3 for K=3, Renders et al. for K=4). The bidirectional index
  is free here: gpualn searches `S = T · revcomp(T)`, so `revcomp(S) = S` and one index serves both
  directions. Up to **20x** the exact engine's kernel time on short-read sets.

Correctness rests on a **superset contract**: the GPU only has to guarantee
`has_hit_gpu(r) ≥ has_hit_bwa(r)`. Every read the GPU flags is then re-aligned on the CPU by the
unmodified `bwt_match_gap`, so the emitted alignments are bwa's own. A false positive costs time;
a false negative would cost correctness, and the startup covering-check plus ~4.3 M reads of
byte-identical validation guard against it.

## Requirements

- NVIDIA GPU, compute capability **8.0+** (developed on RTX 3090 / sm_86), ~4 GB free VRAM for the
  hs37d5 `.bwt`.
- CUDA toolkit 12 or 13 (`nvcc`), zlib, pthreads.

## Build

```sh
make                 # -> ./gpualn       (CUDA; includes the aln-gpu subcommand)
make cpu             # -> ./gpualn-cpu   (no CUDA; plain bwa, builds anywhere)
make CUDA_ARCH=sm_90 # target another architecture
```

## Use

```sh
./gpualn index ref.fa                                    # same index as upstream bwa
./gpualn aln-gpu -l 1024 -n 0.01 -o 2 -t 16 \
    ref.fa reads.fq > out.sai
./gpualn samse ref.fa out.sai reads.fq > out.sam
```

`aln-gpu` takes the same options as `bwa aln` and is a drop-in replacement for it. `-x` picks the
engine — `scheme` (default) or `exact`, the legacy per-read DFS kept for A/B measurement; the `.sai`
is bit-exact either way. `gpualn` is accepted as an alias for the subcommand so existing scripts
keep working. `-t` sets the number of
CPU threads used for the reconcile step, which is the bottleneck once the kernel is fast — give it
every core you have.

### Environment knobs

| Variable | Effect |
| --- | --- |
| `GPUALN_HISTO=1` | print per-stage timings, kernel seconds and hit histograms to stderr |
| `GPUALN_NOORDER=1` | disable longest-first read ordering (for A/B measurement) |
| `GPUALN_BATCH=n` | reads per batched CPU reconcile call (default 8) |

`GPUALN_SCHEME` is obsolete — the engine it selected is now the default. Setting it prints a notice
and changes nothing.

## Measured results

Reference hs37d5, RTX 3090 + 16 CPU threads, all outputs byte-identical to CPU `bwa aln`:

| read set | reads | `-x exact` | `-x scheme` (default) |
| --- | --- | --- | --- |
| 2730 short (30–63 bp) | 671,652 | 22.2 s | **1.1 s** (20.2x) |
| end-to-end `aln` incl. CPU reconcile | 671,652 | 24.3 s | **4.6 s** (5.3x) |

Against single-threaded CPU `bwa aln` on the same set (235 s), the scheme engine end-to-end is
**~51x**. On real whole-library data the speedup tracks the **unmapped fraction**, not read length:
once the kernel is this fast, the CPU reconcile of the reads that *do* hit dominates.

See `cuda/PROGRESS.md` for the full engineering log, including the ideas that were measured and
**rejected** (12-byte node packing, per-lane MLP, `__ldcs`, `-march=native`, e-graph rewriting).

## Layout

- `cuda/` — everything new: device FM-index, both engines, scheme tables, docs, CompileIQ tuner.
- everything else — bwa, essentially unmodified. `bwtgap.c` carries one upstream-safe optimisation
  (prefetching the next queue entry's Occ buckets, +13% on plain `bwa aln`) and a batched entry point.
- `README-bwa.md` — upstream bwa's README, preserved verbatim.

## Licence

GPLv3 (`COPYING`), inherited from the bwa distribution. Individual upstream files keep their own
headers; new files under `cuda/` are GPLv3, Copyright (C) 2026 teepean.
