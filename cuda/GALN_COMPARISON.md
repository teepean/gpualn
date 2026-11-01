# G-ALN: an independent GPU `bwa aln`, and what it means for gpualn

**Status: measured 2026-11-01** (RTX 3090 + Ryzen 9 9950X; sections "Measured head-to-head" onward).
The sections before them are the 2026-10-31 source analysis, unchanged except where marked.
**Submitted the same day:** the prefilter as <https://github.com/bilalbioinfo/galn/pull/3>, the
page-cache problem as <https://github.com/bilalbioinfo/galn/issues/2>.

**Bottom line:** on our box gpualn beats G-ALN 0.2.4 on all three libraries we tried, including a
94%-endogenous one, and G-ALN with gpualn's scheme engine as a prefilter beats both on all three
(2.1-7.9x over G-ALN). The combination is submitted to G-ALN as a pull request (see the end).

G-ALN (<https://github.com/bilalbioinfo/galn>, GPL-3.0, v0.2.4) is a GPU `bwa aln` + `bwa index`
with byte-identical output to BWA 0.7.19, from the Dalén group (Centre for Palaeogenetics):
Sharif, Hedfors, van der Valk, Heintzman, Dalén, *G-ALN: a GPU-accelerated BWA-aln algorithm for
palaeogenomics*, bioRxiv 2026.10.18.759956 (posted 2026-10-21). Same goal and same user community
as gpualn. The preprint says it was "developed through AI-assisted coding (Claude Opus 5.5)".

Sources: the code (read directly, `/home/teemu/sorsa/galn-upstream/galn`, commit `1e27f17`) and the
bioRxiv full text (read via a fetch; the figures below quoted from it, not independently checked).

## The one difference that matters

| | gpualn | G-ALN |
| --- | --- | --- |
| what the GPU computes | **existence only** -- a superset filter (`has_hit`) | the **final alignments** -- `bwt_match_gap` ported line by line |
| what the CPU does | re-runs `bwt_match_gap` for **every** flagged read | only reads the GPU overflows or hands off, from the other end of a shared list |
| no-hit reads | search schemes prove "no hit" cheaply | full faithful bounded search |

That is exactly the lever this project identified and stopped at. Phase 15: on a high-endogenous
sample the CPU reconcile was **98.9%** of the short-read stage (1,094 s vs a 25.7 s kernel), and the
recorded corollary was *"for high-endogenous samples the lever is the reconcile, not the DFS --
either move `bwt_match_gap` for flagged reads onto the GPU, or avoid a full CPU re-alignment per
hit."* G-ALN is the first option, built. Every kernel optimisation since (scheme engine default,
ACF, CAP_SM 1280) has been on the 25.7 s, not the 1,094 s.

Conversely, Phase 15's rule -- *the scheme engine's benefit scales with the fraction of reads that
DO NOT map* (20x kernel / 51x vs CPU at 0.5% mapping) -- is where gpualn may still be ahead: G-ALN
does a full faithful search on every no-hit read. **Expected (2026-10-31):** gpualn wins on
low-endogenous libraries, G-ALN wins on high-endogenous ones. **Measured (2026-11-01): only the first
half held** -- gpualn also wins at 94% endogenous, because even there the unmapped reads are most of
the work (below).

## How G-ALN makes faithful GPU generation pay, where our June prototype did not

Memory note `gpu-align-generation-investigation` (2026-06) recorded GPU alignment generation as
bit-exact but "modest", and warp-cooperative generation as "a dead end": cooperation only pays on
big search trees, but the `c1>1` gate sent exactly the big-tree reads to the CPU -- *"gate and
warp-cooperation are fundamentally opposed"*. G-ALN avoids both halves of that:

1. **No `c1>1` gate at all.** The port preserves upstream's pop/push order exactly (best-first
   buckets, LIFO within a bucket), so order-sensitive multi-mappers come out identical too.
2. **Light tier = our "Route 1".** One thread per read, plain `bwtgap_match_gap`, stack in global
   memory (`--light-cap` 4096 entries). The code comment says ~80% of reads finish here. Our June
   note left the real-index per-thread ratio unmeasured; G-ALN is evidence it is adequate,
   though the preprint does not isolate it (see below).
3. **The warp engine is not tree-parallel** (`src/gapw.cuh`). All 32 lanes run the *same* scalar
   search redundantly in lockstep, so every branch is warp-uniform and no lane is ever idle; the
   lanes split only the memory work -- lanes 0-15 load the 64 B BWT block for `k-1`, lanes 16-31
   for `l`, in one round trip. And it is routed to **high-effort** reads (light-tier overflow): the
   inverse of our gate. It gates on effort, and sends the big trees *to* the warp.
4. **16-byte stack entries** (upstream's are 32): `k`,`l` in 36 bits each, the score implicit in
   the bucket, so a warp can hold more than bwa's default 2 M-entry queue. Any overflow returns an
   error code and the read is re-run from scratch on the CPU, so exactness never depends on a cap.
5. **GPU/CPU share one two-ended list** of heavy reads: warps claim from the back, CPU workers from
   the front, frontiers exchanged ~1/ms with small async copies (no unified memory). Idle CPU
   workers also re-run reads the GPU is still on; first finished copy wins (identical by
   construction).

## Reading the headline number

**6.3x** = 643 s for one full Arrhenius node (4x GH200) vs 4,027 s for `bwa aln -t 256` on a
2x EPYC 9755 node, 50 M mammoth reads (MD090) on the 3.4 Gbp Asian elephant genome,
`-l 1024 -n 0.01 -o 2`. Index: 32x (82 s on one GH200 vs 2,617 s; 19 s in v0.2.4).

That is a **node-to-node** comparison, which is a fair framing for HPC users, but not a GPU
speedup: `--cpu-workers` defaults to *"all cores of the job minus one per GPU"*, so the GPU node's
Grace CPUs (72 cores per GH200) are also searching, with an interleaved CPU search faster per core
than stock bwa. The preprint gives no `--gpu-only` timing and no per-GPU figure, and nothing on a
single consumer GPU with a desktop CPU -- which is our configuration. It measures no accuracy or
reference bias, which is correct for a byte-identical tool: its accuracy *is* `bwa aln`'s, which
`BENCHMARK_ACCURACY.md` measured (93.241% on the simulated set).

## What is and is not new to us

| G-ALN idea | status here |
| --- | --- |
| faithful GPU generation of hits | **the big one.** Not done; Phase 15 named it as the lever |
| CPU interleave, 8 reads/core (`aln_il.h`) | **already done** -- Phase 20 batched reconcile, measured 1.10x at 16 threads, 1.02x at 32. Not a new lever |
| warp lockstep, lanes split occ loads | new mechanism; distinct from the tree-parallel collect we refuted |
| 16 B stack entries | our `SNode` is 24 B because it is bidirectional (`x0,x1` 33 bits each at hs37d5 scale, `x2` < 2^31, `p` 31 bits = 128 bits exactly). 16 B is feasible with zero slack. Shared stack depth was worth 1.27x in Phase 25; 1.5x density at the same shared memory is a candidate, payoff unknown |
| two-ended GPU/CPU list | only relevant once the GPU generates alignments |

A combined design would take both halves: the scheme engine to discard no-hit reads cheaply, then
G-ALN-style faithful generation for the flagged ones on the GPU. Licences are compatible (both
GPL-3.0, both derived from bwa).

Scope, for completeness: G-ALN supports nearly every `bwa aln` option (not `-q`, `-b/-0/-1/-2`,
`-I`, `-Y`, `-B`), gzip input, multiple GPUs with NUMA-aware CPU workers, `samse`/`sampe` writing
BAM, and GPU index construction. gpualn's `aln-gpu` parses `-l -n -o -t -f -S -r -x`.

## On our hardware (static, from their code)

`gapw_smem` is 7,916 B per warp (the width/bid arrays are ~6 KB of it); the warp kernel runs 4
warps per block = 31.7 KB, grid sized by `cudaOccupancyMaxActiveBlocksPerMultiprocessor`. Shared
memory alone allows at most **3 blocks = 12 warps/SM on sm_86** (100 KB/SM) vs 7 blocks = 28 on a
GH200 (228 KB/SM); registers may bind lower. With 82 SMs vs 132 and GDDR6X vs HBM3, the warp tier
should be much weaker on an RTX 3090 than on the hardware it was tuned for. It affects only the
~20% of reads that overflow the light tier, which itself needs no shared memory.

## Side result: Phase 20's `.sai` `-t` mystery, solved

G-ALN's README says `bwa aln` writes its `-t` into the `.sai` header. It does: `bwtaln.c:179`
writes the whole `gap_opt_t`, which contains `int n_threads`. Verified on CPU (sub2k, `-t 2` vs
`-t 4`): `cmp -l` reports **exactly one differing byte, offset 57** = 4 (magic) + 52 (offset of
`n_threads`) + 1, with values `2` and `4`. Every alignment record is thread-count invariant.

So Phase 20's suspect (the per-thread `w` buffer) was wrong, and its consequence 3 -- that
switching to `-t 32` "will change the `.sai` by a byte or so" and needs golden regeneration -- is
only a header change. The `-t 32` reconcile speedup (1.36x) moves no alignment. The header is
68 bytes (4 magic + `sizeof(gap_opt_t)` = 64, both checked with the compiler), so goldens can be
compared past it -- `cmp <(tail -c +69 a.sai) <(tail -c +69 b.sai)`, verified identical for the
`-t 2` / `-t 4` pair -- or regenerated. G-ALN solves the same thing with `-dt` (dummy threads: it
sets only the header value).

## Provenance: did G-ALN build on gpualn? (checked 2026-10-31)

**No sign that it did.** The question arises because both were written with Claude models and
gpualn was public first:

| date | event |
| --- | --- |
| 2026-06-01 | `cuda-port` branch created on public `teepean/bwa` |
| 2026-06-10 | `0f19581` pushed: faithful best-first GPU prototype + generation scoping doc |
| 2026-08-21 | `teepean/gpualn` created, public |
| 2026-10-16 | G-ALN's first public commit -- one squashed import of v0.2.1; its history is not public |
| 2026-10-21 | G-ALN preprint posted |

So the timeline allows it. The code does not show it:

- **Shared text.** Every non-trivial line (>= 30 chars) of G-ALN's own source (29 files;
  `bwa_upstream/` and `third_party/` excluded) against all 71 C/CUDA files in this repo, minus every
  line of upstream bwa 0.7.19 (`cuda/tools/shared_code.py`): **7 lines in common** -- four lines of
  generic CUDA/C boilerplate, a stock `CK` error-check macro that both codebases spell more than one
  way, and two mechanical ports of bwa's own code (`cntl[i] = cntk[i]` is bwa's `memcpy` in
  `bwt_2occ4` unrolled; `primary; /* S^{-1}(0) */` is bwa's `bwt.h` comment). No gpualn identifier
  appears (`SNode`, `CAP_SM`, `d_scheme_has_hit`, `sch_validate`, `gap_run_*`, ...).
- **Shared design: none distinctive.** gpualn's defining parts -- the superset filter with CPU
  reconcile, search schemes on the bidirectional index, shared-memory stacks, the `c1>1` gate -- are
  absent from G-ALN. G-ALN's -- the lockstep warp engine, 16 B chunked buckets, the two-ended list,
  the FNV pop-trace hash -- are absent from gpualn, including from the June prototype (a toy with
  its own suffix-array builder and `uint32` index, fixed per-bin arrays, and a collect-and-sort warp
  kernel). The one idea in both, batched CPU search with prefetch, is QuadRank-style batching that
  Phase 16 itself took from the literature, and is structured differently (our prefetch-then-step
  `gap_run_*` vs their A/B halves split at the occ lookup).
- **The preprint** credits BarraCUDA as the earlier GPU port of bwa aln and does not mention
  gpualn, search schemes, or any other GPU bwa-aln (its reference entries 15-22 were not visible).
- Its README's explanation of the `.sai` `-t` byte is right where our public Phase 20 was wrong;
  working from our notes would have imported the wrong theory.

Limits: this rules out reuse of code and design, not that anyone ever looked at the repository.
And since both were written with Claude on the same upstream code, shared idiom (the `CK` macro,
"byte-identical") is expected even with no contact at all -- only distinctive design would mean
anything, and there is none.

## Measured head-to-head (2026-11-01)

RTX 3090 (24 GB, sm_86) + Ryzen 9 9950X (16 cores / 32 threads), hs37d5, `-l 1024 -n 0.01 -o 2`,
page cache warm. Every configuration uses the whole box (`bwa -t 32`, `gpualn -t 32`; G-ALN always
does). **Every run produced the same `.sai` as stock `bwa aln -t 32`** (whole-file md5; G-ALN with
`-dt 32`). Medians of 3 runs for gpualn / G-ALN / G-ALN + prefilter, single runs otherwise
(run-to-run spread 1-5%; `cuda/tools/bench_galn.sh`).

| | sub100k | p17 1M | p1 1M |
| --- | --- | --- | --- |
| reads / mapped | 100,000 / 0.5% | 999,552 / 94.4% | 1,000,528 / 0.2% |
| `bwa aln -t 32` (0.7.18) | 15.3 s | 56.1 s | 149.1 s |
| gpualn-cpu `aln -t 32` (Phase 16 prefetch) | 13.5 s | 50.6 s | 133.5 s |
| G-ALN `--cpu-only` (interleaved) | 13.1 s | 47.3 s | 134.0 s |
| **gpualn** (scheme engine + CPU reconcile) | **2.22 s** | 17.4 s | 7.77 s |
| G-ALN 0.2.4 (GPU + 31 CPU workers) | 7.03 s | 19.9 s | 55.2 s |
| G-ALN 0.2.4 `--gpu-only` | 20.2 s | 39.3 s | 96.5 s |
| **G-ALN + scheme prefilter** | **1.97 s** | **9.63 s** | **6.95 s** |
| G-ALN + scheme prefilter, `--gpu-only` | 12.9 s | 24.0 s | 9.6 s |

Datasets: `sub100k` = head of AVA1B (30-91 bp, the Phase 0 set; golden md5 `eecf35c1` at `-t 16`);
`p17`/`p1` = random 1 M subsamples of two pipeline short-branch BAMs (30-49 bp, mean 39.5),
shuffled with a seeded RNG (`test_data/p17_shuf1m.fq`, `p1_shuf1m.fq`). The prefilter rows use the
0.2.4 code base; the same patch on 0.2.7 measured 7.1 -> 2.1 s, 20.2 -> 9.9 s, 56.0 -> 7.3 s.

Read it this way:
- **The three CPU paths are within 20% of each other at 32 threads.** G-ALN's interleaved search is
  fastest on p17 (1.19x stock bwa, 1.07x our Phase 16 prefetch build) and ties ours on p1.
- **G-ALN's GPU alone is weak on a 3090**: `--gpu-only` is 1.7-2.9x slower than its default, because
  its 31 CPU workers take the giant reads (`--gpu-only` on sub100k spent 6.8 s on 5 such reads in its
  phase G). The preprint's GH200 nodes have 72 Grace cores per GPU doing the same.
- **gpualn beats G-ALN 0.2.4 on all three**, by 3.2x, 1.14x and 7.1x.

## Where the time goes: unmapped reads, even at 94% endogenous

p17 split by whether bwa finds a hit (`cuda/tools/split_by_sai.py`):

| | 943,731 reads with a hit | 55,821 without |
| --- | --- | --- |
| `bwa aln -t 32` | 15.4 s | **41.4 s** |
| G-ALN 0.2.4 | 6.2 s | **16.7 s** |
| G-ALN + prefilter | 6.8 s | 3.6 s |
| gpualn | 15.9 s | 4.0 s |

The 5.6% of reads that do not map take **73% of bwa's time and 73% of G-ALN's**. A read with no hit
is the worst case for `bwt_match_gap`: nothing ends the search, so every path within `max_diff` is
enumerated. That is what Phase 15's "benefit scales with the non-mapping fraction" understated: the
non-mapping fraction is weighted by ~45x per read (bwa: 742 vs 16 us per read). G-ALN's faithful engines have no way around it;
the scheme engine's whole job is to prove "no hit" cheaply. And the hit reads are where G-ALN is
strong (6.2 s vs gpualn's 15.9 s of CPU reconcile) -- the two tools are good at opposite halves.

## The combination: gpualn's scheme engine as a prefilter inside G-ALN

So the obvious design is both halves, and it is now built: `src/prefilter.cuh` in a G-ALN branch
(`/home/teemu/sorsa/galn-upstream/galn-pf027`, branch `scheme-prefilter-027` on their `dev-0.2.7`;
the patch is also in `cuda/patches/`). Before G-ALN's light tier, per GPU:

1. **Exact-occurrence check**, one thread per read: a read that occurs exactly has a hit, so it is
   passed on without a scheme search. Added after measuring the first version: on the hit-only half
   of p17 the scheme search alone cost 1.65 s (84% of p17's mapped reads occur exactly); with the
   check, 0.54 s.
2. **Scheme search** (the K=1..4 tables and the bidirectional executor of `scheme_engine.cuh`,
   verbatim apart from names) on the rest, one warp per read, longest first. A read whose searches
   all finish without a hit gets an empty `.sai` record and is never searched again.
3. Everything else -- hits, budget or stack overflow, K = 0 or >= 5, length > 127, CUDA errors -- is
   "maybe" and goes through G-ALN unchanged.

Gates, each one measured or required for exactness:
- `BWA_MODE_GAPE` only (with `-e` bwa stops counting gap extensions against `max_diff`; the executor
  counts them, so it would no longer be a superset) -- off with `-e`, verified identical.
- Only reads no longer than `-l`. With seeding on, G-ALN's exact search is already cheap (1.3-2.3 s
  per 1 M reads), and the filter, which ignores the seed, made it slower (p1 at `-l 20`: 1.5 -> 6.3 s
  of search). With the gate, seeded runs are unchanged within run-to-run noise; aDNA settings filter
  every read.
- Index guard: the node keeps interval sizes in 32 bits, so the filter turns off if a base occurs
  2^32 times or more in the BWT (references over ~4.3 Gbp in the worst composition; G-ALN 0.2.6+
  supports up to ~137 Gbp). Human is far below.
- Kernel failures can only yield "maybe": verdicts start as "maybe" on the device and launch errors
  are checked (a launch error is not reported by the stream synchronisation).

Correctness of the final build:
- whole-file md5 = stock `bwa aln` on sub2k/sub10k (= goldens `4b068014`/`54410c75`), sub100k,
  p17 1M, p1 1M, each in default and `--gpu-only` mode, plus the p17 hit/no-hit halves, a file
  where every read is filtered, and an empty file;
- **option matrix 30/30 identical** (`cuda/tools/optmatrix_galn.sh`, on both the 0.2.4 and 0.2.7
  builds): bwa defaults (`-n 0.04 -o 1 -l 32`), `-l 1024`, `-l 16500`, `-l 20 -k 1`, `-n 1 -o 1`,
  `-n 2 -o 0`, `-n 3`, `-n 4`, `-n 0.2`, `-M 1 -O 2 -E 1`, `-e 3`, `-L`, `-m 1000`, `-i 0 -d 100`,
  `-R 1`, each on sub100k (0.5% mapped) and a 100 k slice of shuffled p17 (94% mapped);
- `--verify` (extended to re-run every filtered read on the CPU): **1,143,091 filtered reads
  re-run, 0 with an alignment** (89,203 sub100k + 55,452 p17 + 998,436 p1).

Cost: the worst case is a library in which every read maps -- +0.6 s per ~944 k reads (+10%) on the
3090. Break-even is roughly 0.3% unmapped reads (per read: ~0.6 us of filter on a mapped read
against ~240 us saved on an unmapped one), so for aDNA it is always on the right side.

Not covered: K >= 5 (reads >= 64 bp at `-n 0.01`) and reads > 127 bp go through G-ALN unchanged.
gpualn covers them with its exact DFS existence engine (`dfs_engine.cuh`); porting it is possible but
is a larger piece of code, and on sub100k the long reads cost G-ALN only ~0.4 s.

## Found along the way

- **G-ALN drops the index from the page cache** (since 0.2.6, commit `8990221`; reported as
  [#2](https://github.com/bilalbioinfo/galn/issues/2)). `load_bwt` ends with
  `posix_fadvise(fd, 0, 0, POSIX_FADV_DONTNEED)` ("no page-cache copy"), presumably for GH200 memory.
  So every run re-reads the 3.1 GB `.bwt` from disk (a SATA SSD here), and so does whatever runs next
  on the same index (bwa, `galn samse`/`sampe` via `load_bwt_sa`, gpualn). Controlled check
  (`cuda/tools/galn_pagecache.sh`, sub2k, 0.2.7): load 0.22 s on a warmed cache, then 5.81 and
  5.80 s, with `fincore` showing 0 B of the `.bwt` cached after every run. The same binary with
  `posix_fadvise` made a no-op (`LD_PRELOAD`, `cuda/tools/nofadvise.c`) loads in 0.20 s every run and
  the index stays cached; identical `.sai` either way. `bwa aln -t 32` on sub2k: 0.6 s on a warm
  cache, 6.3 s right after a `galn aln` run.
- **`shuf --random-source=<(yes 42)` does not shuffle.** A constant random source left
  `p17_rand1m.fq` largely in coordinate order (first tenth 100% mapped, last tenth 84.7%), which
  hid every unmapped read from a 100 k head slice. The `*_rand1m.fq` files are kept but are NOT
  shuffled; use `*_shuf1m.fq` (Python `random.Random(42).shuffle`). Correctness results were
  unaffected (every tool read the same file); timings were redone on the shuffled files.
- G-ALN's newer branches (`dev` = 0.2.6, `dev-0.2.7`) are ahead of `main` (0.2.4): GPU index
  building for large genomes, references up to ~137 Gbp (36 -> 38-bit stack fields), NUMA memory
  binding for GH200, `samse -k`. Nothing that moves these numbers: 0.2.7 measured 1-6% slower than
  0.2.4 here with a warm cache.

## Decision

- **The prefilter belongs in G-ALN** -- it is the published, maintained tool for this community,
  and the combination beats both. Submitted 2026-11-01 as
  [#3](https://github.com/bilalbioinfo/galn/pull/3) against their `dev-0.2.7` (fork `teepean/galn`,
  branch `scheme-prefilter-027`, commit `b43ed28` = `cuda/patches/galn-0.2.7-scheme-prefilter.patch`).
- **Pipeline short arm:** gpualn stays until the prefilter is in a G-ALN release (gpualn is faster
  than stock G-ALN on every library we tried); then switch to G-ALN, after a pipeline-level check.
