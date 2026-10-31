# G-ALN: an independent GPU `bwa aln`, and what it means for gpualn

**Status: source-level analysis only (2026-10-31). Nothing has been run** -- the GPU was reserved
for another job. The test plan at the end is what would settle the open questions.

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
does a full faithful search on every no-hit read. **Expected, unmeasured:** gpualn wins on
low-endogenous libraries, G-ALN wins on high-endogenous ones.

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

## Test plan, for when the GPU is free

G-ALN reads the standard bwa index files, so the existing hs37d5 index works as-is.

1. **Byte-identity on our goldens.** `galn aln -dt 16 -l 1024 -n 0.01 -o 2` on sub2k / sub10k /
   sub100k; md5 against `4b068014` / `54410c75` / `eecf35c1`. Also `--gpu-only`.
2. **The two regimes, same machine, same core count:**
   - low-endogenous: sub100k (~0.5% hit) -- gpualn's regime;
   - high-endogenous: the 1 M hit-heavy TAN004/Atacama set from Phases 16/20 -- G-ALN's regime.
   Time gpualn, `galn aln` default, `--gpu-only` and `--cpu-only`, so GPU and CPU contributions
   separate -- the decomposition the preprint does not give.
3. **CPU only, no GPU needed:** `galn aln --cpu-only` returns before any CUDA call (no context is
   created), so its interleaved search can be timed against stock `bwa aln` and our Phase 20
   batch without touching the GPU -- though timings would be noisy while another job loads the box.

Decision it feeds: if G-ALN beats gpualn on high-endogenous data, the pipeline's short arm should
use it there (or gpualn should adopt faithful GPU generation); if gpualn still wins on
low-endogenous data, the routing can depend on expected endogenous content.
