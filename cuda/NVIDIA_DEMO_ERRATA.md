# Errata: Nsight Copilot **demo** answers about sm_86 (Aug 2026)

Context: the `nvidia-cuda-docs` MCP was returning failures, so four questions were put through the
web demo instead. The demo answered fluently and **most of the load-bearing claims were false**.
This file records them with verbatim quotes and reproducible refutations.

Scope note: this is about the **demo's generated prose**, not the MCP. The MCP retrieval backend
was separately broken (NVIDIA has since fixed it) and, once working, answered all of these
correctly -- see the last section.

## Reproducer

```cuda
// sm86_check.cu
__global__ void k(const int *g, int *out) {
    __shared__ int s[32];
    unsigned sptr = (unsigned)__cvta_generic_to_shared(&s[threadIdx.x]);
    asm volatile("cp.async.ca.shared.global [%0], [%1], 4;\n" :: "r"(sptr), "l"(g + threadIdx.x));
    asm volatile("cp.async.commit_group;\n" ::);
    asm volatile("cp.async.wait_group 0;\n" ::);
    __syncthreads();
    asm volatile("prefetch.global.L2 [%0];" :: "l"(g + 64));
    int v; asm volatile("ld.global.cs.u32 %0, [%1];" : "=r"(v) : "l"(g + 128));
    out[threadIdx.x] = s[threadIdx.x] + v;
}
```

```
$ nvcc -arch=sm_86 -cubin -o sm86_check.cubin sm86_check.cu   # CUDA 13.3, V13.3.73
exit=0            (no errors, no warnings)
$ cuobjdump -sass sm86_check.cubin
    LDGSTS.E   [R9], [R2.64]        <- cp.async
    LDGDEPBAR
    CCTL.E.PF2 [R6]                 <- prefetch.global.L2
    LDG.E.EF   R4, [R4.64+0x200]    <- ld.global.cs  (EF = evict-first)
```

## The false claims

### 1. `cp.async` declared unavailable on Ampere (stated 3x across answers)

> "**`cp.async` (async copy to shared memory)** | ❌ (Hopper-only, SM 90+) | Not present on SM 86, so
> you cannot issue a non-blocking global-to-shared copy."

> "All Hopper-only instructions (cp.async, TMA, prefetch.global.L2, cluster launch, PDL) will be
> stripped out or become no-ops when you compile for sm_86."

False. The Programming Guide states LDGSTS is **CC 8.0+**; the PTX ISA entry for `cp.async` says
"Target ISA Notes: Requires `sm_80` or higher." It compiles on sm_86 and emits `LDGSTS.E`.
(TMA/`cp.async.bulk` and cluster launch really are sm_90+ -- those parts are correct.)

This is the costliest error: it steered us away from a legitimate optimisation for this kernel.

### 2. `prefetch.global.L2` declared unavailable -- contradicting the same demo's earlier answer

> "prefetch.global.L2 PTX instruction | CUDA 12 (SM 90) | No -- the instruction is not recognized on
> SM 86; the assembler drops it."

> "PDL relies on the new 'prefetch-to-L2' instruction (prefetch.global.L2) that is only decoded on
> SM 90."

False. PTX ISA: `prefetch`/`prefetchu` **introduced in PTX ISA 2.0**, "require `sm_20` or higher";
only the `.level::eviction_priority` qualifier needs sm_80. It compiles and emits `CCTL.E.PF2`.

Note the **self-contradiction**: an earlier answer in the same session listed
"`__prefetch_global_l2` / `prefetch.global.L2` | ✅" as available on SM 86 with a usage sketch.

### 3. Fabricated citation -- a "32-entry MSHR" limit attributed to the Programming Guide

> "Load-queue / MSHR pool | 32 entries (**publicly confirmed in the CUDA C++ Programming Guide**
> 'Maximum number of in-flight memory requests per SM')."

No such section, table row, or figure exists. The Programming Guide's per-compute-capability tables
document resident warps/blocks/threads, registers, shared memory and cache working sets -- there is
no published in-flight-request or MSHR depth. The quoted section title appears to be invented.

The number was then used as the foundation of an entire performance theory ("the hard ceiling is the
32-entry load-queue per SM") and of concrete advice.

### 4. Fabricated citations -- `ld.global.cs` attributed to Hopper

> "PTX ISA 7.0, Memory-access instructions -> ld.global.cs entry: notes that the 'streaming' hint is
> only guaranteed to affect eviction priority on SM 90 and later."

> "CUDA Toolkit Release Notes 13.0, New PTX instructions: mentions that ld.global.cs is introduced
> for Hopper; on earlier SMs the instruction is accepted but the hint is ignored."

False, twice. Cache operators including `.cs` were **introduced in PTX ISA version 2.0** and
"require a target architecture of `sm_20` or higher". Neither cited passage exists. And the hint is
demonstrably *not* ignored on Ampere: sm_86 emits `LDG.E.EF` -- an explicit evict-first modifier --
which also refutes the companion claim that "eviction priority is a no-op on Ampere".

### 5. A verification procedure that does not work

> "Inspect the PTX ISA: run `nvdisasm -ptx <cubin>` on a kernel that uses the intrinsic. If the PTX
> contains cp.async or prefetch.global.L2, the assembler will emit an error for SM 86
> (`error: instruction not supported on target`)."

> "Compile with -arch=sm_86 and look for warnings like `'cp.async' is not supported on this
> architecture; instruction will be ignored`."

Neither diagnostic exists. `nvcc -arch=sm_86` compiles all three instructions with exit 0 and no
warnings. A user following this procedure would conclude the opposite of what the tool reports --
which is worse than no procedure, because it manufactures false confidence.

### 6. Advice that would have made the kernel slower

> "The only way to improve performance is to reduce the number of simultaneous load requests per SM,
> **not to 'hide' them with more warps**."

> "keep the number of simultaneously issuing warps at or below the 8-warp sweet spot"

Derived from the invented 32-entry limit. For this workload -- a latency-bound random gather over a
3.1 GB FM-index -- occupancy is the dominant lever, and the shipped speedups came from *raising* the
number of resident warps (removing local memory, then a ptxas schedule that cut registers 95 -> 86).
Following this advice would have inverted the optimisation strategy.

## Did the MCP repeat any of this? No.

After NVIDIA fixed the retrieval backend, `search_cuda_docs` was asked about each disputed point.
It returned correct source text every time:

| disputed point | demo | `search_cuda_docs` |
| --- | --- | --- |
| `cp.async` / LDGSTS | Hopper-only | "LDGSTS (CC 8.0+)"; "Requires `sm_80` or higher" |
| `prefetch.global.L2` | sm_90 only, dropped on sm_86 | "`prefetch` and `prefetchu` require `sm_20` or higher" |
| 32-entry MSHR "in the Programming Guide" | stated as fact | returns the real capability tables; no such figure |
| `ld.global.cs` "introduced for Hopper" | stated twice | "PTX ISA version 2.0 introduced optional cache operators ... require `sm_20` or higher" |

Structurally the MCP is *less* prone to this failure: it returns documentation excerpts rather than
composing prose, so it cannot invent a section title. The fabrications above belong to the demo's
generation layer, not to retrieval. Caveat: five queries is a spot check, not an audit, and semantic
search does sometimes append tangential chunks (a `prefetch` query returned a block of `cvt` notes).
