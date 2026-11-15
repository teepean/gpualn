# Measured accuracy: bwa aln vs bwa mem vs the length-split hybrid

First time this project has measured **accuracy against known truth** rather than inferring it from
raw mapped counts or citing Oliva et al.'s published numbers. 300,000 simulated damaged aDNA reads,
20-169 bp in 15 equal bins, drawn from chr22 carrying NA19471 (LWK) SNVs, La Brana damage profile,
mapped to hs37d5. Harness and all caveats: `cuda/tools/benchmark/README.md`.

Definitions: a read is **correct** when it is mapped to the true chromosome with
`|POS - true| == 0`. **sensitivity** = mapped/total, **precision** = correct/mapped,
**accuracy** = correct/total.

## Headline

| strategy | mapped % | precision % | **accuracy %** |
| --- | --- | --- | --- |
| BWA1 `bwa aln -l 1024 -n 0.01 -o 2` | 99.966 | 93.272 | **93.241** |
| **HYBRID c=64** (aln <64, mem9 >=64) | 99.999 | 93.252 | **93.251** |
| BWA9 `bwa mem -k 19 -r 2.5 -L 15` | 92.885 | 94.536 | **87.810** |
| `bwa mem -k 19 -r 2.5` (BWA9 as published) | 92.885 | 94.384 | **87.669** |
| `bwa mem` (default) | 92.885 | 94.389 | **87.674** |

After `MAPQ >= 25`:

| strategy | mapped % | precision % | **accuracy %** |
| --- | --- | --- | --- |
| BWA1 aln | 86.601 | 99.990 | **86.593** |
| HYBRID c=64 | 85.833 | 99.983 | **85.818** |
| BWA9 mem | 82.021 | 99.944 | **81.975** |
| mem default | 81.707 | 99.811 | **81.552** |

**The hybrid is as accurate as all-aln unfiltered (93.251 vs 93.241) and 5.4 points more accurate
than mem alone**, while costing a fraction of all-aln's runtime on real libraries.

Note the trap this measurement springs: **mem has HIGHER precision than aln (94.54 vs 93.27) and
LOWER accuracy by 5.4 points**, because it refuses to map 7% of the reads -- essentially the whole
20-29 bp bin. Precision alone, or mapped-count alone, ranks these backwards. This is the concrete
form of the warning in [[gpu-short-long-cutoff]] that raw mapped-count "can't tell true from false".

## Accuracy by read length (MAPQ >= 0)

| bin | BWA1 aln | BWA9 mem | mem default |
| --- | --- | --- | --- |
| 20-29 | **75.975** | 0.000 | 0.000 |
| 30-39 | **83.855** | 78.145 | 78.080 |
| 40-49 | **88.400** | 88.280 | 88.105 |
| 50-59 | **91.670** | 91.575 | 91.485 |
| 60-69 | 93.775 | **93.850** | 93.725 |
| 70-79 | **94.865** | 94.775 | 94.610 |
| 80-89 | 95.565 | **95.605** | 95.445 |
| 90-99 | **96.165** | 96.085 | 95.880 |
| 100-109 | **96.340** | 96.170 | 95.990 |
| 110-119 | 96.655 | **96.670** | 96.500 |
| 120-129 | 96.805 | **96.870** | 96.710 |
| 130-139 | **97.055** | 97.030 | 96.915 |
| 140-149 | 96.940 | **97.205** | 97.055 |
| 150-159 | 97.200 | **97.415** | 97.265 |
| 160-169 | 97.350 | **97.470** | 97.345 |

`bwa mem` maps **0%** of 20-29 bp reads -- the structural reason BWA-backtrack exists for aDNA. The
aln advantage is **+76.0 points at 20-29 bp and +5.7 at 30-39**, then collapses to <=0.12 points
from 40 bp up, and mem edges ahead above ~110 bp.

## The cutoff question, finally answered with truth

`hybrid_curve.py` gives the exact hybrid accuracy at every cutoff (one aln run + one mem run;
validated against a real split run, see the harness README).

| cutoff | accuracy %, MAPQ>=0 | accuracy %, MAPQ>=25 |
| --- | --- | --- |
| 0 (all mem) | 87.810 | 81.975 |
| 30 | 92.875 | 85.762 |
| 40 | 93.255 | 85.881 |
| 50 | 93.263 | 85.688 |
| 60 | 93.270 | 85.784 |
| **64 (production)** | **93.266** | **85.818** |
| 70 | 93.265 | 85.865 |
| 75 | 93.266 | 85.906 |
| 100 | 93.273 | 86.059 |
| 170 (all aln) | 93.241 | 86.593 |

Two regimes, and they disagree:

- **No MAPQ filter**: the curve is FLAT from 40 to 100 -- every cutoff in that range is within
  **0.01 points**. The cutoff genuinely does not matter. Only the decision to keep <40 bp on aln
  matters, and every sensible cutoff already does that.
- **MAPQ >= 25**: accuracy rises monotonically with the cutoff, and all-aln wins by 0.78 points
  over c=64. Raising the cutoff *helps* here.

Why the difference: at MAPQ>=25 aln both **retains more reads** than mem in almost every bin AND
holds precision >=99.98% in every bin, against mem's 99.67% at 30-39 and 99.76% at 40-49. So aln's
advantage under filtering is real calibration, not overconfidence -- checked explicitly, because an
accuracy gain from a miscalibrated MAPQ would have been an artifact.

**Practical reading.** `adna_aligner.sh` applies no MAPQ filter itself, so for the pipeline as it
stands, c=64 is already optimal to within 0.01 points and [[gpu-short-long-cutoff]]'s "keep auto/64"
conclusion is **confirmed on ground truth** -- it just turns out to be right for a different reason
than precision: the cutoff is simply irrelevant in that range. If a downstream analysis filters
MAPQ (most aDNA work does, at 25 or 30), a higher cutoff buys up to +0.78 accuracy points, paid for
in aln runtime on 64-169 bp reads.

## What reproduced and what did not

- **Reproduced**: well-tuned `bwa aln` is the more accurate choice for aDNA read lengths; the
  advantage is concentrated in short reads; MAPQ filtering sharply improves precision for every
  mapper (94.86% -> 99.98% on the clean validation set).
- **Reproduced, direction only**: their follow-up claim that BWA9 (`-k 19 -r 2.5`) is *worse* than
  default mem below 70 bp. With `-L 15` removed to match their definition, BWA9 is worse or equal in
  every bin below 70 bp -- but by <=0.035 points, far less than the paper implies. Our pipeline's
  extra `-L 15` is worth **+0.14 accuracy points** over published BWA9, so keep it.
- **Not reproducible at all**: their code. All three published scripts fail to parse or run; the VCF
  they name has been superseded. Details in the harness README.
- **Not tested here**: reference bias at SNPs, D-statistics, PCA -- the other halves of their paper.

## Free side-benefit

`gpualn` was md5-identical to CPU `bwa aln` on this dataset (`043cb47f...`), which is an independent
bit-exactness check on a read population nothing like the aDNA libraries used so far: simulated,
damaged, uniformly spread over 20-169 bp.

# Reference bias at heterozygous SNPs

The metric that matters most for aDNA downstream work. 1.5 M simulated reads (the read count Oliva
et al. used), chr22 + NA19471, La Brana damage, hs37d5. The simulator draws each read's haplotype
50/50, so at a het site the **expected alt-allele fraction is exactly 0.500**; anything below is
reference bias. `cuda/tools/benchmark/ref_bias.py`.

Transitions (A<->G, C<->T) are reported separately from transversions because C>T / G>A damage
masquerades as a real allele at transition sites and cannot at transversion sites. **Transversions
are the damage-immune read on mapper bias** and all headline numbers below use them.

## The control comes first

A "perfect mapper" measurement, counting alleles at the sites each read TRULY spans using the
encoded truth, isolates bias from damage and sequencing error alone:

| | ref | alt | alt_frac | bias |
| --- | --- | --- | --- | --- |
| transversion | 23,093 | 23,057 | **0.49961** | +0.00039 |
| transition | 55,416 | 55,226 | 0.49914 | +0.00086 |

0.49961 against an expected 0.5, n = 46,150, SE 0.0023. The simulation carries no built-in bias, so
any deviation in the alignments is attributable to the mapper.

## Measured bias (transversions)

| strategy | no filter | MAPQ>=25 | MAPQ>=25, trim 5 |
| --- | --- | --- | --- |
| control (perfect mapper) | +0.00039 | — | — |
| BWA1 `aln` | +0.01312 | +0.00952 | +0.00857 |
| **HYBRID c=64** | **+0.01196** | **+0.00889** | **+0.00735** |
| `mem -k19 -r2.5 -L15` | +0.01482 | +0.01181 | +0.01047 |
| `mem -L15` (no -k/-r) | — | +0.01188 | — |
| `mem` default | — | +0.01342 | — |

**Reproduces Oliva et al.'s MAPQ finding**: filtering reduces bias for every mapper (aln 1.31 ->
0.95, mem 1.48 -> 1.18, hybrid 1.20 -> 0.89 points). End-trimming 5 bp reduces it further.

## Significance: paired bootstrap, not eyeballing

Differences here are 0.1-0.3 points against an unpaired SE of ~0.24 points, so the ordering cannot
be read off the table. All strategies see the same reads, so the comparison is paired and the het
SITE is the resampling unit (`bootstrap_bias.py`, B = 2000 over 10,804 sites):

| config | hybrid - aln | hybrid - mem9 | aln - mem9 |
| --- | --- | --- | --- |
| no filter | **-0.00116** sig | **-0.00287** sig | **-0.00171** sig |
| MAPQ>=25 | -0.00063 **ns** | **-0.00292** sig | **-0.00228** sig |
| MAPQ>=25 trim5 | **-0.00123** sig | **-0.00312** sig | **-0.00189** sig |

- **aln has significantly less reference bias than mem in all three configurations.**
- **The hybrid has significantly less bias than mem in all three.**
- The hybrid's edge over all-aln is small (~0.1 points) and significant in 2 of 3 configurations --
  worth stating as "at least as good as aln", not as a win. On first inspection the hybrid looked
  lowest everywhere; the bootstrap showed one of those orderings is noise.

## `-k 19 -r 2.5` does nothing; `-L 15` does all the work

Decomposing the mem options (1.5 M reads):

| options | accuracy % | bias (MAPQ>=25) |
| --- | --- | --- |
| `mem` default | 87.764 | +0.01342 |
| `mem -L 15` | **87.891** | **+0.01188** |
| `mem -k 19 -r 2.5 -L 15` | 87.882 | +0.01181 |
| `mem -k 19 -r 2.5` | 87.669 (300k set) | — |

`-k 19` **is bwa mem's own default**, so BWA9's `-k 19` is a no-op; `-r 2.5` (against a default of
1.5) moves accuracy by -0.009 and bias by 0.0001 -- nothing. Adding `-L 15` is worth **+0.13
accuracy points and 0.15 bias points**.

So the DNAharvester-inherited `-k 19 -r 2.5` in `MEM_OPTS` can be dropped with no measurable loss,
and `-L 15` is the part worth keeping. This also explains why Oliva et al.'s follow-up claim about
BWA9 is so weak in our hands: the parameters it turns on barely do anything.

## Where the bias actually lives: 30-49 bp, and it is enormous for mem-family mappers

The aggregate numbers above hide the effect that matters. Transversion bias per length bin,
MAPQ>=25 (`cuda/tools/benchmark/bias_by_len.py`):

| bin | BWA1 aln | mem default | minibwa |
| --- | --- | --- | --- |
| 20-29 | +0.08500 | *unmapped* | *unmapped* |
| **30-39** | **+0.03682** | **+0.20064** | **+0.34600** |
| **40-49** | **+0.00129** | **+0.02982** | **+0.10415** |
| 50-59 | +0.00563 | +0.01699 | +0.02937 |
| 60-69 | +0.01066 | +0.02278 | +0.02649 |
| 70-79 | +0.01601 | +0.02088 | +0.02274 |
| 80-89 | +0.00711 | +0.01333 | +0.01333 |
| 90-99 | -0.00158 | +0.00196 | +0.00286 |
| 100-169 | -0.012 .. +0.021 | -0.008 .. +0.022 | -0.008 .. +0.021 |

At **30-39 bp, `bwa mem` loses 20% of alt alleles and minibwa loses 35%** -- an alt-allele fraction
of 0.299 and 0.154 against a true 0.500. `bwa aln` loses 3.7%. By 90 bp all three converge and the
residual is noise (it goes slightly negative in places).

**Mechanism.** mem and minibwa require a minimum alignment score. A read carrying the non-reference
allele has at least one mismatch, which on a 30-40 bp read is enough to push it under the floor --
so the alt-bearing copy is dropped while the ref-bearing copy is kept. minibwa is worse because its
floor is higher relative to short reads (`-s 30` x `-A 2` = score 60, and a perfect 30 bp read
scores exactly 60). Supporting evidence: minibwa maps fewer short reads than mem (92.1% vs 93.8% at
30-39 bp). `bwa aln`'s bounded-edit-distance search has no score floor, so it keeps them.

A soft-clipping explanation was proposed first and **refuted by measurement**: minibwa soft-clips
0.42% of reads against mem default's 0.39% -- indistinguishable. (`-L 15` does suppress clipping as
intended, 0.39% -> 0.01%, which is consistent with it reducing mem's bias.)

**This is why the hybrid works.** The catastrophic bins are exactly the ones a length split routes
to aln, and it is why the hybrid tracks all-aln rather than mem. It also means the choice of
long-read mapper barely matters for bias -- above ~90 bp every mapper tested is within noise -- while
the choice of SHORT-read mapper matters enormously. Any aDNA pipeline that sends 30-50 bp reads to
a mem-family aligner is discarding roughly a fifth to a third of its non-reference alleles in that
band.

# Alternative long-arm mappers: minibwa and strobealign

Two newer short-read mappers tested as drop-in replacements for `bwa mem`, both standalone and as
the hybrid's long arm. Same 1.5 M reads, same reference, defaults as a mem replacement (neither has
published aDNA settings).

- **minibwa** 2606.15357 (lh3), the declared successor to bwa-mem: BWT seeding plus minimap2-style
  chaining and SIMD extension. Needs its own index (`minibwa index`, 6.3 GB `.mbw` + 785 MB `.l2b`).
- **strobealign** 0.19.0 (ksahlin), randstrobe seeding, now rewritten in Rust. Its README states it
  is "most suited for read lengths between 100 and 500 bp" -- most aDNA reads are outside that range,
  which is a caveat on its results here, not a defect.

## Standalone (all 20-169 bp reads)

| strategy | mapped % | precision % | accuracy % | bias (tv, MAPQ>=25) |
| --- | --- | --- | --- | --- |
| BWA1 `aln` | 99.962 | 93.370 | **93.334** | **+0.00952** |
| `mem -k19 -r2.5 -L15` | 92.912 | 94.586 | 87.882 | +0.01181 |
| `mem` default | 92.913 | 94.459 | 87.764 | +0.01342 |
| minibwa | 92.724 | 94.429 | 87.558 | +0.01711 |
| strobealign | 97.543 | 93.383 | 91.088 | +0.01953 |

strobealign **maps far more reads than mem (97.5% vs 92.9%) and is 3.2 accuracy points better**,
because it manages short reads that mem drops outright. But it has the **worst reference bias of
everything tested**.

## As the hybrid's long arm (aln below 64 bp, X at or above)

| hybrid | accuracy % | bias (tv, MAPQ>=25) |
| --- | --- | --- |
| **aln + `mem -k19 -r2.5 -L15`** | **93.345** | **+0.00889** |
| aln + minibwa | 93.243 | +0.01075 |
| aln + strobealign | 93.101 | +0.01394 |

Paired bootstrap over het sites, all differences **significant**:
`aln+mem9 - aln+minibwa = -0.00186` [-0.00229, -0.00141];
`aln+mem9 - aln+strobe = -0.00505` [-0.00626, -0.00376];
`aln+minibwa - aln+strobe = -0.00319` [-0.00444, -0.00195].

**Verdict: `bwa mem` with `-L 15` remains the best long arm. Neither newer mapper improves on it,
on accuracy or on bias.**

## The non-obvious part

**Standalone ranking does not predict hybrid-arm ranking.** strobealign beats mem by 3.2 accuracy
points standalone, yet makes the hybrid 0.24 points *worse* than mem does. The reason: strobealign's
entire standalone advantage is in short reads mem cannot map -- and in the hybrid those reads go to
`aln` anyway, which handles them better than either. What remains is strobealign's weaker long-read
behaviour and its higher bias. Evaluating a candidate long arm on whole-library numbers would pick
the wrong one.

## Speed, with a caveat

1.5 M reads, 16 threads: `bwa mem` ~2.7 s, minibwa 2.57 s, strobealign 20.7 s -- but strobealign's
figure includes building its index on the fly on every run, so it is not a mapping-only time, and
none of the three is in its design regime on 20-169 bp reads. Treat these as "no mapper here is a
throughput problem", not as a benchmark. minibwa's claimed 3x over bwa-mem does not appear on this
workload.

## minibwa has no clipping-penalty option

mem's `-L` is the single option that mattered (+0.13 accuracy, 0.15 bias points). minibwa's
alignment options are `-A -B -O -E -s` with **no clipping penalty at all**, so that lever is simply
unavailable there. Worth knowing before adopting it for aDNA.

## Validated on bwa-mem3, and applied to the production pipeline

The benchmark above used `bwa mem` 0.7.18-r1243 -- `/usr/bin/bwa`, what the `-a mem` path runs --
but `adna_aligner.sh -a mem3` runs **bwa-mem3** with an extra flag (`--supp-rep-hard-cap 5`).
Re-ran the decomposition on the actual production binary and flag, same 1.5 M reads:

| variant | mapped % | precision % | accuracy % | bias (tv, MAPQ>=25) |
| --- | --- | --- | --- | --- |
| mem3 `-k19 -r2.5 -L15 --supp-rep-hard-cap 5` | 92.912 | 94.586 | 87.882 | +0.01181 |
| mem3 `-L15 --supp-rep-hard-cap 5` | 92.913 | 94.595 | **87.891** | +0.01188 |
| bwa mem `-k19 -r2.5 -L15` (reference) | 92.912 | 94.586 | 87.882 | +0.01181 |

bwa-mem3 reproduces bwa mem to three decimals, and `--supp-rep-hard-cap 5` changes nothing on
single-end data. Dropping `-k 19 -r 2.5` gains +0.009 accuracy points; the bias difference is
**-0.00008, 95% CI [-0.00034, +0.00018] -- not significant** by paired bootstrap. So the change is
accuracy-neutral-to-marginally-positive and bias-neutral.

**Applied** to `/home/dnastorage/aDNApipeline/adna_aligner.sh` on 2026-10-14:
`MEM_OPTS="-k 19 -r 2.5 -L 15"` -> `MEM_OPTS="-L 15"`, with the two help-text defaults and the now
-stale fq2bam comments updated to match (fq2bam's `--bwa-options` never accepted `-k`/`-r`, so it
no longer has to drop anything). Timestamped backup alongside the script.

## bwa-mem3 v0.14.1: 2x faster, byte-identical output

The production binary was `0.2.1-a24de32` (built 2026-05-20). Upstream has moved 407 commits and 15
releases since, so the whole decomposition was re-run on **v0.14.1** (2026-10-11). The question is
not "is the new one faster" -- it is **"does any of that development change an alignment"**, because
accuracy is the deciding factor here and a 2x speedup that moves a single base is not a speedup.

### Accuracy and bias, same 1.5 M reads, same options

| variant | mapped % | precision % | accuracy % | bias (tv, MAPQ>=25) |
| --- | --- | --- | --- | --- |
| 0.2.1 `-L 15 --supp-rep-hard-cap 5` | 92.913 | 94.595 | **87.891** | +0.01188 |
| **v0.14.1 `-L 15 --supp-rep-hard-cap 5`** | 92.913 | 94.595 | **87.891** | +0.01188 |
| v0.14.1 `--compat=bwa-mem -L 15` | 92.913 | 94.595 | **87.891** | +0.01188 |
| v0.14.1 `--fast -L 15 --supp-rep-hard-cap 5` | 92.899 | 94.398 | 87.694 | +0.01487 |

### Byte-identity, verified on three independent inputs

`cmp` of the alignment records (`@PG` excluded -- it embeds argv) between 0.2.1 and v0.14.1 at the
production option set:

| input | reads | mapped | result |
| --- | --- | --- | --- |
| `sim15M.fq` simulated damaged aDNA, 20-169 bp | 1,500,000 | 92.9% | **byte-identical** |
| `test_data/sub100k.fq` real aDNA (atacama soil), 30-91 bp | 100,000 | 0.5% | **byte-identical** |
| `HRR1274614_f1` real human, 150 bp | 200,000 | 99.9% | **byte-identical** |

Three regimes on purpose: the simulated set for the damage/length profile, the soil sample because
almost everything is unmapped there, and the 150 bp human set because almost everything maps and it
is the only one of the three that produces supplementary records.

`--compat=bwa-mem` is also exactly what it claims: its records are **byte-identical to
`/usr/bin/bwa mem -L 15`** over all 1.5 M reads. The default path differs from `--compat` only by an
added `HN:i:` tag (1,393,697 records); POS, MAPQ, CIGAR and NM are untouched, which is why all three
rows above score the same. v0.12.0's chain-weight-filter change (#489) moved the *default* path onto
bwa-mem's "leave it unmapped" behaviour, so `--compat=bwa-mem` is now a no-op for mapping on this
workload and only `--compat=bwa-mem2` resurrects those alignments.

### Speed, 1.5 M reads, 16 threads, best of 3 (all within 1.5%)

| binary | wall | index load | mapping only |
| --- | --- | --- | --- |
| `bwa` 0.7.18-r1243 | 13.33 s | 0.40 s | 12.93 s |
| bwa-mem3 0.2.1 | 6.01 s | 2.20 s | 3.84 s |
| **bwa-mem3 v0.14.1** | **3.02 s** | **0.72 s** | **2.31 s** |
| bwa-mem3 v0.14.1 `--fast` | 2.16 s | 0.72 s | 1.44 s |

**2.0x end-to-end over 0.2.1 and 4.4x over `bwa mem`, at byte-identical output.** Note where a third
of the gain comes from: **index load dropped 2.20 s -> 0.72 s** (3.1x). That is per invocation, and
`adna_aligner.sh` invokes the aligner once per sample, so it is a real saving on short libraries
rather than a benchmark artifact.

### `--fast` is rejected

v0.12.0 added a `--fast` preset (`--max-extend-chains 20 --adaptive-band --extend-tie-frac 0.95`
...), honestly documented upstream as NOT byte-identical. It costs **0.197 accuracy points** and
**+0.00299 reference bias**, the latter significant by paired bootstrap over 10,904 het sites:
`v0.14.1 - v0.14.1 --fast = -0.00299`, 95% CI `[-0.00391, -0.00207]`. For 1.4x it gives back a fifth
of the 1.0-point margin this project spent the whole benchmark establishing. Not worth it here.
Useful to know it is available for a throughput-bound, non-population-genetic workload.

### The upgrade is drop-in

- **No index rebuild needed** -- v0.14.1 reads the existing hs37d5 index as-is (verified).
- `index --max-memory` (the flag `adna_aligner.sh:217` passes) still exists.
- Both index call sites in the script stay correct: it builds the classic `bwa` index
  unconditionally for the `bwa aln` short arm, and the mem3 index only for `-a mem3`.
- `BWA_MEM3` is already an overridable env var (`adna_aligner.sh:51`), so switching is one path.
- One difference if the mem3 index is ever *rebuilt*: v0.14.1 no longer writes `.0123` by default
  (`--emit-unpacked-ref` is now opt-in) and does not need it to map. That saves 5.8 GB for hs37d5.
  Nothing in the pipeline reads `.0123` -- fq2bam uses the classic `bwa` index -- so this is safe,
  but do not "fix" a missing `.0123` after a rebuild.

**Applied** to `/home/dnastorage/aDNApipeline/adna_aligner.sh` on 2026-10-14: `BWA_MEM3` now
defaults to `/home/teemu/sorsa/bwa-mem3-v0141/bwa-mem3`, with the measurement and the "do not add
`--fast`" warning recorded at the assignment. The stale `# bwa-mem3 with the local hard-cap fix`
comment at the `mem3)` call site was corrected -- the cap is upstream (PRs #100/#101) and was
re-verified to behave identically on v0.14.1 (119 -> 152 supplementary records forced to MAPQ 0 on
the 200 k human set, same counts as 0.2.1). Timestamped backup alongside the script.

Two other scripts in that directory still hardcode the 0.2.1 path and were left alone:
`paired_mem3.sh:14` (also single-end -- collapsed reads, no `-p` -- so the same verification covers
it; note it runs `-k 19 -r 2.5` and no `-L 15`, an option set that predates the MEM_OPTS work) and
`test_supp_rep_hard_cap.sh:24` (a test harness, where pinning the old binary is arguably correct).

# strobealign-HRY: tested 2026-11-15, not adopted

Tolstoganov & Sahlin, *Hybrid RY-mers increase aDNA alignment sensitivity* (bioRxiv
10.64898/2026.11.05.764514, v2 2026-11-10) adds an aDNA mode to strobealign. Its seeds are k-mers
whose first `--ry-len` (8) positions are hashed by purine/pyrimidine class only, so C>T / G>A damage
near read ends no longer breaks the seed, while the rest of the k-mer stays nucleotide-exact. It
claims higher accuracy than `bwa aln -n 0.01 -o 2 -l 1024` (the same options we use) at 4-6x its
speed. Tested here as published: tag `hry` (`4be2dcc`, strobealign 0.18.0-alpha, Rust, `cargo build
--release`) with the authors' settings from their benchmark repo (`adna_bench`,
`simx/config/config.yaml`): `--adna --ry-len 8 -M 500 --ssw -k 17 -s 15`.

## How their benchmark differs from this one

| | the paper | here |
| --- | --- | --- |
| reads | 40 bp only (gargammel), 80% endogenous | 20-169 bp in 15 bins, all endogenous (+ real libraries for speed) |
| divergence | 0.5 / 2 / 5% Mason substitutions + indels written into one haploid copy of fruit fly, Zymo, CHM13 | real NA19471 SNVs on chr22 (heterozygous and homozygous), mapped to hs37d5 |
| "correct" | primary overlaps the true interval, and alignment score >= 60 under strobealign's scoring (bwa's recomputed from MD/CIGAR) | unclipped start == true start (a 100 bp tolerance changes accuracy by 0.03 points) |
| MAPQ-filtered results | none | all tables at MAPQ>=0 and >=25, plus precision-matched |
| reference bias | **not measured** -- their design has no heterozygous sites, so there is no allele balance to check | transversion bias at 10,811 het sites, paired bootstrap |

Also seen in `adna_bench`: its `compile_strobealign` rule builds with plain `cargo build` (dev
profile, `opt-level = 1`). If the paper's timings came from that rule, they understate HRY's speed;
the speed numbers below use a release build.

## Accuracy (1.5 M simulated reads, unclipped start)

| strategy | mapped % | precision % | accuracy % | accuracy, MAPQ>=25 | best accuracy at precision >= 99.9% / 99.99% |
| --- | --- | --- | --- | --- | --- |
| **BWA1 `aln -l 1024 -n 0.01 -o 2`** | 99.962 | 93.370 | **93.334** | **86.731** | **89.836 / 86.731** |
| bwa-mem3 v0.14.1 `-L 15` (pipeline long arm) | 92.913 | 94.595 | 87.891 | 81.827 | 84.405 / never reaches 99.99 |
| **strobealign-HRY, paper settings** | 98.286 | 94.361 | **92.744** | **80.075** | **81.769 / 72.065** |
| strobealign-DK, paper settings (`-l 1 -u 1 -M 500 -k 17 -s 15`) | 98.704 | 94.081 | 92.862 | 80.246 | |
| strobealign, defaults (same binary) | 97.553 | 93.501 | 91.213 | 81.195 | |
| strobealign-HRY, release default (`--adna` only) | 96.963 | 93.622 | 90.779 | 78.357 | |

HRY with `-r 50` (the short-read profile; the set's estimated length is 94 bp) gave byte-identical
records to the paper-settings run: with `-k`/`-s` fixed, the profile does not matter.

Per length bin, accuracy % (MAPQ>=0 | MAPQ>=25):

| bin | aln | mem3 | HRY | DK | strobealign | HRY - aln |
| --- | --- | --- | --- | --- | --- | --- |
| 20-29 | **75.79** \| **56.70** | 0 \| 0 | 73.61 \| 46.57 | 72.48 \| 42.33 | 53.89 \| 46.21 | -2.18 \| -10.14 |
| 30-39 | **83.96** \| **73.38** | 78.42 \| 69.14 | 83.20 \| 65.35 | 82.82 \| 64.21 | 79.53 \| 65.88 | -0.76 \| -8.03 |
| 40-49 | **88.67** \| 78.85 | 88.43 \| **79.97** | 87.99 \| 72.98 | 88.24 \| 72.50 | 86.63 \| 72.16 | -0.69 \| -5.87 |
| 50-59 | **91.99** \| **83.19** | 91.74 \| 81.49 | 91.19 \| 77.02 | 91.62 \| 77.00 | 90.84 \| 76.16 | -0.79 \| -6.17 |
| 60-69 | **93.87** \| **86.49** | 93.84 \| 85.27 | 92.93 \| 79.62 | 93.49 \| 79.85 | 93.13 \| 79.10 | -0.93 \| -6.88 |
| 70-99 | **95.09-96.17** \| **88.9-91.4** | within 0.1 \| 0.8-1.1 lower | 0.6-0.8 lower \| ~7 lower | | | -0.64..-0.82 \| ~-7.0 |
| 100-169 | 96.5-97.4 \| **92.2-93.9** | +-0.15 \| 0.7-1.4 lower | -0.47..+0.04 \| 5.5-6.6 lower | | | -0.47..+0.04 \| -5.5..-6.6 |

- **HRY is below `bwa aln` in 14 of 15 bins** (the 15th: +0.04 at 150-159 bp). In the paper's own
  read length, 40-49 bp, it is 0.69 points below.
- **At MAPQ>=25 it is 5.5-10 points below in every bin.** This is not MAPQ scale: at equal precision
  aln keeps 89.8% vs HRY's 81.8% (99.9%) and 86.7% vs 72.1% (99.99%). There HRY is below `bwa mem`
  too (84.4% at 99.9%).
- **What the RY seeds do buy**: the shortest reads. At 20-29 bp HRY maps 73.6% correctly, against
  53.9% for default strobealign and 0% for mem; vs the paper's own non-aDNA baseline DK the gain is
  +1.1 points at 20-29 and +0.4 at 30-39, and DK is ahead from 40 bp on (overall 92.86 vs 92.74).
- The paper's one non-win is also human: CHM13 at 2%, where bwa aln is 0.1 points ahead. Its
  lowest divergence is 0.5%, about five times the SNV rate here, and its margins grow with
  divergence -- a regime that matters for non-model and microbial references, not for human aDNA on
  hs37d5.

## Reference bias (transversions, MAPQ>=25; perfect-mapper control +0.00039)

| strategy | bias | per bin 30-39 / 40-49 / 60-69 |
| --- | --- | --- |
| **aln** | **+0.00952** | +0.037 / +0.001 / +0.011 |
| mem3 `-L 15` | +0.01188 | +0.198 / +0.028 / +0.018 |
| **HRY** | **+0.01819** | +0.132 / +0.029 / +0.032 |
| DK | +0.01769 | +0.130 / +0.035 / +0.032 |
| strobealign default | +0.01977 | |

At 20-29 bp (n ~ 220-400) HRY loses 42.5% of alt alleles (aln 8.5%). Paired bootstrap over 10,811
het sites (B = 2000): `aln - HRY = -0.00867 [-0.01048, -0.00685]`, `HRY - mem3 = +0.00631
[+0.00467, +0.00802]`, both significant. **HRY nearly doubles aln's reference bias and has more than
mem**. From 50 to 89 bp it is 2-5x aln's and above mem3's, so it is not only the short-read score
floor that hurts mem-family mappers. Hashing the first 8 seed positions by R/Y class tolerates
transitions only; a transversion SNP breaks the seed exactly as in plain strobealign, and HRY's bias
is close to DK's (+0.0182 vs +0.0177).

## As either arm of the length-split hybrid (cutoff 64)

| hybrid (short arm + long arm) | accuracy % | accuracy, MAPQ>=25 | bias, MAPQ>=25 |
| --- | --- | --- | --- |
| **aln + bwa-mem3 (production)** | **93.347** | **85.962** | **+0.00913** |
| HRY + bwa-mem3 | 93.030 | 83.773 | +0.01397 |
| aln + HRY | 93.061 | 82.264 | +0.01291 |

At cutoff 40, HRY + mem3 is 93.117 / 84.679 vs 93.313 / 85.890. Bootstrap: `production - (HRY +
mem3) = -0.00484 [-0.00556, -0.00414]`, `production - (aln + HRY) = -0.00378 [-0.00534, -0.00236]`,
both significant. **The production hybrid wins whichever arm HRY replaces.** (The production bias
here is +0.00913 rather than the +0.00889 measured with `mem -k 19 -r 2.5 -L 15` before the
2026-10-14 MEM_OPTS change; the long arm is now what the pipeline runs.)

## Real libraries: speed (strictly one run at a time, `-t 32`, page cache warm)

| | sub100k | p17 1M | p1 1M |
| --- | --- | --- | --- |
| reads / mapped by bwa aln | 100,000 / 0.5% | 999,552 / 94.4% | 1,000,528 / 0.2% |
| `bwa aln -t 32` (today; 2026-11-01: 15.3 / 56.1 / 149.1) | 15.40 s | 57.05 s | 138.69 s |
| **strobealign-HRY, saved index** (median of 3) | **9.14 s** | **8.91 s** | **14.38 s** |
| strobealign-HRY, index built in the run | 32.10 s | 28.93 s | 33.03 s |
| gpualn (2026-11-01, GPU held today) | 2.22 s | 17.4 s | 7.77 s |
| G-ALN + scheme prefilter (2026-11-01) | 1.97 s | 9.63 s | 6.95 s |
| peak memory: bwa aln / HRY | 3.7 / 20.9 GB | 6.7 / 21.6 GB | 3.8 / 21.4 GB |

HRY's index for hs37d5 takes 41.7 s and 23.4 GB to build and is 21.4 GB on disk (`.sti`); a saved
run spends 2.2 s reading it and 3.8-9.2 s mapping. Saved-index and built-in-run outputs are
byte-identical on all three libraries. So: **1.7x / 6.4x / 9.6x over CPU bwa aln**, consistent with
the paper's 4-6x; **faster than gpualn on the 94%-endogenous p17 (8.9 vs 17.4 s), but 4.1x slower on
sub100k and 1.9x slower on p1**, where gpualn's filter proves most reads unalignable. Against G-ALN +
prefilter, HRY is slightly faster on p17 (8.9 vs 9.6 s) and slower on the other two.

## Real libraries: HRY reports almost every read as mapped

| | sub100k | p17 | p1 |
| --- | --- | --- | --- |
| bwa aln mapped / MAPQ>=25 | 0.51 / 0.21% | 94.42 / 68.96% | 0.20 / 0.14% |
| HRY mapped / MAPQ>=25 | **99.72** / 0.35% | 95.45 / 63.27% | **99.99** / 0.36% |
| HRY with AS>=60 / + MAPQ>=25 | 0.28 / 0.17% | 90.15 / 63.26% | 0.18 / 0.10% |

strobealign always reports its best local alignment. On p1, 99.8% of HRY's mapped records are
soft-clipped with an alignment score below 60: a 17-20 bp exact hit with the rest of the read
clipped. **MAPQ does not remove them**: of HRY's 3,587 MAPQ>=25 records on p1, 71% have AS < 60
(e.g. `17M28S`, MAPQ 25), and only 1,022 of them are reads bwa aln maps at all. The simulated set
(all endogenous) cannot show this. The paper's evaluation applies AS >= 60, so its numbers never
include these records; a pipeline would have to apply the same filter (`samtools view -q` is not
enough). That filter costs HRY 2.6 accuracy points at 20-39 bp on the simulated set (78.40 ->
75.81%; aln 79.87%), nothing from 40 bp up, and on p17 it leaves HRY with 5.7 points fewer MAPQ>=25
reads than aln (63.26 vs 68.96%), in line with the 5.9-8.0 point gap at 30-49 bp on the simulated
set.

## Verdict

**Not adopted.** On human aDNA against hs37d5, strobealign-HRY is less accurate than `bwa aln` in
every length bin, far less at any MAPQ or precision cut, and has about twice its reference bias. As
the short arm it makes the hybrid worse on both accuracy and bias, and as the long arm too. Its
real-library speed (9-14 s per 1 M reads on 32 cores, 21 GB RAM, 21 GB index) beats CPU bwa aln
6-10x but not gpualn on low-endogenous libraries, and its output needs an alignment-score filter
before it can be counted. It is the strongest seed-based option tested for very short reads
(73.6% at 20-29 bp, where mem maps none), and a reasonable choice where the paper aims it:
divergent or metagenomic references, and CPU-only machines where bwa aln's runtime is the
bottleneck.

## Two scorer fixes made for this test

- `score_sam.py` / `hybrid_curve.py` compared SAM `POS` to the true start, so a read with a leading
  soft clip was always wrong. They now compare the unclipped start (`--raw-pos` restores the old
  rule). aln never clips and mem3 `-L 15` clips 0.007% of reads, so their numbers are unchanged;
  strobealign clips 0.6-0.7% (1.1-1.4% below 64 bp) and gains 0.12 points (default 91.094 ->
  91.213). The earlier verdict on strobealign as a long arm stands.
- `ref_bias.py` / `bias_by_len.py` counted alleles at chr22 coordinates for reads mapped to *any*
  chromosome. They now skip records off chr22 (`--chrom ''` restores the old rule). At MAPQ>=25
  nothing changes to 5 decimals (aln +0.00952, mem3 +0.01188); unfiltered aln moves +0.01312 ->
  +0.01306.

Reproduce: `cuda/tools/benchmark/hry_run.sh` (simulated set; scoring with `score_sam.py`,
`mapq_curve.py`, `hybrid_curve.py`, `merge_hybrid.py`, `ref_bias.py --dump`, `bootstrap_bias.py`,
`bias_by_len.py`, `clipstats.py`) and `hry_speed.sh` (real libraries). The simulated set was
regenerated with seed 7 and reproduces the earlier aln and mem numbers exactly; it now lives in
`/home/teemu/sorsa/adna_simbench/work` (outside the repo, survives reboots).
