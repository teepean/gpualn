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
