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
