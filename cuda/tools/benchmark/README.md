# Accuracy benchmark: aln vs mem vs the length-split hybrid, on simulated aDNA with known truth

Answers a question the production pipeline could not: **which mapper is actually more accurate per
read length, and where should the hybrid cut over?** Real data cannot answer it -- there is no truth,
and raw mapped-count cannot tell a true mapping from a confident wrong one.

Modelled on the design of Oliva et al. 2021 (*Briefings in Bioinformatics*,
[doi:10.1093/bib/bbab076](https://doi.org/10.1093/bib/bbab076)) and their
[Benchmark-aDNA-Mapping](https://github.com/AdrienOliva/Benchmark-aDNA-Mapping) repo, but NOT their
code -- see "Why not their code" below.

## Pipeline

```sh
samtools faidx hs37d5.fa 22 > chr22.fa
bcftools view -s NA19471 -v snps -m2 -M2 -Ou ALL.chr22.v5b.vcf.gz \
  | bcftools query -f '%POS\t%REF\t%ALT\t[%GT]\n' > NA19471.chr22.snv.tsv
curl -sSLO https://raw.githubusercontent.com/grenaud/gargammel/master/examplesMapDamage/results_LaBrana/misincorporation.txt

python3 -I simulate_adna.py --ref chr22.fa --variants NA19471.chr22.snv.tsv \
        --damage misincorporation.txt --n 300000 --out sim300k.fq --seed 1

bwa aln -l 1024 -n 0.01 -o 2 -t 16 $REF sim300k.fq > a.sai
bwa samse $REF a.sai sim300k.fq > sim.aln.sam            # BWA1
bwa mem -k 19 -r 2.5 -L 15 -t 16 $REF sim300k.fq > sim.mem9.sam   # the hybrid's long arm

python3 -I score_sam.py   sim.aln.sam --tol 0 --mapq 25
python3 -I hybrid_curve.py --aln sim.aln.sam --mem sim.mem9.sam --mapq 25
```

`simulate_adna.py` encodes the truth in the read name as
`r<idx>:<chrom>:<1-based leftmost ref pos>:<strand>:<len>:<hap>`. That is directly comparable to
SAM `POS` for an end-to-end alignment on either strand, so no strand- or clipping-dependent
correction is needed.

`hybrid_curve.py` computes the hybrid's accuracy at **every** cutoff from one aln SAM plus one mem
SAM: the pipeline routes each read to exactly one mapper by length, and single-end aln and mem are
both per-read independent. Verified against a real split run -- identical at MAPQ>=25, and 44 reads
in 300,000 apart unfiltered, all MAPQ-0 multi-mappers where bwa's choice among equally-best hits
depends on read order.

## Validation of the harness itself

Simulated a clean set (`--no-damage --qual 60`) and mapped it. Of the wrong alignments at >=120 bp,
266 were MAPQ 0 and only **2** were MAPQ>=30; precision at MAPQ>=30 is 99.95-100%. So the truth
encoding is right and the residual error is genuine multi-mapping in chr22 segmental duplications
and hs37d5 decoys. Do not skip this step if you change the simulator.

## Why not their code

All three published scripts fail to run, verified by parsing only:

| file | failure |
| --- | --- |
| `Analysis/CalculateMappingDistance.R` | **R parse error, line 55**: `reads on foward strand with clipping` is a comment missing its `#`, so the file never parses |
| `Simulations/CreateReads/[Pipeline]-...py` | **SyntaxError, line 40**: tabs mixed with spaces; also calls `get_ipython()` and `/gen_reads.sh` (absolute path at the filesystem root) |
| `gen_reads.sh` | needs Mitty + the `1kg-pcr-free.pkl` read model, and the VCF it names (`...v5a...`) no longer exists -- 1000G now ships `v5b` |

So strict reproduction is impossible regardless of tooling. This harness reproduces their
experimental *design* with runnable code.

## Deliberate deviations

- **SNVs only, no indels.** Indels shift coordinates between haplotype and reference and blur the
  ground truth the benchmark rests on. Reference bias at SNPs is unaffected.
- **No Mitty / gargammel.** Damage comes from the same La Brana mapDamage table they used, applied
  as C>T from the 5' end and G>A from the 3' end (what that table measures). Their command passed
  `single` to `deamSim`; this is the double-stranded reading.
- Base qualities are flat Q37 with a matching error rate, not Mitty's `1kg-pcr-free` model.
- hs37d5 (what the pipeline actually uses) rather than `human_g1k_v37`; chr22 coordinates match, but
  hs37d5's decoys legitimately attract some reads.
- 300,000 reads over 15 length bins (20-169 bp), not their 1.5 M.
- Scores mapping accuracy / precision / sensitivity only -- **not** reference bias at SNPs, nor the
  downstream D-statistics and PCA effects their paper also reports.
- Every read is endogenous (drawn from chr22). Real aDNA libraries are mostly exogenous, so these
  numbers describe accuracy on mappable reads and say nothing about runtime on real libraries.
