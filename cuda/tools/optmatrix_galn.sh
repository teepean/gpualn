#!/bin/bash
# Byte-identity of G-ALN + prefilter against stock bwa aln across bwa option sets.
#
#   optmatrix_galn.sh OUTDIR "DATASETS"
#
# For every option set: stock `bwa aln -t 32` and `galn aln -dt 32` (prefilter build), whole-file
# md5 compared (same -t/-dt, so the headers match too). One line per pair in OUTDIR/matrix.tsv:
#   dataset  optset  bwa_ms  galn_ms  bwa_md5  galn_md5  SAME|DIFF  prefilter-line
set -u
OUT=$1 DATASETS=$2
ROOT=/home/teemu/sorsa/backtrackcuda
IDX=${IDX:-/home/dnastorage/aDNApipeline/hs37d5.fa}
BWA=${BWA:-/usr/bin/bwa}
GALN_PF=${GALN_PF:-/home/teemu/sorsa/galn-upstream/build-prefilter/galn}
mkdir -p "$OUT"
# name|options.  bwa's defaults are -n 0.04 -o 1 -e -1 -l 32 -k 2 -m 2000000 -M 3 -O 11 -E 4.
OPTSETS=(
  "default|"
  "adna|-l 1024 -n 0.01 -o 2"
  "adna-l16500|-l 16500 -n 0.01 -o 2"
  "seed20|-l 20 -k 1 -n 0.01 -o 2"
  "n1|-n 1 -o 1"
  "n2-nogap|-n 2 -o 0"
  "n3|-n 3 -o 2 -l 1024"
  "n4|-n 4 -o 2 -l 1024"
  "n0.2|-n 0.2 -o 1 -l 1024"
  "penalties|-l 1024 -n 0.01 -o 2 -M 1 -O 2 -E 1"
  "e3-nonGAPE|-l 1024 -n 0.01 -o 2 -e 3"
  "logGap|-l 1024 -n 0.01 -o 2 -L"
  "m1000|-l 1024 -n 0.01 -o 2 -m 1000"
  "i0-d100|-l 1024 -n 0.01 -o 2 -i 0 -d 100"
  "R1|-l 1024 -n 0.01 -o 2 -R 1"
)
for ds in $DATASETS; do
  fq=$ROOT/test_data/$ds.fq
  for os in "${OPTSETS[@]}"; do
    name=${os%%|*} opts=${os#*|}
    b=$OUT/$ds.$name.bwa.sai g=$OUT/$ds.$name.galn.sai
    t0=$(date +%s%N); $BWA aln -t 32 $opts "$IDX" "$fq" > "$b" 2> "$OUT/$ds.$name.bwa.log"; t1=$(date +%s%N)
    $GALN_PF aln -v -dt 32 $opts "$IDX" "$fq" > "$g" 2> "$OUT/$ds.$name.galn.log"; t2=$(date +%s%N)
    mb=$(md5sum < "$b" | cut -c1-8) mg=$(md5sum < "$g" | cut -c1-8)
    same=DIFF; [ "$mb" = "$mg" ] && same=SAME
    pf=$(grep -o "prefilter: [0-9]* of [0-9]* reads have no alignment, [0-9]* occur exactly" "$OUT/$ds.$name.galn.log" || echo "prefilter off")
    printf "%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n" "$ds" "$name" $(( (t1-t0)/1000000 )) $(( (t2-t1)/1000000 )) "$mb" "$mg" "$same" "$pf" | tee -a "$OUT/matrix.tsv"
  done
done
