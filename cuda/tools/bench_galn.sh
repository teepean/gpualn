#!/bin/bash
# Head-to-head: stock bwa aln vs gpualn vs G-ALN (bilalbioinfo/galn), same box, same index, same reads.
#
#   bench_galn.sh OUTDIR REP "DATASETS" "CONFIGS"
#
# DATASETS: names of test_data/<name>.fq; CONFIGS: any of the keys in the case below.
# Every config uses all 32 hardware threads (bwa -t 32, gpualn -t 32; G-ALN always uses the whole
# box, -dt only sets the .sai header).  Appends one line per run to OUTDIR/results.tsv:
#   dataset config rep exit_code wall_ms md5(file) md5(records, i.e. without the 68-byte header)
# The header holds bwa's -t (n_threads, byte 57), so whole-file md5s are comparable only between
# runs that write the same thread count; the record md5 is comparable across all of them.
set -u
OUT=$1 REP=$2 DATASETS=$3 CONFIGS=$4
ROOT=/home/teemu/sorsa/backtrackcuda
IDX=${IDX:-/home/dnastorage/aDNApipeline/hs37d5.fa}
BWA=${BWA:-/usr/bin/bwa}
GPUALN=${GPUALN:-$ROOT/gpualn}
GPUALN_CPU=${GPUALN_CPU:-$ROOT/gpualn-cpu}
GALN=${GALN:-/home/teemu/sorsa/galn-upstream/build-sm86/galn}
GALN_PF=${GALN_PF:-/home/teemu/sorsa/galn-upstream/build-prefilter/galn}   # G-ALN + gpualn's scheme prefilter
OPTS="-l 1024 -n 0.01 -o 2"
T=${T:-32}
mkdir -p "$OUT"

for ds in $DATASETS; do
  fq=$ROOT/test_data/$ds.fq
  for cfg in $CONFIGS; do
    case $cfg in
      bwa)        cmd=("$BWA" aln -t $T $OPTS "$IDX" "$fq") ;;
      gpualn-cpu) cmd=("$GPUALN_CPU" aln -t $T $OPTS "$IDX" "$fq") ;;
      gpualn)     cmd=("$GPUALN" aln-gpu -t $T $OPTS "$IDX" "$fq") ;;
      galn)       cmd=("$GALN" aln -v -dt $T $OPTS "$IDX" "$fq") ;;
      galn-gpu)   cmd=("$GALN" aln -v --gpu-only -dt $T $OPTS "$IDX" "$fq") ;;
      galn-cpu)   cmd=("$GALN" aln -v --cpu-only -dt $T $OPTS "$IDX" "$fq") ;;
      galnpf)     cmd=("$GALN_PF" aln -v -dt $T $OPTS "$IDX" "$fq") ;;
      galnpf-gpu) cmd=("$GALN_PF" aln -v --gpu-only -dt $T $OPTS "$IDX" "$fq") ;;
      galnpf-verify) cmd=("$GALN_PF" aln -v --verify -dt $T $OPTS "$IDX" "$fq") ;;
      *) echo "unknown config $cfg" >&2; exit 1 ;;
    esac
    sai=$OUT/$ds.$cfg.r$REP.sai log=$OUT/$ds.$cfg.r$REP.log
    echo "${cmd[*]}" > "$log"
    t0=$(date +%s%N)
    "${cmd[@]}" > "$sai" 2>> "$log"
    rc=$?
    t1=$(date +%s%N)
    ms=$(( (t1 - t0) / 1000000 ))
    m_all=$(md5sum < "$sai" | cut -c1-8)
    m_rec=$(tail -c +69 "$sai" | md5sum | cut -c1-8)
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$ds" "$cfg" "$REP" "$rc" "$ms" "$m_all" "$m_rec" | tee -a "$OUT/results.tsv"
  done
done
