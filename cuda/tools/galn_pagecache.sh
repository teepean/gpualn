#!/bin/bash
# G-ALN 0.2.6+ drops the .bwt from the page cache after loading it (POSIX_FADV_DONTNEED at the end of
# load_bwt; reported as bilalbioinfo/galn#2). Measures the index load with and without that call, and
# what it does to the next program on the same index (bwa aln).
#
#   galn_pagecache.sh OUTDIR
#
# "Without" is the same binary with posix_fadvise made a no-op by an LD_PRELOAD shim (nofadvise.c).
# Fixed in the 0.2.7 release (79f951f): the cache is kept unless GALN_DROP_BWT_CACHE=1, checked last.
# Uses the GPU, unless EXTRA=--cpu-only (load_bwt is the same path; then no CUDA call at all).
set -u
O=$1; mkdir -p "$O"
IDX=${IDX:-/home/dnastorage/aDNApipeline/hs37d5.fa}
FQ=${FQ:-/home/teemu/sorsa/backtrackcuda/test_data/sub2k.fq}
GALN=${GALN:-/home/teemu/sorsa/galn-upstream/build-027/galn}
BWA=${BWA:-/usr/bin/bwa}
BOPTS="-l 1024 -n 0.01 -o 2"      # bwa's options; galn gets EXTRA on top
OPTS="$BOPTS ${EXTRA:-}"
SHIM=$O/nofadvise.so
gcc -shared -fPIC -O2 -o "$SHIM" "$(dirname "$0")/nofadvise.c" || exit 1
res()  { LC_ALL=C fincore -n -o RES "$IDX.bwt" | tr -d ' '; }
load() { grep -o "load [0-9.]* s" "$1"; }
tm()   { local s=$(date +%s%N); "$@" > /dev/null 2>&1; echo "$(( ($(date +%s%N) - s) / 1000000 )) ms"; }

echo "== galn as built"
cat "$IDX.bwt" > /dev/null; echo "cache warmed: .bwt resident $(res)"
for i in 1 2 3; do
  "$GALN" aln -v $OPTS "$IDX" "$FQ" > "$O/galn.$i.sai" 2> "$O/galn.$i.log"
  echo "run $i: $(load "$O/galn.$i.log"); .bwt resident afterwards $(res)"
done
echo "== galn with posix_fadvise made a no-op (LD_PRELOAD)"
cat "$IDX.bwt" > /dev/null; echo "cache warmed: .bwt resident $(res)"
for i in 1 2 3; do
  LD_PRELOAD=$SHIM "$GALN" aln -v $OPTS "$IDX" "$FQ" > "$O/shim.$i.sai" 2> "$O/shim.$i.log"
  echo "run $i: $(load "$O/shim.$i.log"); .bwt resident afterwards $(res)"
done
echo "== the next program: bwa aln on the same index"
cat "$IDX.bwt" > /dev/null
echo "bwa aln, warm cache:           $(tm "$BWA" aln -t 32 $BOPTS "$IDX" "$FQ")"
echo "bwa aln, warm cache again:     $(tm "$BWA" aln -t 32 $BOPTS "$IDX" "$FQ")"
"$GALN" aln $OPTS "$IDX" "$FQ" > /dev/null 2>&1
echo "bwa aln, right after galn aln: $(tm "$BWA" aln -t 32 $BOPTS "$IDX" "$FQ")"
echo "== galn with GALN_DROP_BWT_CACHE=1 (0.2.7+: the old behaviour on request)"
cat "$IDX.bwt" > /dev/null; echo "cache warmed: .bwt resident $(res)"
GALN_DROP_BWT_CACHE=1 "$GALN" aln -v $OPTS "$IDX" "$FQ" > "$O/drop.sai" 2> "$O/drop.log"
echo "run: $(load "$O/drop.log"); .bwt resident afterwards $(res)"
cat "$IDX.bwt" > /dev/null
echo "== .sai md5s (all six should match)"
md5sum "$O"/galn.*.sai "$O"/shim.*.sai | cut -c1-8 | sort | uniq -c
