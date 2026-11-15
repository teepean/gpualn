#!/bin/bash
# Real-library speed: strobealign-HRY (paper settings) vs CPU bwa aln, -t 32, strictly one run at a time.
#   hry_speed.sh OUTDIR
set -u
O=$1; mkdir -p "$O"; cd "$O" || exit 1
BREF=/home/dnastorage/aDNApipeline/hs37d5.fa          # bwa index
HREF=/home/teemu/sorsa/adna_simbench/ref/hs37d5.fa    # symlink; HRY writes its .sti next to it
SA=/home/teemu/sorsa/hry/strobealign-hry/target/release/strobealign
TD=/home/teemu/sorsa/backtrackcuda/test_data
HOPT="--adna --ry-len 8 -M 500 --ssw -k 17 -s 15 -r 50"
T=32
secs() { awk '/Elapsed \(wall clock\)/{n=split($NF,a,":"); s=0; for(i=1;i<=n;i++) s=s*60+a[i]; printf "%.2f", s}' "$1"; }
rss()  { awk '/Maximum resident/{printf "%.1f GB", $NF/1048576}' "$1"; }

echo "== HRY index build (-r 50 -> short-read profile)"
if [ ! -s "$HREF.sti" ]; then
  /usr/bin/time -v $SA --create-index -t $T $HOPT "$HREF" > /dev/null 2> index.log
  echo "index: $(secs index.log) s, peak $(rss index.log), $(du -h "$HREF.sti" | cut -f1) on disk"
else echo "index exists"; fi

for lib in sub100k p17_shuf1m p1_shuf1m; do
  FQ=$TD/$lib.fq
  echo "== $lib ($(( $(wc -l < $FQ) / 4 )) reads)"
  cat $BREF.bwt $BREF.sa $BREF.pac $FQ > /dev/null
  /usr/bin/time -v bwa aln -t $T -l 1024 -n 0.01 -o 2 $BREF $FQ > $lib.aln.sai 2> $lib.aln.log
  echo "bwa aln -t $T:            $(secs $lib.aln.log) s, peak $(rss $lib.aln.log)"
  cat "$HREF.sti" $FQ > /dev/null
  for i in 1 2 3; do
    /usr/bin/time -v $SA -v -t $T $HOPT --use-index "$HREF" $FQ > $lib.hry.sam 2> $lib.hry.$i.log
    echo "HRY --use-index run $i:   $(secs $lib.hry.$i.log) s, peak $(rss $lib.hry.$i.log)"
  done
  /usr/bin/time -v $SA -v -t $T --adna --ry-len 8 -M 500 --ssw -k 17 -s 15 "$HREF" $FQ > $lib.hryfly.sam 2> $lib.hryfly.log
  echo "HRY index on the fly, -r estimated: $(secs $lib.hryfly.log) s, peak $(rss $lib.hryfly.log)"
done
echo "== output identity: --use-index -r 50 vs on-the-fly with estimated -r (records only)"
for lib in sub100k p17_shuf1m p1_shuf1m; do
  a=$(grep -v '^@' $lib.hry.sam | md5sum | cut -c1-8); b=$(grep -v '^@' $lib.hryfly.sam | md5sum | cut -c1-8)
  echo "$lib: $a $b $([ $a = $b ] && echo same || echo DIFFERENT)"
done
echo done
