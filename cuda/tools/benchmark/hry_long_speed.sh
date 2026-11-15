#!/bin/bash
# Long-arm speed: the pipeline's bwa-mem3 long arm vs strobealign-HRY (paper settings) on reads >= 64 bp
# of real libraries. The short arm (gpualn) is the same in both combos, so this is the whole difference.
# Strictly one run at a time, -t 32, indexes warmed into the page cache first, 3 runs per tool.
#   hry_long_speed.sh OUTDIR LIB1.fq [LIB2.fq ...]
set -u
O=$1; shift; mkdir -p "$O"; cd "$O" || exit 1
BREF=/home/dnastorage/aDNApipeline/hs37d5.fa          # bwa / bwa-mem3 index
HREF=/home/teemu/sorsa/adna_simbench/ref/hs37d5.fa    # HRY .sti next to this symlink
SA=/home/teemu/sorsa/hry/strobealign-hry/target/release/strobealign
M3=/home/teemu/sorsa/bwa-mem3-v0141/bwa-mem3
HOPT="--adna --ry-len 8 -M 500 --ssw -k 17 -s 15 -r 50"
T=32
export LC_ALL=C
secs() { awk '/Elapsed \(wall clock\)/{n=split($NF,a,":"); s=0; for(i=1;i<=n;i++) s=s*60+a[i]; printf "%.2f", s}' "$1"; }
rss()  { awk '/Maximum resident/{printf "%.1f GB", $NF/1048576}' "$1"; }
stats() {  # primary records: mapped, MAPQ>=25, MAPQ>=30, and (HRY) AS>=60
  awk -F'\t' '!/^@/ && !and($2,0x900){n++; if(!and($2,4)){m++; if($5>=25)q25++; if($5>=30)q30++;
      as=-1; for(i=12;i<=NF;i++) if(substr($i,1,5)=="AS:i:"){as=substr($i,6)+0; break}
      if(as>=60){a++; if($5>=25)aq++}}}
    END{printf "reads %d | mapped %.2f%% MAPQ>=25 %.2f%% MAPQ>=30 %.2f%% | AS>=60 %.2f%% +MAPQ>=25 %.2f%%\n",
        n, 100*m/n, 100*q25/n, 100*q30/n, 100*a/n, 100*aq/n}' "$1"; }

for FQ in "$@"; do
  lib=$(basename "$FQ" .fq)
  echo "== $lib ($(( $(wc -l < "$FQ") / 4 )) reads)"
  cat $BREF.bwt.2bit.64 $BREF.0123 $BREF.pac $BREF.ann $BREF.amb "$FQ" > /dev/null
  for i in 1 2 3; do
    /usr/bin/time -v $M3 mem -t $T -L 15 --supp-rep-hard-cap 5 $BREF "$FQ" > $lib.mem3.sam 2> $lib.mem3.$i.log
    echo "bwa-mem3 -L 15 run $i:   $(secs $lib.mem3.$i.log) s, peak $(rss $lib.mem3.$i.log)"
  done
  cat "$HREF.sti" "$FQ" > /dev/null
  for i in 1 2 3; do
    /usr/bin/time -v $SA -v -t $T $HOPT --use-index "$HREF" "$FQ" > $lib.hry.sam 2> $lib.hry.$i.log
    echo "HRY --use-index run $i:  $(secs $lib.hry.$i.log) s, peak $(rss $lib.hry.$i.log)"
  done
  echo "mem3: $(stats $lib.mem3.sam)"
  echo "HRY:  $(stats $lib.hry.sam)"
done
echo done
