#!/bin/bash
# Map the simulated set (sim15M.fq in WORKDIR) with bwa aln, the pipeline's bwa-mem3 long arm,
# strobealign-HRY at the paper's settings and three other strobealign configurations. One SAM and
# one /usr/bin/time -v log per strategy; runs strictly one after another. Accuracy, not speed: the
# strobealign rows build their index on every run. Real-library speed: hry_speed.sh.
#   hry_run.sh WORKDIR
set -u
W=$1; cd "$W" || exit 1
REF=${REF:-/home/dnastorage/aDNApipeline/hs37d5.fa}
FQ=$W/sim15M.fq
SA=${SA:-/home/teemu/sorsa/hry/strobealign-hry/target/release/strobealign}   # tag hry = 4be2dcc
M3=${M3:-/home/teemu/sorsa/bwa-mem3-v0141/bwa-mem3}
T=16
info() { printf "%-8s wall %s  peak %s kB\n" "$1" "$(awk '/Elapsed/{print $NF}' "$2")" "$(awk '/Maximum resident/{print $NF}' "$2")"; }
run() {   # name, command...
  local name=$1; shift
  [ -s "big.$name.sam" ] && { echo "$name: exists, skipped"; return; }
  /usr/bin/time -v "$@" > "big.$name.sam.tmp" 2> "big.$name.log" && mv "big.$name.sam.tmp" "big.$name.sam"
  info "$name" "big.$name.log"
}
cat $REF.bwt $REF.sa $REF.pac > /dev/null
if [ ! -s big.aln.sam ]; then
  /usr/bin/time -v bwa aln -l 1024 -n 0.01 -o 2 -t $T $REF $FQ > big.aln.sai 2> big.aln.log
  bwa samse $REF big.aln.sai $FQ > big.aln.sam 2> big.samse.log
  info aln big.aln.log
fi
run mem    $M3 mem -t $T -L 15 --supp-rep-hard-cap 5 $REF $FQ                   # pipeline long arm
run hry    $SA -v -t $T --adna --ry-len 8 -M 500 --ssw -k 17 -s 15 $REF $FQ     # paper's HRY
run hry50  $SA -v -t $T --adna --ry-len 8 -M 500 --ssw -k 17 -s 15 -r 50 $REF $FQ  # = hry (checked)
run hrydef $SA -v -t $T --adna $REF $FQ                                          # release default
run dk     $SA -v -t $T -l 1 -u 1 -M 500 -k 17 -s 15 $REF $FQ                    # paper's DK
run strobe $SA -v -t $T $REF $FQ                                                 # plain strobealign
