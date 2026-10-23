#!/bin/bash
# Regression check for the Fortran signal code: build mktemplate and compute a
# few m_aa points with short VEGAS runs for several settings, then print the
# checksum of each output table. Two runs on different source versions must
# print identical output if the changes do not affect the results (comments,
# removal of unused declarations, ...).
#
# usage (inside the Docker image tmdp, repo mounted at /work):
#   scripts/check_regression.sh > before.txt
#   ... change the sources ...
#   scripts/check_regression.sh > after.txt;  diff before.txt after.txt
set -euo pipefail
REPO=${REPO:-/work}
make -s -C "$REPO/fortran" clean >/dev/null 2>&1 || true
make -s -C "$REPO/fortran" >/dev/null
EXE=$REPO/fortran/build/mktemplate.exe
W=$(mktemp -d); cd "$W"
run() { # $1 = tag, rest = namelist entries
  tag=$1; shift
  { echo "&TMPL"; for kv in "$@"; do echo "  $kv,"; done
    echo "  MAAMIN=330D0, MAAMAX=350D0, DMAA=5D0, NCALL=4000, ITMX=3,"
    echo "  OUTFILE='$tag.dat' /"; } > $tag.nml
  "$EXE" < $tag.nml > $tag.log 2>&1
  echo "$tag md5 $(md5sum < $tag.dat | cut -c1-32)  $(awk 'NR==3{print $2}' $tag.dat)"
}
run matched_nlo  "RS=13000D0" "MODE=3" "IORDG=1"
run matched_lo   "RS=13000D0" "MODE=3" "IORDG=0"
run oneloop      "RS=13000D0" "MODE=0"
run fcc_cuts     "RS=100000D0" "MODE=3" "PTRATIO=0.4D0"
rm -rf "$W"
