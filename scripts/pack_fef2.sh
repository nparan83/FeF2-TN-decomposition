#!/bin/bash
# pack_fef2.sh - FeF2 (JPCS) data and code package for GitHub. Copy + compress only: no calculation,
# nothing in the FeF2 directory is changed or deleted. Run on the login node (a few minutes).
#   bash pack_fef2.sh            -> /arf/scratch/nparan/repo_fef2/FeF2-TN-decomposition/  and  ...tar.gz
# Layout of the package:
#   scripts/  every *.py and *.sh of the FeF2 top directory (DFT generation, fit, QC, MC, EFT, Plan B)
#   data/     every *.dat, *.txt, *.log of the top directory (fef2_J.dat, fit logs, EFT table, MC summaries, hp_scf.dat ...)
#   mcB/ mcDC/  Monte Carlo directories (job files, raw Binder data, FSS results, scripts), files as they are
#   calc/<dir>.tar.gz  one archive per calculation directory: *.in *.out *.dat *.txt *.log (no wavefunctions,
#             charge densities, tmp_* or *.save); top-level *.in/*.out (job logs, hp.x) in calc/top_level.tar.gz
# Excluded directories (EXCL): S0 and U3.0..U7.0 = superseded v1/v2 runs in the wrong orbital state; alt = separate
#   study (not part of this paper); pseudo = PSlibrary files (named in the README instead).
# GitHub web upload limit is 25 MB per file: larger archives are split into 23 MB parts (cat parts > file).
set -u
SRC=${PACK_SRC:-/arf/scratch/nparan/FeF2}
TOP=${PACK_TOP:-/arf/scratch/nparan/repo_fef2}
NAME=FeF2-TN-decomposition
OUT=$TOP/$NAME
LIM=$((24*1024*1024))
EXCL=" S0 U3.0 U4.0 U5.0 U6.0 U7.0 alt pseudo "
cd "$SRC" || { echo "### $SRC YOK"; exit 1; }
mkdir -p "$TOP"; rm -rf "$OUT" "$TOP/$NAME.tar.gz"; mkdir -p "$OUT"/{scripts,data,calc}
: > "$TOP/skipped.txt"; : > "$TOP/inventory.txt"
pack () {   # $1 archive name, then the file list on stdin (NUL separated)
  tar czf "$OUT/calc/$1.tar.gz" --null -T -
  s=$(stat -c %s "$OUT/calc/$1.tar.gz")
  if [ "$s" -gt "$LIM" ]; then
    split -b 23m -d "$OUT/calc/$1.tar.gz" "$OUT/calc/$1.tar.gz.part"; rm "$OUT/calc/$1.tar.gz"
    echo "split: calc/$1.tar.gz ($((s/1048576)) MB) -> $(ls "$OUT"/calc/$1.tar.gz.part* | wc -l) parts" >> "$TOP/skipped.txt"
  fi
}
# 1) calculation directories
for d in */; do
  d=${d%/}; echo "$(du -sh "$d" 2>/dev/null | cut -f1)  $d" >> "$TOP/inventory.txt"
  case "$EXCL" in *" $d "*) echo "excluded dir: $d" >> "$TOP/skipped.txt"; continue ;; esac
  case "$d" in mcB|mcDC) continue ;; tmp_*|*.save|repo*) echo "excluded dir: $d" >> "$TOP/skipped.txt"; continue ;; esac
  n=$(find "$d" \( -name 'tmp_*' -o -name '*.save' \) -prune -o -type f \( -name '*.in' -o -name '*.out' -o -name '*.dat' \
        -o -name '*.txt' -o -name '*.log' \) -size -100M -print | wc -l)
  [ "$n" -gt 0 ] || { echo "empty (no in/out/dat/txt/log): $d" >> "$TOP/skipped.txt"; continue; }
  find "$d" \( -name 'tmp_*' -o -name '*.save' \) -prune -o -type f \( -name '*.in' -o -name '*.out' -o -name '*.dat' \
        -o -name '*.txt' -o -name '*.log' \) -size -100M -print0 | pack "$d"
  find "$d" \( -name 'tmp_*' -o -name '*.save' \) -prune -o -type f -size +100M -print | sed 's/^/over 100 MB, not packed: /' >> "$TOP/skipped.txt"
done
# 2) top level: scripts, derived data, inputs/outputs
find . -maxdepth 1 -type f \( -name '*.py' -o -name '*.sh' \) -exec cp {} "$OUT/scripts/" \;
find . -maxdepth 1 -type f \( -name '*.dat' -o -name '*.txt' -o -name '*.log' \) -size -24M -exec cp {} "$OUT/data/" \;
find . -maxdepth 1 -type f \( -name '*.dat' -o -name '*.txt' -o -name '*.log' \) -size +24M | sed 's/^/top-level over 24 MB, not copied: /' >> "$TOP/skipped.txt"
n=$(find . -maxdepth 1 -type f \( -name '*.in' -o -name '*.out' \) | wc -l)
[ "$n" -gt 0 ] && find . -maxdepth 1 -type f \( -name '*.in' -o -name '*.out' \) -print0 | pack top_level
# 3) Monte Carlo directories, files as they are (no SLURM restart or core files)
for d in mcB mcDC; do
  [ -d "$d" ] || { echo "missing: $d" >> "$TOP/skipped.txt"; continue; }
  mkdir -p "$OUT/$d"
  find "$d" -maxdepth 1 -type f -size -24M \( -name '*.py' -o -name '*.sh' -o -name '*.dat' -o -name '*.txt' -o -name '*.log' -o -name '*.out' \) \
       -exec cp {} "$OUT/$d/" \;
  find "$d" -maxdepth 1 -type f -size +24M | sed 's/^/over 24 MB, not copied: /' >> "$TOP/skipped.txt"
done
# 4) manifest and archive
(cd "$OUT" && find . -type f ! -name MANIFEST.txt -print0 | sort -z | xargs -0 md5sum > MANIFEST.txt)
cd "$TOP" && tar czf "$NAME.tar.gz" "$NAME"
echo "=== INVENTORY (top-level directories of $SRC)"; cat "$TOP/inventory.txt"
echo "=== SKIPPED / EXCLUDED / SPLIT"; cat "$TOP/skipped.txt"
echo "=== calc/ archives"; ls -l "$OUT/calc" | awk 'NR>1{printf "%9.2f MB  %s\n", $5/1048576, $9}'
echo "=== counts: scripts $(ls "$OUT/scripts" | wc -l), data $(ls "$OUT/data" | wc -l), mcB $(ls "$OUT/mcB" 2>/dev/null | wc -l), mcDC $(ls "$OUT/mcDC" 2>/dev/null | wc -l), files total $(find "$OUT" -type f | wc -l)"
echo "=== largest files"; find "$OUT" -type f -printf '%s %P\n' | sort -n | tail -5 | awk '{printf "%9.2f MB  %s\n", $1/1048576, $2}'
echo "=== package: $(du -sh "$OUT" | cut -f1) folder, $(du -sh "$TOP/$NAME.tar.gz" | cut -f1) tar.gz -> $TOP/$NAME.tar.gz"
echo "### PACK TAMAM $(date)"
