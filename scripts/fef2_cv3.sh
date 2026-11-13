#!/bin/bash
#SBATCH -J fef2_cv3
#SBATCH -p hamsi
#SBATCH -A nparan
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=08:00:00
#SBATCH --output=fef2_cv3_%j.out
# Convergence tests in the TARGET orbital at the v3 geometry (U=4): cutoff 70/90/110 Ry, k 8x8x12,
# AFM and FM; commensurate-k/FFT check 2x1x2 vs 1x1x1. Each input tries m until the orbital matches.
export OMP_NUM_THREADS=1
if [ -z "$LOCALTEST" ]; then
module purge
module load comp/oneapi/2023
module load lib/hdf5/1.14.3-oneapi-2023.0
module load apps/espresso/7.2-oneapi-2023.0
fi
cd ${SLURM_SUBMIT_DIR:-.}
. ./fef2_run.sh
[ -s u_relaxed_v3.txt ] && [ -s v3_U4.0/orb_target.dat ] || { echo "### v3 ciktilari YOK"; exit 1; }
[ -d CV3 ] || python3 fef2_gen.py special 4.0 $(cat u_relaxed_v3.txt) CV3 > /dev/null
cp v3_U4.0/orb_target.dat v3_U4.0/orb_morder.txt CV3/
cd CV3 && X=../fef2_orbfix.py
for n in cv_e70_afm cv_e70_fm cv_e90_afm cv_e90_fm cv_e110_afm cv_e110_fm cv_k8812_afm cv_k8812_fm chk212_afm chk212_fm; do
  if python3 $X check $n; then echo "ok   $n"; continue; fi
  case $n in chk212_*) L=$(python3 $X order | cut -d' ' -f1-2) ;; *) L=$(python3 $X order) ;; esac
  ok=0
  for m in $L; do
    [ -f ${n}_m$m.in ] || python3 $X mkin $n $m > /dev/null
    run_one ${n}_m$m
    if python3 $X check ${n}_m$m; then python3 $X accept $n ${n}_m$m; echo "fix  $n m=$m"; ok=1; break; fi
  done
  [ $ok = 1 ] || echo "FAIL $n"
done
cd .. && python3 fef2_conv.py CV3
echo "### CV3 TAMAM $(date)"
