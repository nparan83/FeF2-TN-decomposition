#!/bin/bash
#SBATCH -J fef2_v3
#SBATCH -p hamsi
#SBATCH -A nparan
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=3-00:00:00
#SBATCH --array=1-5
#SBATCH --output=fef2_v3_%A_%a.out
# v3 production, one array task per U, geometry from S1 (u_relaxed_v3.txt). Configs in priority order
# (first NMAX, default 14). Orbital target from c111_ud m-scan; each 24-atom config is screened at
# 50/600 Ry to find the m that reaches the target, then run at full settings and checked. Resumable.
# v3b: a 24-atom config with no successful screening fails at once; full runs = screened m + at most 1 more.
export OMP_NUM_THREADS=1
if [ -z "$LOCALTEST" ]; then
module purge
module load comp/oneapi/2023
module load lib/hdf5/1.14.3-oneapi-2023.0
module load apps/espresso/7.2-oneapi-2023.0
fi
cd ${SLURM_SUBMIT_DIR:-.}
. ./fef2_run.sh
python3 -c "import numpy" || { echo "### NUMPY YOK"; exit 1; }
ULIST=(3.0 4.0 5.0 6.0 7.0)
U=${ULIST[$((${SLURM_ARRAY_TASK_ID:-2}-1))]}
[ -s u_relaxed_v3.txt ] || { echo "### u_relaxed_v3.txt YOK - once fef2_s1 bitmeli"; exit 1; }
D=v3_U$U
[ -d $D ] || python3 fef2_gen.py $U $(cat u_relaxed_v3.txt) 90 1080 $D > /dev/null
cd $D && python3 ../fef2_orbfix.py reorder ${NMAX:-14}
[ -f c111_ud_m1.in ] || python3 ../fef2_orbfix.py plan
for m in 1 2 3 4 5; do run_one c111_ud_m$m; done
python3 ../fef2_orbfix.py target > orb_target.log || { echo "### TARGET YOK"; exit 1; }
X=../fef2_orbfix.py
for f in $(awk '!/^#/{print $1}' configs.dat); do
  if python3 $X check $f; then echo "ok   $f"; continue; fi
  M=""
  case $f in c111_*|c112_*|c211_*) ;; *)
    for m in $(python3 $X order); do
      [ -f ${f}_s$m.in ] || python3 $X mkin $f $m s > /dev/null
      run_one ${f}_s$m
      if python3 $X check ${f}_s$m; then M=$m; echo "scr  $f m=$m"; break; fi
    done ;;
  esac
  done_=0
  case $f in
    c111_*|c112_*|c211_*) LIST=$(python3 $X order) ;;
    *) if [ -n "$M" ]; then LIST="$M $(python3 $X order | tr ' ' '\n' | grep -vx $M | head -1)"; else LIST=""; fi ;;
  esac
  for m in $LIST; do
    [ -f ${f}_m$m.in ] || python3 $X mkin $f $m > /dev/null
    run_one ${f}_m$m
    if python3 $X check ${f}_m$m; then python3 $X accept $f ${f}_m$m; echo "fix  $f m=$m"; done_=1; break; fi
  done
  [ $done_ = 1 ] || echo "FAIL $f"
done
echo "### V3 U=$U TAMAM $(date)"
