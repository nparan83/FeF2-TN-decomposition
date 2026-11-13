#!/bin/bash
#SBATCH -J fef2_orb
#SBATCH -p hamsi
#SBATCH -A nparan
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=3-00:00:00
#SBATCH --array=1-5
#SBATCH --output=fef2_orb_%A_%a.out
# Orbital fix, one array task per U: scan starting eigenvector m=1..5 in c111_ud, lowest energy = target
# orbital; every config is checked and rerun (m in try order) until all Fe sit in the target. Resumable.
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
cd U$U || { echo "### U$U YOK"; exit 1; }
python3 ../fef2_orbfix.py plan
for m in 1 2 3 4 5; do run_one c111_ud_m$m; done
python3 ../fef2_orbfix.py target > orb_target.log || { echo "### TARGET YOK"; exit 1; }
if [ -z "$GEOM_OK" ] && ! python3 ../fef2_orbfix.py geomcheck; then
  echo "### S0 GEOMETRISI HEDEF YORUNGEDE DEGIL - durdu (orb_target.dat'i sohbete yapistir)"; exit 2; fi
for f in $(awk '!/^#/{print $1}' configs.dat); do
  if python3 ../fef2_orbfix.py check $f; then echo "ok   $f"; continue; fi
  fixed=0
  for m in $(python3 ../fef2_orbfix.py order); do
    [ -f ${f}_m$m.in ] || python3 ../fef2_orbfix.py mkin $f $m
    run_one ${f}_m$m
    if python3 ../fef2_orbfix.py check ${f}_m$m; then
      python3 ../fef2_orbfix.py accept $f ${f}_m$m; echo "fix  $f m=$m"; fixed=1; break; fi
  done
  [ $fixed = 1 ] || echo "FAIL $f"
done
echo "### ORB U=$U TAMAM $(date)"
