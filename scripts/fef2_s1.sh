#!/bin/bash
#SBATCH -J fef2_s1
#SBATCH -p hamsi
#SBATCH -A nparan
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=08:00:00
#SBATCH --output=fef2_s1_%j.out
# v3 step 1: re-relax F position (AFM unit cell, U=4) with the minority electron started in each
# eigenvector m=1..5; the lowest-energy relaxed state defines the common geometry (u_relaxed_v3.txt).
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
[ -s u_relaxed.txt ] || { echo "### u_relaxed.txt YOK"; exit 1; }
python3 fef2_gen.py special 4.0 $(cat u_relaxed.txt) S1 > /dev/null
cd S1
for m in 1 2 3 4 5; do
  [ -f relax_afm_m$m.in ] || python3 ../fef2_orbfix.py mkin relax_afm $m > /dev/null
  run_one relax_afm_m$m
done
python3 ../fef2_orbfix.py relaxpick || exit 1
echo "### S1 TAMAM $(date)"
