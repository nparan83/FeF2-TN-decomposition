#!/bin/bash
#SBATCH -J fef2_uhp
#SBATCH -p hamsi
#SBATCH -A nparan
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=2-00:00:00
#SBATCH --output=fef2_uhp_%j.out
# Production at the linear-response U (set UHP explicitly, e.g. sbatch --export=ALL,UHP=4.30 fef2_uhp.sh).
# v3 geometry, first NMAX (14) priority configs. Target orbital from the c111_ud m-scan; every other config:
# probe (1 SCF iteration) -> starting eigenvector m per Fe -> full run (fef2_pa.py search). Resumable.
export OMP_NUM_THREADS=1
if [ -z "$LOCALTEST" ]; then
module purge
module load comp/oneapi/2023
module load lib/hdf5/1.14.3-oneapi-2023.0
module load apps/espresso/7.2-oneapi-2023.0
fi
cd ${SLURM_SUBMIT_DIR:-.}
. ./fef2_run.sh
[ -n "$UHP" ] || { echo "### UHP verilmedi"; exit 1; }
[ -s u_relaxed_v3.txt ] || { echo "### u_relaxed_v3.txt YOK"; exit 1; }
D=v3_U$UHP
[ -d $D ] || python3 fef2_gen.py $UHP $(cat u_relaxed_v3.txt) 90 1080 $D > /dev/null
cd $D && python3 ../fef2_orbfix.py reorder ${NMAX:-14}
[ -f c111_ud_m1.in ] || python3 ../fef2_orbfix.py plan
for m in 1 2 3 4 5; do run_one c111_ud_m$m; done
python3 ../fef2_orbfix.py target > orb_target.log || { echo "### TARGET YOK"; exit 1; }
L=$(awk '!/^#/{print $2, $1}' configs.dat | sort -n -k1,1 | awk '{print $2}')
for f in $L; do
  if python3 ../fef2_orbfix.py check $f; then echo "ok   $f"; continue; fi
  python3 ../fef2_pa.py search $f
done
echo "### UHP U=$UHP TAMAM $(date)"
