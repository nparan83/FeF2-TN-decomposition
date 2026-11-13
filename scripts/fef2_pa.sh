#!/bin/bash
#SBATCH -J fef2_pa
#SBATCH -p hamsi
#SBATCH -A nparan
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=3-00:00:00
#SBATCH --array=2-5
#SBATCH --output=fef2_pa_%A_%a.out
# Per-Fe orbital control, one array task per U (2=U4 ... 5=U7). Configs without a passing NAME.out are
# redone (default PA_MAXFE=4: 4-Fe cells + c221_uudddduu; PA_MAXFE=8: all): probe (1 SCF iteration)
# -> starting eigenvector m per Fe -> full run (greedy fallback). First a validation: a passing config recomputed with one species per Fe must give
# the same energy (< 0.01 meV/Fe) in the target orbital, else the task stops. Order: 4-Fe cells first,
# then 8-Fe (c221_uudddduu first: it makes set B solvable at U=6). Resumable.
export OMP_NUM_THREADS=1
if [ -z "$LOCALTEST" ]; then
module purge
module load comp/oneapi/2023
module load lib/hdf5/1.14.3-oneapi-2023.0
module load apps/espresso/7.2-oneapi-2023.0
fi
cd ${SLURM_SUBMIT_DIR:-.}
python3 -c "import numpy" || { echo "### NUMPY YOK"; exit 1; }
ULIST=(3.0 4.0 5.0 6.0 7.0)
U=${ULIST[$((${SLURM_ARRAY_TASK_ID:-3}-1))]}
cd v3_U$U || { echo "### v3_U$U YOK"; exit 1; }
X=../fef2_orbfix.py
python3 ../fef2_pa.py val ${PA_VAL:-c112_uudd} || { echo "### DOGRULAMA BASARISIZ - durdu"; exit 1; }
L=$(awk -v mx=${PA_MAXFE:-4} '!/^#/ && ($2<=mx || $1=="c221_uudddduu"){print $2, ($1=="c221_uudddduu"?0:1), $1}' configs.dat | sort -n -k1,1 -k2,2 | awk '{print $3}')
for f in ${PA_LIST:-$L}; do
  if python3 $X check $f; then continue; fi
  python3 ../fef2_pa.py search $f
done
echo "### PA U=$U TAMAM $(date)"
