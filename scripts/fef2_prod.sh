#!/bin/bash
#SBATCH -J fef2_prod
#SBATCH -p hamsi
#SBATCH -A nparan
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=16:00:00
#SBATCH --array=1-5
#SBATCH --output=fef2_prod_%A_%a.out
# Production: 27 collinear configurations at U = 3, 4, 5, 6, 7 eV (one array task per U), common geometry from S0
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
[ -s u_relaxed.txt ] || { echo "### u_relaxed.txt YOK - once fef2_s0 bitmeli"; exit 1; }
u=$(cat u_relaxed.txt)
[ -d U$U ] || python3 fef2_gen.py $U $u ${ECW:-90} ${ECR:-1080} U$U
cd U$U
for f in $(awk '!/^#/{print $1}' configs.dat); do run_one $f; done
echo "### U=$U TAMAM $(date)"
