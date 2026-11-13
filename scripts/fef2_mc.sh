#!/bin/bash
#SBATCH -J fef2_mc
#SBATCH -p hamsi
#SBATCH -A nparan
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=08:00:00
#SBATCH --output=fef2_mc_%j.out
# Classical MC (Blume-Capel S=2 and Heisenberg |S|=2 + D) for the neutron set and every U in fef2_J.dat.
# Resumable: resubmit to continue. Needs fef2_mc.py, fef2_mcrun.py, (fef2_J.dat). numpy only.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd ${SLURM_SUBMIT_DIR:-.}
python3 -c "import numpy as n; v=tuple(int(x) for x in n.__version__.split('.')[:2]); assert v>=(1,17), v" || { echo "### NUMPY YOK/ESKI"; exit 1; }
echo "Basladi: $(date)"
NPROC=${NPROC:-56} python3 fef2_mcrun.py > fef2_mcrun.log 2>&1
tail -3 fef2_mcrun.log
python3 fef2_mcsum.py > fef2_mcsum.txt 2>&1; cat fef2_mcsum.txt
echo "Bitti: $(date)"
