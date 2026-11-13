#!/bin/bash
#SBATCH -J fef2_mcf
#SBATCH -p hamsi
#SBATCH -A nparan
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=12:00:00
#SBATCH --output=fef2_mcf_%j.out
# Precision Binder crossings (stage 3): neutron + linear-response-U sets, L = 16/24/32, 4 replicas,
# T0 +- 3 % in 0.5 % steps. Needs fef2_mc.py, fef2_mcrun.py, fef2_mcfine.py, fef2_J.dat, fef2_mcsum.txt.
# Resumable: resubmit to continue. numpy only.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd ${SLURM_SUBMIT_DIR:-.}
echo "Basladi: $(date)"
NPROC=${NPROC:-56} python3 fef2_mcfine.py > fef2_mcfine.log 2>&1
tail -12 fef2_mcfine.log
echo "Bitti: $(date)"
