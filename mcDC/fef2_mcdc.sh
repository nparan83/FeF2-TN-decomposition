#!/bin/bash
#SBATCH -J fef2_mcdc
#SBATCH -p hamsi
#SBATCH -A nparan
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=1-00:00:00
#SBATCH --output=fef2_mcdc_%j.out
# Heisenberg-only Binder crossings with D' = 0.875*D (El Mendili & Zhitomirsky D_C), |S| = 2; T x 1.5 afterwards.
# Resumable: resubmit to continue.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd ${SLURM_SUBMIT_DIR:-.}
echo "Basladi: $(date)"
MCF_SETS=neutron,dft_U4.51_B14 NREP=12 MCF_WIN=0.03 NPROC=56 python3 fef2_mcfine.py > mcdc_a.log 2>&1
MCF_SETS=dft_U3.0_B14,dft_U4.0_B14,dft_U5.0_B14,dft_U6.0_B14,dft_U7.0_B14 NREP=4 MCF_WIN=0.03 NPROC=56 python3 fef2_mcfine.py > mcdc_b.log 2>&1
python3 fef2_mcfine.py sum > fef2_mcdc_sum.txt; cat fef2_mcdc_sum.txt
echo "Bitti: $(date)"
