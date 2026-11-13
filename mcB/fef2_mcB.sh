#!/bin/bash
#SBATCH -J fef2_mcB
#SBATCH -p hamsi
#SBATCH -A nparan
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=1-00:00:00
#SBATCH --output=fef2_mcB_%x_%j.out
# Plan-B Monte Carlo (run inside mcB/). Submit:  sbatch -J <name> --export=ALL,JOB=<job> fef2_mcB.sh
#   JOB = tau | focus | scan5 | d0 | val | syn   -> mcB_<JOB>.txt -> mcB_<JOB>.dat  (tau: fef2_tau.py -> tau.dat)
# Resumable: resubmit the same command to continue. numpy only.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd ${SLURM_SUBMIT_DIR:-.}
echo "Basladi: $(date)  JOB=$JOB"
[ -s mcB_$JOB.txt ] || { echo "### mcB_$JOB.txt YOK (once python3 fef2_mkjobsB.py)"; exit 1; }
if [ "$JOB" = tau ]; then NPROC=${NPROC:-56} python3 fef2_tau.py mcB_tau.txt tau.dat > tau.log 2>&1; tail -20 tau.log
else NPROC=${NPROC:-56} python3 fef2_mcB.py mcB_$JOB.txt mcB_$JOB.dat > mcB_$JOB.log 2>&1; tail -3 mcB_$JOB.log; fi
echo "Bitti: $(date)"
