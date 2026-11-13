#!/bin/bash
#SBATCH -J fef2_s0
#SBATCH -p hamsi
#SBATCH -A nparan
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=08:00:00
#SBATCH --output=fef2_s0_%j.out
# Step 0: relax F position (AFM unit cell, U=4) -> convergence tests + commensurate-k check
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
U0=4.0
python3 fef2_gen.py special $U0 0.30 S0 > /dev/null
cd S0
run_one relax_afm
grep -q "End of BFGS Geometry Optimization" relax_afm.out || { echo "### RELAX BITMEDI"; exit 1; }
u=$(sed -n '/Begin final coordinates/,/End final coordinates/p' relax_afm.out | awk '$1=="F"{print $2; exit}')
echo "$u" > ../u_relaxed.txt; echo "u(F) relaxed = $u"
cd .. && python3 fef2_gen.py special $U0 $u S0 > /dev/null && cd S0
run_one cv_e90_afm
python3 -c "
import re; t=open('cv_e90_afm.out').read(); i=open('cv_e90_afm.in').read()
f=lambda k: re.search(k+r' grid:.*FFT dimensions: \(\s*(\d+),\s*(\d+),\s*(\d+)\)', t).groups()
e=lambda k: re.search(k+r' = (\S+)', i).group(1)
print(e('ecutwfc'), e('ecutrho'), *f('Dense '), *f('Smooth'))" > ../fft_base.txt
echo "FFT base: $(cat ../fft_base.txt)"
cd .. && python3 fef2_gen.py special $U0 $u S0 > /dev/null && cd S0
for n in cv_e70_afm cv_e70_fm cv_e90_afm cv_e90_fm cv_e110_afm cv_e110_fm cv_k8812_afm cv_k8812_fm cv_orb4_afm chk212_afm chk212_fm; do
  run_one $n
done
echo "### S0 TAMAM $(date)"
