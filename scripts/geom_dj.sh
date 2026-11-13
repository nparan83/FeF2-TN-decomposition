#!/bin/bash
#SBATCH -J fef2_gdj
#SBATCH -p hamsi
#SBATCH -A nparan
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=02:00:00
#SBATCH --output=fef2_gdj_%j.out
# Geometry check of the fixed-geometry hp.x loop: FM-AFM (= 8(J2+J5)S^2 per Fe) at U=4.51 with u relaxed at
# U=4.51 (geomchk/u_relaxed_v3.txt) vs the production geometry (v3_U4.51, u relaxed at U=4). Probe orbital control.
export OMP_NUM_THREADS=1
module purge
module load comp/oneapi/2023
module load lib/hdf5/1.14.3-oneapi-2023.0
module load apps/espresso/7.2-oneapi-2023.0
cd ${SLURM_SUBMIT_DIR:-.}/geomchk
cp ../fef2_pa.py .
[ -d G451 ] || python3 fef2_gen.py 4.51 $(cat u_relaxed_v3.txt) 90 1080 G451 > /dev/null
cd G451 && cp ../../v3_U4.51/orb_target.dat ../../v3_U4.51/orb_morder.txt .
for n in c111_ud c111_uu; do python3 ../fef2_pa.py search $n; done
python3 - << 'PY'
import re
RY = 13605.693122994
E = lambda f: float(re.findall(r"^!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", open(f).read(), re.M)[-1])
a = (E("c111_uu.out") - E("c111_ud.out")) / 2 * RY
b = (E("../../v3_U4.51/c111_uu.out") - E("../../v3_U4.51/c111_ud.out")) / 2 * RY
print("FM-AFM meV/Fe: u(U=4.51) %.4f   u(U=4, production) %.4f   relative change %+.3f %%" % (a, b, 100 * (a / b - 1)))
PY
echo "### GDJ TAMAM $(date)"
