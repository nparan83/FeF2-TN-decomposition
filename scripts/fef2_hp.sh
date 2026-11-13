#!/bin/bash
#SBATCH -J fef2_hp
#SBATCH -p hamsi
#SBATCH -A nparan
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=1-00:00:00
#SBATCH --output=fef2_hp_%j.out
# Self-consistent linear-response U (hp.x, DFPT, ortho-atomic) for FeF2, AFM unit cell, v3 geometry.
# Loop: U_in -> c111_ud m-scan (target = lowest energy) -> scf (smearing) -> scf (fixed occ., 2-step as in
# HP/example02) -> hp.x q=2x2x2 -> U_out; U_in := U_out until |dU| < 0.05 eV (max 5). Then q=3x3x4 at the
# last U_in (q convergence). Table: hp_scf.dat. Resumable per iteration directory hp_U<U_in>.
export OMP_NUM_THREADS=1
if [ -z "$LOCALTEST" ]; then
module purge
module load comp/oneapi/2023
module load lib/hdf5/1.14.3-oneapi-2023.0
module load apps/espresso/7.2-oneapi-2023.0
fi
cd ${SLURM_SUBMIT_DIR:-.}
. ./fef2_run.sh
[ -s u_relaxed_v3.txt ] || { echo "### u_relaxed_v3.txt YOK"; exit 1; }
UIN=${UIN:-5.00}; Q1=${Q1:-"2 2 2"}; Q2=${Q2:-"3 3 4"}; HPCUT=${HPCUT:-90,1080}
hp_one () {   # $1 = U_in, $2 = q mesh, $3 = tag ; prints U_out
  local D=hp_U$1
  [ -d $D ] || python3 fef2_gen.py $1 $(cat u_relaxed_v3.txt) ${HPCUT%,*} ${HPCUT#*,} $D > /dev/null
  cd $D
  [ -f c111_ud_m1.in ] || python3 ../fef2_orbfix.py plan
  for m in 1 2 3 4 5; do run_one c111_ud_m$m > /dev/null; done
  python3 ../fef2_orbfix.py target > orb_target.log || { cd ..; echo "TARGETYOK"; return; }
  local m0=$(head -1 orb_target.dat | sed 's/.*(m=\([0-9]\)).*/\1/')
  if [ ! -s hpscf2.out ] || ! grep -q "JOB DONE" hpscf2.out; then
    sed "s/prefix = .*/prefix = 'hp'/; s#outdir = .*#outdir = './tmp_hp'#" c111_ud_m$m0.in > hpscf.in
    $MPIRUN pw.x -nk $NK -in hpscf.in > hpscf.out 2>&1
    python3 - << 'PY'
import re
t = open("hpscf.in").read()
t = t.replace("  occupations = 'smearing'\n  smearing = 'mv'\n  degauss = 0.005\n", "  occupations = 'fixed'\n  tot_magnetization = 0\n")
t = re.sub(r"  starting_ns_eigenvalue.*\n", "", t).replace("startingwfc = 'atomic'", "startingwfc = 'file'\n  startingpot = 'file'")
open("hpscf2.in", "w").write(t)
PY
    $MPIRUN pw.x -nk $NK -in hpscf2.in > hpscf2.out 2>&1
  fi
  python3 ../fef2_orbfix.py check hpscf || { cd ..; echo "ORBYOK"; return; }
  if [ ! -s hp_$3.out ] || ! grep -q "JOB DONE" hp_$3.out; then
    rm -rf tmp_hp/HP
    printf "&inputhp\n  prefix = 'hp'\n  outdir = './tmp_hp'\n  nq1 = %s, nq2 = %s, nq3 = %s\n  conv_thr_chi = 1.0d-6\n  iverbosity = 2\n/\n" $2 > hp_$3.in
    $MPIRUN hp.x -nk $NK -in hp_$3.in > hp_$3.out 2>&1
    cp hp.Hubbard_parameters.dat Hubbard_parameters_$3.dat 2>/dev/null
  fi
  local uo=$(awk '$3=="Fe1" && $7=="3d" {print $8; exit}' Hubbard_parameters_$3.dat 2>/dev/null)
  local e=$(grep "^!" hpscf.out | tail -1 | awk '{print $5}')
  echo "$1 $3 m=$m0 E=$e U_out=${uo:-NA}" >> ../hp_scf.dat
  cd ..; echo ${uo:-NA}
}
for it in $(seq 1 ${HPMAXIT:-5}); do
  ULAST=$UIN; UO=$(hp_one $UIN "$Q1" q1)
  echo "iter $it U_in=$UIN U_out=$UO $(date +%H:%M)"
  case $UO in NA|TARGETYOK|ORBYOK) echo "### HATA: $UO"; exit 1 ;; esac
  if python3 -c "import sys; sys.exit(0 if abs($UO-$UIN) < 0.05 else 1)"; then break; fi
  UIN=$(python3 -c "print('%.2f' % $UO)")
done
UO2=$(hp_one $ULAST "$Q2" q2)
echo "q-convergence at U_in=$ULAST: U_out(q1)=$UO  U_out(q2)=$UO2"
echo "### HP TAMAM $(date)"
