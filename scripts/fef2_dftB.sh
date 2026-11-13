#!/bin/bash
#SBATCH -J fef2_dftB
#SBATCH -p hamsi
#SBATCH -A nparan
#SBATCH -N 1
#SBATCH -n 56
#SBATCH --time=3-00:00:00
#SBATCH --output=fef2_dftB_%x_%j.out
# Plan-B DFT control calculations, U = 4.51 eV unless stated, v3 geometry (u = u_relaxed_v3.txt), production
# settings (90/1080 Ry, conv 1e-10, MV 0.005 Ry), per-Fe orbital control (fef2_pa.py search). Resumable.
# Submit from the FeF2 root:  sbatch -J <name> --export=ALL,MODE=<mode> fef2_dftB.sh
#   MODE=lat     FM-AFM of the unit cell with cited room-temperature lattice constants:
#                v7_lat_SR  a = 4.6966, c = 3.3091 A (Stout & Reed 1954) ; v7_lat_BK  a = 4.6945, c = 3.3097 A (Baur & Khan 1971)
#   MODE=csub    J1..J7 from the 8-configuration subset at c x 0.99 (v7_c099) and c x 1.01 (v7_c101), a fixed;
#                commensurate FFT grids from a unit-cell probe at each c; reference: same subset in v3_U4.51
#   MODE=s113    1x1x3 cell (k 6x6x3 = 6x6x9 per unit cell), 8 configurations: J2 and J5 separately (v7_113)
#   MODE=atomsub non-orthogonalized atomic projectors: 8-config subset + FM-AFM at the self-consistent U of
#                atom/hp_scf.dat (run atom/fef2_hp.sh first; directory atom/v7a_U<U>)
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
[ -s u_relaxed_v3.txt ] || { echo "### u_relaxed_v3.txt YOK"; exit 1; }
UF=$(cat u_relaxed_v3.txt); ECW=${ECW:-90}; ECR=${ECR:-1080}
SUB="c111_uu c111_ud c112_uudd c112_uuud c211_uudd c211_uuud c212_uudddduu c221_uudddduu"
target_dir () {      # inside a config dir: c111_ud m-scan -> orb_target.dat (skipped if present)
  [ -s orb_target.dat ] && return 0
  [ -f c111_ud_m1.in ] || python3 ../fef2_orbfix.py plan
  for m in 1 2 3 4 5; do run_one c111_ud_m$m; done
  python3 ../fef2_orbfix.py target > orb_target.log || { echo "### TARGET YOK $(pwd)"; return 1; }
}
do_list () {         # inside a config dir: every config -> check, else per-Fe probe search
  for f in "$@"; do
    if python3 ../fef2_orbfix.py check $f; then echo "ok   $f"; continue; fi
    python3 ../fef2_pa.py search $f
  done
}
fftprobe () {        # $1 = probe dir (sibling of the config dir), env FEF2_A/FEF2_C set -> $2 = fft_base file
  [ -s $2 ] && return 0
  FEF2_FFTBASE=none python3 fef2_genB.py 4.51 $UF $ECW $ECR $1 c111_ud > /dev/null
  ( cd $1 && sed 's/electron_maxstep = 300/electron_maxstep = 1/' c111_ud.in > p.in && $MPIRUN pw.x -nk $NK -in p.in > p.out 2>&1; rm -rf tmp_* )
  python3 - $1/p.out $2 $ECW $ECR << 'PY'
import re, sys
t = open(sys.argv[1]).read()
d = re.search(r"Dense\s+grid:.*?FFT dimensions:\s*\(\s*(\d+),\s*(\d+),\s*(\d+)\)", t)
s = re.search(r"Smooth\s+grid:.*?FFT dimensions:\s*\(\s*(\d+),\s*(\d+),\s*(\d+)\)", t) or d
if not d: sys.exit("### FFT boyutu bulunamadi: " + sys.argv[1])
open(sys.argv[2], "w").write("%s %s %s %s %s %s %s %s\n" % ((sys.argv[3], sys.argv[4]) + d.groups() + s.groups()))
print("fft_base", sys.argv[2], open(sys.argv[2]).read().strip())
PY
}
case $MODE in
lat)
  for X in "SR 4.6966 3.3091" "BK 4.6945 3.3097"; do
    set -- $X; D=v7_lat_$1
    [ -d $D ] || FEF2_A=$2 FEF2_C=$3 FEF2_FFTBASE=none python3 fef2_genB.py 4.51 $UF $ECW $ECR $D c111_ud c111_uu
    ( cd $D && target_dir && do_list c111_uu && python3 ../fef2_jsub.py FA | tee jsub_FA.txt )
  done
  echo "reference (production lattice 4.6974 / 3.3082):"; ( cd v3_U4.51 && python3 ../fef2_jsub.py FA ) ;;
csub)
  for X in "c099 0.99" "c101 1.01"; do
    set -- $X; D=v7_$1; CC=$(python3 -c "print('%.6f' % (3.3082 * $2))")
    export FEF2_C=$CC
    fftprobe ${D}_fft ${D}_fft_base.txt || exit 1
    [ -d $D ] || FEF2_FFTBASE=${D}_fft_base.txt python3 fef2_genB.py 4.51 $UF $ECW $ECR $D $SUB
    unset FEF2_C
    ( cd $D && target_dir && do_list $SUB && python3 ../fef2_jsub.py B $SUB | tee jsub_B.txt && python3 ../fef2_jsub.py FA | tee jsub_FA.txt )
  done
  echo "reference (production c, same 8 configs):"; ( cd v3_U4.51 && python3 ../fef2_jsub.py B $SUB | tee jsub_B_sub8.txt ) ;;
s113)
  D=v7_113
  [ -d $D ] || FEF2_CELLS=113 FEF2_KBASE=6,6,9 python3 fef2_genB.py 4.51 $UF $ECW $ECR $D
  [ -s $D/orb_target.dat ] || cp v3_U4.51/orb_target.dat v3_U4.51/orb_morder.txt $D/
  L=$(awk '!/^#/{print $1}' $D/configs.dat)
  ( cd $D && do_list $L && python3 ../fef2_jsub.py 113 | tee jsub_113.txt && python3 ../fef2_jsub.py FA c113_uuuuuu c113_ududud | tee jsub_FA.txt ) ;;
atomsub)
  cd atom || { echo "### atom/ YOK"; exit 1; }
  UA=$(awk '$2=="q2"{u=$5} END{sub("U_out=","",u); print u}' hp_scf.dat)
  [ -n "$UA" ] || { echo "### atom/hp_scf.dat'ta q2 satiri yok (fef2_hp.sh bitmedi)"; exit 1; }
  UA=$(python3 -c "print('%.2f' % $UA)"); D=v7a_U$UA; echo "atomic projectors: U = $UA eV"
  [ -d $D ] || python3 fef2_genB.py $UA $UF $ECW $ECR $D $SUB
  ( cd $D && target_dir && do_list $SUB && python3 ../fef2_jsub.py B $SUB | tee jsub_B.txt && python3 ../fef2_jsub.py FA | tee jsub_FA.txt ) ;;
*) echo "### MODE = lat | csub | s113 | atomsub"; exit 1 ;;
esac
echo "### DFTB $MODE TAMAM $(date)"
