# fef2_run.sh - common functions (sourced by job scripts). ASCII only.
MPIRUN=${MPIRUN:-"mpirun"}
NK=${NK:-4}
run_one () {   # $1 = input name without .in ; skips finished runs, retries once with cg, beta 0.2, ndim 12
  local n=$1
  if [ -f $n.out ] && grep -q "JOB DONE" $n.out && grep -q "convergence has been achieved" $n.out; then
    echo "skip $n"; return 0; fi
  echo "run  $n  $(date +%H:%M:%S)"
  $MPIRUN pw.x -nk $NK -in $n.in > $n.out 2>&1
  if ! grep -q "convergence has been achieved" $n.out || ! grep -q "JOB DONE" $n.out; then
    echo "   retry $n with cg"
    mv $n.out ${n}_david_fail.out
    sed "s/^&ELECTRONS/\&ELECTRONS\n  diagonalization = 'cg'\n  mixing_ndim = 12/; s/mixing_beta = 0.3/mixing_beta = 0.2/; s/electron_maxstep = 300/electron_maxstep = 500/" $n.in > ${n}_cg.in
    $MPIRUN pw.x -nk $NK -in ${n}_cg.in > $n.out 2>&1
  fi
  rm -rf tmp_${n}_U*
}
