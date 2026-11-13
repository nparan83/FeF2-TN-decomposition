# fef2_mcB.py - Plan-B Monte Carlo driver: one job file -> one output file. Same engine (fef2_mc.py run(),
# heat bath / adaptive Metropolis, 16-colour update, NEQ = 4000 sweeps, NM = 8000 measurements every 2 sweeps)
# and same output columns as fef2_mcfine.py:  set model L rep T m m2 U4 dU4 e C chi acc
# Job file, one line per parameter set (J in meV, S = 2 convention of fef2_J.dat / 4):
#   name model slen D J1 J2 J3 J4 J7 T0 win nT Ls nrep rep0
#   T grid = T0 * (1 + linspace(-win, win, nT)); Ls comma separated; replicas rep0 .. rep0+nrep-1.
#   model: bc | heis (D = 0 -> full staggered vector order parameter) | ising1 (+-1, test)
# Seeds: crc32("B name model L rep T") - deterministic, independent of the stage-3 seeds.
# Usage: python3 fef2_mcB.py JOBFILE OUTFILE      (resumable: finished (set, model, L, rep, T) are skipped)
# ASCII only. numpy only.
import os, sys, zlib
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fef2_mc import run
NEQ = int(os.environ.get("NEQ", "4000")); NM = int(os.environ.get("NM", "8000")); NPROC = int(os.environ.get("NPROC", "56"))

def jobs(fn):
    out = []
    for l in open(fn):
        w = l.split()
        if not w or l.startswith("#"): continue
        name, model, slen, D = w[0], w[1], float(w[2]), float(w[3])
        J = dict(zip(("J1", "J2", "J3", "J4", "J7"), map(float, w[4:9])))
        T0, win, nT = float(w[9]), float(w[10]), int(w[11])
        Ls = [int(x) for x in w[12].split(",")]; nrep, rep0 = int(w[13]), int(w[14])
        for T in np.round(T0 * (1 + np.linspace(-win, win, nT)), 4):
            for L in Ls:
                for rep in range(rep0, rep0 + nrep): out.append((name, model, slen, D, J, L, rep, float(T)))
    return out

def job(a):
    name, model, slen, D, J, L, rep, T = a
    seed = zlib.crc32(("B %s %s %d %d %.4f" % (name, model, L, rep, T)).encode())
    r = run(L, model, J, D, T, NEQ, NM, seed, slen=slen)
    return "%-16s %-6s %3d %2d %9.4f %.6f %.6f %.6f %.6f %.6f %.6e %.6e %.3f" % (name, model, L, rep, T, r["m"], r["m2"], r["U4"], r["dU4"], r["e"], r["C"], r["chi"], r["acc"])

if __name__ == "__main__":
    jf, out = sys.argv[1], sys.argv[2]; tasks = jobs(jf); dn = set()
    if os.path.exists(out):
        for l in open(out):
            if not l.startswith("#"): w = l.split(); dn.add((w[0], w[1], int(w[2]), int(w[3]), "%.4f" % float(w[4])))
    else: open(out, "w").write("# set model L rep T m m2 U4 dU4 e C chi acc   (fef2_mcB.py %s, NEQ %d, NM %d)\n" % (jf, NEQ, NM))
    todo = [t for t in tasks if (t[0], t[1], t[5], t[6], "%.4f" % t[7]) not in dn]
    todo.sort(key=lambda t: -t[5]); print("%s: %d runs, %d to do" % (jf, len(tasks), len(todo)), flush=True)
    with Pool(NPROC) as p, open(out, "a") as g:
        for line in p.imap_unordered(job, todo): g.write(line + "\n"); g.flush()
    print("### MCB TAMAM %s" % jf, flush=True)
