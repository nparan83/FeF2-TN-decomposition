# fef2_mcfine.py - precision Binder-crossing MC for selected FeF2 parameter sets (stage 3).
# Sets (MCF_SETS, default neutron + the linear-response-U set) x models bc, heis; L = 16, 24, 32;
# 13 temperatures in T0*(1 + MCF_SHIFT +- MCF_WIN) (T0 from fef2_mcfine_T0.txt, else mean Binder crossing in fef2_mcsum.txt);
# NREP independent replicas (deterministic seeds). Output fef2_mcfine.dat; "sum" prints crossings per
# L pair (v6: local crossing of the replica-mean U4 curves, jackknife error; MCF_WIN, MCF_SHIFT set the window).
# Same Hamiltonian, D and units as fef2_mc.py / fef2_mcrun.py. ASCII only. numpy only.
import os, re, sys, zlib
import numpy as np
from multiprocessing import Pool
src = open("fef2_mcrun.py").read().split('\nif __name__ == "__main__":')[0]
exec(src)                                    # MC, run(), sets(), D, NPROC, QUICK
OUTF = "fef2_mcfine.dat"
T0F = {}                                     # optional window centres: fef2_mcfine_T0.txt lines 'set model T0'
if os.path.exists("fef2_mcfine_T0.txt"):
    for l in open("fef2_mcfine_T0.txt"):
        w = l.split()
        if len(w) == 3 and not l.startswith("#"): T0F[(w[0], w[1])] = float(w[2])
SETS = os.environ.get("MCF_SETS", "neutron,dft_U4.51_B14").split(",")
LS = (6, 8, 10) if QUICK else (16, 24, 32)
NREP = int(os.environ.get("NREP", "2" if QUICK else "4"))
WIN = float(os.environ.get("MCF_WIN", "0.03"))        # half-width of the T window (fraction of T0)
NEQ, NM = (100, 200) if QUICK else (4000, 8000)
def t0(name, model):
    for l in open("fef2_mcsum.txt"):
        w = l.split()
        if len(w) > 1 and w[0] == name and w[1] == model:
            c = [float(x) for x in re.findall(r"\d+/\d+:(\d+\.\d+)", l)]
            if c: return float(np.mean(c))
    return None
def job(a):
    name, model, L, rep, T, J = a
    seed = zlib.crc32(("fine %s %s %d %d %.4f" % (name, model, L, rep, T)).encode())
    r = run(L, model, J, D, T, NEQ, NM, seed)
    return "%-16s %-5s %3d %2d %9.4f %.6f %.6f %.6f %.6f %.6f %.6e %.6e %.3f" % (name, model, L, rep, T, r["m"], r["m2"], r["U4"], r["dU4"], r["e"], r["C"], r["chi"], r["acc"])
def cross(T, a, b):
    # Binder crossing: first temperature (scanning upward) where d = U4(L2) - U4(L1) changes from + (ordered side,
    # larger L closer to 2/3) to - (disordered side); linear interpolation between the two bracketing points.
    # (v6: replaces the straight-line fit over the whole window, which was biased low because U4 is flat on the
    #  ordered side.)
    d = b - a
    for i in range(len(T) - 1):
        if d[i] > 0 and d[i + 1] <= 0 and a[i] < 0.66:
            return T[i] + (T[i + 1] - T[i]) * d[i] / (d[i] - d[i + 1]), 0.5 * (a[i] + a[i + 1])
    return None, None
def summary():
    # all sets/models/replicas present in fef2_mcfine.dat; crossing of the replica-mean curves, jackknife error
    # (leave one replica out); EDGE = crossing in the outermost interval of the window (widen/shift and rerun).
    R = {}
    for l in open(OUTF):
        if l.startswith("#"): continue
        w = l.split(); R.setdefault((w[0], w[1], int(w[2])), {}).setdefault(round(float(w[4]), 4), {})[int(w[3])] = float(w[7])
    print("%-16s %-5s %-7s %s" % ("set", "model", "pair", "T_cross (K) +- jackknife (n replicas) | U4 at crossing, window"))
    for s, m in sorted({(k[0], k[1]) for k in R}, key=lambda x: (x[0] != "neutron", x[0], x[1])):
        Ls = sorted(k[2] for k in R if k[0] == s and k[1] == m)
        for L1, L2 in zip(Ls[:-1], Ls[1:]):
            A, B = R[(s, m, L1)], R[(s, m, L2)]; T = np.array(sorted(set(A) & set(B)))
            reps = sorted(set.intersection(*[set(A[t]) & set(B[t]) for t in T]))
            mean = lambda X, rr: np.array([np.mean([X[t][r] for r in rr]) for t in T])
            tc, u = cross(T, mean(A, reps), mean(B, reps)); jk = []
            for r0 in reps:
                rr = [r for r in reps if r != r0]
                if rr:
                    c, _ = cross(T, mean(A, rr), mean(B, rr))
                    if c is not None: jk.append(c)
            n = len(jk); err = np.sqrt((n - 1) / n * np.sum((np.array(jk) - np.mean(jk)) ** 2)) if n > 1 else float("nan")
            edge = " EDGE" if tc is not None and not (T[1] < tc < T[-2]) else ""
            print("%-16s %-5s %2d/%-4d %8.2f +- %.2f (%d) | U4*=%s  T %.2f-%.2f%s" % (s, m, L1, L2, tc if tc else float("nan"), err, len(reps),
                  "%.3f" % u if u else "-", T[0], T[-1], edge))
if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "sum": summary(); sys.exit(0)
    S = sets(); tasks = []
    for name in SETS:
        if name not in S: print("### set yok:", name, "mevcut:", " ".join(S)); sys.exit(1)
        for model in ("heis",):
            T0 = T0F.get((name, model)) or t0(name, model)
            if T0 is None: print("### T0 yok:", name, model); sys.exit(1)
            Ts = np.round(T0 * (1 + float(os.environ.get("MCF_SHIFT", "0")) + np.linspace(-WIN, WIN, 13)), 4)
            print(name, model, "T0 = %.2f K, grid %.2f - %.2f K" % (T0, Ts[0], Ts[-1]), flush=True)
            for L in LS:
                for rep in range(NREP):
                    for T in Ts: tasks.append((name, model, L, rep, float(T), S[name]))
    dn = set()
    if os.path.exists(OUTF):
        for l in open(OUTF):
            if not l.startswith("#"): w = l.split(); dn.add((w[0], w[1], int(w[2]), int(w[3]), "%.4f" % float(w[4])))
    else: open(OUTF, "w").write("# set model L rep T m m2 U4 dU4 e C chi acc\n")
    todo = [t for t in tasks if (t[0], t[1], t[2], t[3], "%.4f" % t[4]) not in dn]
    todo.sort(key=lambda t: -t[2]); print("runs: %d total, %d to do" % (len(tasks), len(todo)), flush=True)
    with Pool(NPROC) as p, open(OUTF, "a") as g:
        for line in p.imap_unordered(job, todo): g.write(line + "\n"); g.flush()
    print("### MCFINE TAMAM", flush=True); summary()
