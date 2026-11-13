# fef2_tau.py - integrated autocorrelation times at T_N for the production Monte Carlo (Plan B, referee M1).
# Same engine and settings as production (fef2_mc.py MC class; NEQ = 4000 equilibration sweeps), but after
# equilibration NSW = 16000 sweeps (= production length: 8000 measurements x 2 sweeps) with a measurement
# after EVERY sweep. tau_int (in sweeps) of |m|, m^2 and e: Sokal automatic window, W >= 6 tau.
# Production sampling (one measurement per 2 sweeps): independent samples per replica = NSW / (2 tau).
# Equilibration check: z = (mean|m| second half - first half) / its error (with tau); |z| < 3 expected.
# Job file lines:  name model slen D J1 J2 J3 J4 J7 T Ls nrep
# Usage: python3 fef2_tau.py JOBFILE OUTFILE ;  python3 fef2_tau.py sum OUTFILE
# ASCII only. numpy only.
import os, sys, zlib
import numpy as np
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fef2_mc import MC
NEQ = int(os.environ.get("NEQ", "4000")); NSW = int(os.environ.get("NSW", "16000")); NPROC = int(os.environ.get("NPROC", "56"))

def tau_int(x, c=6.0):
    x = np.asarray(x, float) - np.mean(x); n = len(x)
    if np.allclose(x, 0): return 0.5, 1
    f = np.fft.rfft(x, 2 * n); a = np.fft.irfft(f * np.conj(f))[:n]; a /= a[0]
    t = 0.5
    for W in range(1, n):
        t += a[W]
        if W >= c * t: return t, W
    return t, n

def job(a):
    name, model, slen, D, J, L, rep, T = a
    sim = MC(L, model, J, D, T, zlib.crc32(("TAU %s %s %d %d %.4f" % (name, model, L, rep, T)).encode()), slen)
    for i in range(NEQ): sim.sweep(adapt=(i < NEQ // 2))
    ms, es = [], []
    for i in range(NSW):
        sim.sweep(); m, e = sim.measure(); ms.append(abs(m)); es.append(e)
    ms = np.array(ms); es = np.array(es)
    ta, W = tau_int(ms); t2, _ = tau_int(ms ** 2); te, _ = tau_int(es)
    h = NSW // 2; err = np.std(ms, ddof=1) * np.sqrt(2 * ta / h) * np.sqrt(2)
    z = (ms[h:].mean() - ms[:h].mean()) / err if err > 0 else 0.0
    return "%-16s %-6s %3d %2d %9.4f %9.2f %9.2f %9.2f %6d %8.1f %6.2f" % (name, model, L, rep, T, ta, t2, te, W, NSW / (2 * ta), z)

def summary(fn):
    R = {}
    for l in open(fn):
        if l.startswith("#"): continue
        w = l.split(); R.setdefault((w[0], w[1]), {}).setdefault(int(w[2]), []).append([float(x) for x in w[5:11]])
    print("# set model L  n  tau|m|  tau_m2  tau_e (sweeps, mean over replicas)  Neff/replica(prod)  max|z_eq|")
    for (s, m), d in sorted(R.items()):
        Ls = sorted(d); lt = []
        for L in Ls:
            a = np.array(d[L]); lt.append(np.log(a[:, 0].mean()))
            print("%-16s %-6s %3d %2d %8.1f %8.1f %8.1f %10.1f %8.2f" % (s, m, L, len(a), a[:, 0].mean(), a[:, 1].mean(), a[:, 2].mean(), a[:, 4].mean(), np.abs(a[:, 5]).max()))
        if len(Ls) > 1:
            z = np.polyfit(np.log(Ls), lt, 1)[0]; print("%-16s %-6s dynamic exponent z (tau|m| ~ L^z): %.2f" % (s, m, z))

if __name__ == "__main__":
    if sys.argv[1] == "sum": summary(sys.argv[2]); sys.exit(0)
    jf, out = sys.argv[1], sys.argv[2]; tasks = []
    for l in open(jf):
        w = l.split()
        if not w or l.startswith("#"): continue
        J = dict(zip(("J1", "J2", "J3", "J4", "J7"), map(float, w[4:9])))
        for L in [int(x) for x in w[10].split(",")]:
            for rep in range(int(w[11])): tasks.append((w[0], w[1], float(w[2]), float(w[3]), J, L, rep, float(w[9])))
    dn = set()
    if os.path.exists(out):
        for l in open(out):
            if not l.startswith("#"): w = l.split(); dn.add((w[0], w[1], int(w[2]), int(w[3])))
    else: open(out, "w").write("# set model L rep T tau|m| tau_m2 tau_e window Neff_prod z_eq   (NEQ %d, NSW %d)\n" % (NEQ, NSW))
    todo = sorted([t for t in tasks if (t[0], t[1], t[5], t[6]) not in dn], key=lambda t: -t[5])
    print("%s: %d runs, %d to do" % (jf, len(tasks), len(todo)), flush=True)
    with Pool(NPROC) as p, open(out, "a") as g:
        for line in p.imap_unordered(job, todo): g.write(line + "\n"); g.flush()
    print("### TAU TAMAM", flush=True); summary(out)
