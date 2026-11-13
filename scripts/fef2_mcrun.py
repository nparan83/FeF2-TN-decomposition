# fef2_mcrun.py - MC driver for FeF2: parameter sets x models, stage 1 (coarse T scan, L=12) then
# stage 2 (L = 12, 16, 24 around the chi peak). Checkpointed: finished runs are skipped on restart.
# Sets: "neutron" (Sears et al., arXiv:2601.04303: J1=-0.030, J2=0.46, J4=0.005, (J7a+J7b)/2=0.011 meV) and
# one set per U line of fef2_J.dat (set B preferred, else A; file stores J_n*S^2 -> divided by 4 for S=2).
# Set name carries fit set + number of configs (dft_U5.0_B7): a refit with more configs gives a new name, no stale reuse.
# D = 0.82 meV for every set (neutron value; DFT here has no SOC -> experimental parameter, not fitted).
# Output fef2_mc.dat: set model L T <|m|> <m2> U4 dU4 <e> C chi acc. ASCII only. numpy only.
import os, re, sys, itertools, zlib
import numpy as np
from multiprocessing import Pool
sys.argv = sys.argv[:1] + [a for a in sys.argv[1:]]
exec(open("fef2_mc.py").read().split('if __name__ == "__main__"')[0])
QUICK = os.environ.get("QUICK"); NPROC = int(os.environ.get("NPROC", "56")); D = 0.82
OUT = "fef2_mc.dat"
def sets():
    S = {"neutron": {"J1": -0.030, "J2": 0.46, "J3": 0.0, "J4": 0.005, "J7": 0.011}}
    if os.path.exists("fef2_J.dat"):
        by = {}
        for l in open("fef2_J.dat"):
            if l.startswith("#"): continue
            w = l.split(); J = {k: float(v) / 4.0 for k, v in re.findall(r"(J\d)=(-?[\d.]+)\(", l)}
            J.setdefault("J7", 0.0)
            if w[0] not in by or w[1].startswith("B"): by[w[0]] = (w[1][0] + w[2], J)      # e.g. B14 = set B, 14 configs
        for U, (sn, J) in sorted(by.items(), key=lambda x: float(x[0])): S["dft_U%s_%s" % (U, sn)] = J
    return S
def done():
    d = set()
    if os.path.exists(OUT):
        for l in open(OUT):
            if not l.startswith("#"): w = l.split(); d.add((w[0], w[1], int(w[2]), "%.3f" % float(w[3])))
    return d
def job(a):
    name, model, L, T, J, neq, nm = a
    seed = zlib.crc32(("%s %s %d %.3f" % (name, model, L, T)).encode())   # deterministic (hash() is salted per process)
    r = run(L, model, J, D, T, neq, nm, seed)
    return "%-16s %-5s %3d %9.3f %.6f %.6f %.6f %.6f %.6f %.6e %.6e %.3f" % (name, model, L, T, r["m"], r["m2"], r["U4"], r["dU4"], r["e"], r["C"], r["chi"], r["acc"])
def go(tasks):
    dn = done(); todo = [t for t in tasks if (t[0], t[1], t[2], "%.3f" % t[3]) not in dn]
    print("runs: %d total, %d to do" % (len(tasks), len(todo)), flush=True)
    todo.sort(key=lambda t: -t[2])                       # big L first (load balance)
    with Pool(NPROC) as p, open(OUT, "a") as g:
        for line in p.imap_unordered(job, todo):
            g.write(line + "\n"); g.flush()
def read():
    R = {}
    for l in open(OUT):
        if l.startswith("#"): continue
        w = l.split(); R.setdefault((w[0], w[1], int(w[2])), []).append([float(x) for x in w[3:]])
    return {k: np.array(sorted(v)) for k, v in R.items()}
if __name__ == "__main__":
    S = sets(); models = ("bc", "heis")
    if not os.path.exists(OUT): open(OUT, "w").write("# set model L T m m2 U4 dU4 e C chi acc   (H = sum J S.S - D sum Sz^2, meV, K)\n")
    for k, v in S.items(): print(k, v)
    L1, Ts1 = (8, np.arange(30, 250, 20.)) if QUICK else (12, np.arange(30, 241, 10.))
    ne1, nm1 = (100, 200) if QUICK else (1000, 2000)
    go([(n, m, L1, T, J, ne1, nm1) for n, J in S.items() for m in models for T in Ts1])
    R = read(); tasks = []
    Ls = (8, 12) if QUICK else (12, 16, 24)
    ne2, nm2 = (150, 300) if QUICK else (2000, 4000)
    for n, J in S.items():
        for m in models:
            a = R[(n, m, L1)]; a = a[(a[:, 0] >= Ts1[0] - 1e-6) & (np.isin(np.round(a[:, 0], 3), np.round(Ts1, 3)))]
            Ts = a[np.argmax(a[:, 7]), 0]                    # chi peak of stage 1
            for L in Ls:
                for T in np.round(np.linspace(0.88 * Ts, 1.12 * Ts, 5 if QUICK else 13), 3): tasks.append((n, m, L, float(T), J, ne2, nm2))
    go(tasks)
    print("### MC TAMAM", flush=True)
