# fef2_jsub.py - least-squares exchange fit on a chosen set of configurations (Plan-B control calculations).
# Run inside a config directory (configs.dat, NAME.in, NAME.out, orb_target.dat); parse/sublat come from
# ../fef2_fit.py. A configuration enters only if converged and every Fe passes the production orbital criterion
# (|overlap of the minority eigenvector with orb_target.dat| > 0.95, as fef2_orbfix.py check).
# Convention as production: E/Fe = E0 + sum_n c_n (J_n S^2), c_n = bonds per Fe x sigma_i sigma_j (configs.dat);
# J(S=2) = J S^2 / 4; J > 0 AFM.
#   python3 ../fef2_jsub.py B   [names]   E0 + J1, J2, J3, J4, J7   (production set B; in these cells J2 = J2 + J5)
#   python3 ../fef2_jsub.py 113 [names]   1x1x3 cell: E0 + X, J2, J5, where X = J1 + 4 J4 + J6 (aliased in-cell)
#   python3 ../fef2_jsub.py FA  [FM AFM]  FM - AFM energy (meV/Fe; default c111_uu c111_ud) and (J2+J5) = /(8 S^2)
# Default names: every config in configs.dat. Prints J S^2, J(S=2), leave-one-out std, residuals, max |TS|.
# ASCII only. numpy only.
import sys, os, re
import numpy as np
RY = 13605.693122994
src = open("../fef2_fit.py").read(); exec(src[src.index("def parse"):src.index("rows = []")])
COLS = {"B": (("J1", [0]), ("J2", [1]), ("J3", [2]), ("J4", [3]), ("J7", [6])),
        "113": (("X=J1+4J4+J6", [0]), ("J2", [1]), ("J5", [4]))}

def load(names):
    cf = {l.split()[0]: (int(l.split()[1]), np.array(l.split()[2:9], float)) for l in open("configs.dat") if not l.startswith("#")}
    ref = {int(l.split()[1]): np.array(l.split()[2:7], float) for l in open("orb_target.dat") if l.startswith("sub")}
    out, bad = [], []
    for n in (names or list(cf)):
        if n not in cf or not os.path.exists(n + ".out"): bad.append((n, "missing")); continue
        r = parse(n + ".out"); sl = sublat(n + ".in")
        ok = r["conv"] and r["done"] and r["E"] is not None and len(r["atoms"]) == len(sl)
        ov = [abs(float(np.dot(a["minvec"], ref[s]))) for a, s in zip(r["atoms"], sl)] if ok else []
        if not ok or min(ov) <= 0.95: bad.append((n, "orbital/convergence")); continue
        out.append((n, cf[n][0], cf[n][1], r["E"] / cf[n][0] * RY, r["TS"] * RY if r["TS"] is not None else float("nan"), min(ov), [abs(a["mag"]) for a in r["atoms"]]))
    return out, bad

def lsq(X, y):
    c, *_ = np.linalg.lstsq(X, y, rcond=None); return c, y - X @ c

if __name__ == "__main__":
    mode = sys.argv[1]; names = sys.argv[2:]
    if mode == "FA":
        fm, afm = (names + ["c111_uu", "c111_ud"])[:2] if len(names) < 2 else names[:2]
        rows, bad = load([fm, afm])
        if len(rows) < 2: print("### FA: eksik/basarisiz:", bad); sys.exit(1)
        d = rows[0][3] - rows[1][3]
        print("FA %s - %s = %.4f meV/Fe   (J2+J5) S^2 = %.4f meV   J2+J5 (S=2) = %.5f meV   TS(FM)-TS(AFM) = %.2e meV/cell" % (
            fm, afm, d, d / 8, d / 32, rows[0][4] - rows[1][4])); sys.exit(0)
    rows, bad = load(names); cols = COLS[mode]
    for n, why in bad: print("skip %-24s %s" % (n, why))
    X = np.array([[1.0] + [sum(r[2][i] for i in ix) for lab, ix in cols] for r in rows]); y = np.array([r[3] for r in rows]); y0 = y.mean()
    p = X.shape[1]; sv = np.linalg.svd(X[:, 1:] - X[:, 1:].mean(0), compute_uv=False) if len(rows) > 1 else [0]
    if len(rows) < p or min(sv) < 1e-6: print("### %s: N=%d < %d parametre veya rank eksik (smin %.2e)" % (mode, len(rows), p, min(sv))); sys.exit(1)
    c, res = lsq(X, y - y0)
    loo = []
    if len(rows) > p:
        for i in range(len(rows)):
            m = np.arange(len(rows)) != i
            if np.linalg.matrix_rank(X[m]) == p: loo.append(lsq(X[m], y[m] - y0)[0])
    sd = np.std(loo, axis=0) if len(loo) > 2 else np.full(p, np.nan)
    print("# fef2_jsub %s  dir %s  N=%d  params=%d  smin=%.3f  rms=%.4f max=%.4f meV/Fe  min overlap %.3f  |m| %.2f-%.2f  max|TS| %.2e meV/cell" % (
        mode, os.path.basename(os.getcwd()), len(rows), p, min(sv), np.sqrt((res ** 2).mean()), np.abs(res).max(),
        min(r[5] for r in rows), min(min(r[6]) for r in rows), max(max(r[6]) for r in rows), max(abs(r[4]) for r in rows)))
    for k, (lab, ix) in enumerate(cols):
        print("  %-12s J*S^2 = %9.5f +- %.5f meV    J(S=2) = %9.5f meV" % (lab, c[k + 1], sd[k + 1], c[k + 1] / 4))
    for r, e in zip(rows, res): print("  res %-24s %8.4f meV/Fe" % (r[0], e))
