# fef2_fss.py - finite-size-scaling estimate of T_N from Binder-cumulant data (Plan B; replaces the
# first-sign-change estimator of fef2_mcfine.py sum, which is fragile where U4(L2)-U4(L1) is within noise).
# Model: U4(T, L) = P_k(x),  x = (T - Tc) L^(1/nu) / T0, all sizes fitted together, nu fixed:
#   3D Ising nu = 0.6300 (bc, ising1, heis with D > 0)  |  3D Heisenberg nu = 0.7112 (heis with D = 0: --nu 0.7112)
# Only points with lo < U4 < hi enter (critical region). Weights: replica standard error (floor 2e-3).
# Tc: chi2 scan over T0*(1 +- 2.5 %) + golden-section refinement. U* = P_k(0).
# Statistical error: bootstrap over replicas (replicas resampled independently at each T, L), NB samples.
# Systematic error: largest shift among variants (k = 2, 3, 4 with the window narrowed/widened by 10-20 %;
#   default window: k=2 0.24-0.60, k=4 0.105-0.618, k=3 0.24-0.618, k=3 0.105-0.555) and the
#   two-largest-sizes-only fit. Total = quadrature sum.
# Input: any number of files with the fef2_mcfine.dat columns  set model L rep T m m2 U4 dU4 e C chi acc
# Usage: python3 fef2_fss.py MODEL file1.dat [file2.dat ...] [--nu 0.63] [--sets a,b] [--nb 200] [--Ls 16,24,32]
#                            [--lo 0.15 --hi 0.60]   (heis with D = 0, 3-component order parameter: U4 runs from
#                            4/9 (disordered) to 2/3 (ordered) -> use --nu 0.7112 --lo 0.50 --hi 0.655)
# ASCII only. numpy only.
import sys
import numpy as np
from collections import defaultdict
LO, HI, K = 0.15, 0.60, 3

def load(files, model, sets=None, Ls=None):
    R = defaultdict(lambda: defaultdict(list))          # (set, L) -> T -> [U4 per replica]
    for f in files:
        for l in open(f):
            if l.startswith("#"): continue
            w = l.split()
            if len(w) < 8 or w[1] != model: continue
            if sets and w[0] not in sets: continue
            L = int(w[2])
            if Ls and L not in Ls: continue
            R[(w[0], L)][round(float(w[4]), 4)].append(float(w[7]))
    return R

def points(R, s, Ls, rng=None):
    T, L, U, E = [], [], [], []
    for Lx in Ls:
        for t, v in R[(s, Lx)].items():
            v = np.array(v)
            if len(v) < 2: continue
            e = max(np.std(v, ddof=1) / np.sqrt(len(v)), 2e-3)
            if rng is not None: v = rng.choice(v, len(v))
            T.append(t); L.append(Lx); U.append(v.mean()); E.append(e)
    return np.array(T), np.array(L, float), np.array(U), np.array(E)

def chi2(tc, T, L, U, E, t0, nu, k):
    x = (T - tc) * L ** (1.0 / nu) / t0
    c = np.polyfit(x, U, k, w=1.0 / E)
    return (((np.polyval(c, x) - U) / E) ** 2).sum(), c

def fit(T, L, U, E, t0, nu, k=K, lo=LO, hi=HI):
    m = (U > lo) & (U < hi)
    if m.sum() < k + 3 or len(set(L[m])) < 2: return None
    T, L, U, E = T[m], L[m], U[m], E[m]
    grid = t0 * (1 + np.linspace(-0.025, 0.025, 201))
    c2 = [chi2(t, T, L, U, E, t0, nu, k)[0] for t in grid]
    i = int(np.argmin(c2)); a, b = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
    g = (np.sqrt(5) - 1) / 2
    for _ in range(40):
        x1, x2 = b - g * (b - a), a + g * (b - a)
        if chi2(x1, T, L, U, E, t0, nu, k)[0] < chi2(x2, T, L, U, E, t0, nu, k)[0]: b = x2
        else: a = x1
    tc = 0.5 * (a + b); cc, coef = chi2(tc, T, L, U, E, t0, nu, k)
    edge = not (grid[2] < tc < grid[-3])
    return tc, cc / max(1, len(U) - k - 2), len(U), np.polyval(coef, 0.0), edge

def guess(R, s, Ls, lo, hi):
    # T0: first sign change of U4(Lmax) - U4(Lmin) on heating; else median T of points in the window
    A, B = R[(s, Ls[0])], R[(s, Ls[-1])]; T = sorted(set(A) & set(B))
    d = [np.mean(B[t]) - np.mean(A[t]) for t in T]
    for i in range(len(T) - 1):
        if d[i] > 0 and d[i + 1] <= 0: return T[i] + (T[i + 1] - T[i]) * d[i] / (d[i] - d[i + 1])
    allT = [t for Lx in Ls for t, v in R[(s, Lx)].items() if lo < np.mean(v) < hi]
    return float(np.median(allT)) if allT else None

if __name__ == "__main__":
    a = sys.argv[1:]; model = a[0]; files = []; nu = 0.6300; sets = None; nb = 200; Ls = None; lo, hi = LO, HI; i = 1
    while i < len(a):
        if a[i] == "--nu": nu = float(a[i + 1]); i += 2
        elif a[i] == "--sets": sets = a[i + 1].split(","); i += 2
        elif a[i] == "--nb": nb = int(a[i + 1]); i += 2
        elif a[i] == "--Ls": Ls = [int(x) for x in a[i + 1].split(",")]; i += 2
        elif a[i] == "--lo": lo = float(a[i + 1]); i += 2
        elif a[i] == "--hi": hi = float(a[i + 1]); i += 2
        else: files.append(a[i]); i += 1
    R = load(files, model, sets, Ls)
    print("# fef2_fss.py  model=%s  nu=%.4f  window %.3f-%.3f  files=%s" % (model, nu, lo, hi, " ".join(files)))
    dd = hi - lo; VAR = ((2, lo + 0.2 * dd, hi), (4, lo - 0.1 * dd, hi + 0.04 * dd), (3, lo + 0.2 * dd, hi + 0.04 * dd), (3, lo - 0.1 * dd, hi - 0.1 * dd))
    print("# %-15s %-6s %-10s %9s %6s %6s %6s %6s %4s %6s %9s %s" % ("set", "model", "sizes", "Tc(K)", "stat", "syst", "total", "chi2/f", "n", "U*", "Tc(2L)", "flags"))
    for s in sorted({k[0] for k in R}, key=lambda x: (x != "neutron", x)):
        Lall = sorted({k[1] for k in R if k[0] == s})
        if len(Lall) < 2: continue
        t0 = guess(R, s, Lall, lo, hi)
        if t0 is None: print("%-17s no data in window" % s); continue
        T, L, U, E = points(R, s, Lall); r = fit(T, L, U, E, t0, nu, K, lo, hi)
        if r is None: print("%-17s too few points" % s); continue
        tc, cf, n, us, edge = r
        rng = np.random.default_rng(12345); bs = []
        for _ in range(nb):
            Tb, Lb, Ub, Eb = points(R, s, Lall, rng); rb = fit(Tb, Lb, Ub, Eb, t0, nu, K, lo, hi)
            if rb is not None: bs.append(rb[0])
        stat = float(np.std(bs)) if len(bs) > 10 else float("nan")
        var = [fit(T, L, U, E, t0, nu, k, l1, h1) for k, l1, h1 in VAR]
        T2, L2, U2, E2 = points(R, s, Lall[-2:]); r2 = fit(T2, L2, U2, E2, t0, nu, K, lo, hi)
        shifts = [abs(v[0] - tc) for v in var if v is not None] + ([abs(r2[0] - tc)] if r2 is not None else [])
        syst = max(shifts) if shifts else float("nan"); tot = float(np.hypot(stat, syst))
        flags = ("EDGE " if edge else "") + ("CHI2 " if cf > 4 else "")
        print("%-17s %-6s %-10s %9.3f %6.3f %6.3f %6.3f %6.2f %4d %6.3f %9s %s" % (s, model, "/".join(map(str, Lall)), tc, stat, syst, tot, cf, n, us,
              "%.3f" % r2[0] if r2 is not None else "-", flags), flush=True)
