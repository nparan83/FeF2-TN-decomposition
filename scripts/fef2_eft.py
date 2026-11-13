# fef2_eft.py - effective-field theory (Honmura-Kaneyoshi differential-operator EFT with the standard
# decoupling approximation) for the FeF2 Blume-Capel model, same Hamiltonian as fef2_mc.py:
#   H = sum_bonds J_n S_i S_j - D sum_i S_i^2,  S in {-2..2},  J > 0 AFM,  two sublattices (bct).
# Exact identity:  <S_i^k> = < F_k(x_i) >,  x_i = -sum_j J_ij S_j,  F_k(x) = sum_s s^k e^{b(xs + D s^2)} / Z.
# EFT: <F_k(x_i)> = < prod_j exp(-J_ij S_j d/dy) > F_k(y)|_{y=0}, decoupled as prod_j < exp(-J_ij S_j d/dy) >.
# For a spin taking 5 values, exp(a S) = sum_{n=0..4} A_n(a) S^n, so each factor is fixed by the moments
# (m, q, m3, r) = (<S>, <S^2>, <S^3>, <S^4>). This is identical to averaging F_k over INDEPENDENT neighbour
# spins whose 5-point distribution reproduces those moments - evaluated here exactly by convolution.
# T_N: largest eigenvalue of the linearised odd-moment map (m, m3) at the paramagnetic fixed point = 1.
# Also: single-site MFA (as fef2_mcsum.py). Validation: spin-1/2 Ising (values +-1), square & sc.
# ASCII only. numpy only.
import sys, os, re
import numpy as np
KB = 0.08617333262
SHELL_Z = {"J1": (2, +1), "J2": (8, -1), "J3": (4, +1), "J4": (8, +1), "J7": (4, +1)}   # (z, same(+1)/other(-1) sublattice)

def dist_from_moments(vals, mom):
    # probabilities p(s) on the support 'vals' with sum p s^n = mom[n] (n = 0..len-1)
    V = np.vander(vals, len(vals), increasing=True).T
    return np.linalg.solve(V, mom)

def field_dist(vals, p_same, p_other, shells):
    # distribution of x = -sum_j J_j S_j over all neighbours (independent, as in the decoupling approximation).
    # EFT_GRID=0: exact enumeration (slow for many shells). Default: FFT convolution on a uniform grid of
    # spacing EFT_DX meV (each delta split linearly between two grid points: mean exact, spurious variance
    # <= z*dx^2/4 per shell; with dx = 2e-4 meV this is ~1e-7 of (kT)^2).
    if os.environ.get("EFT_GRID", "1") == "0":
        xs, ws = np.array([0.0]), np.array([1.0])
        for J, z, sgn in shells:
            p = p_same if sgn > 0 else p_other
            for _ in range(z):
                nx = (xs[:, None] - J * vals[None, :]).ravel(); nw = (ws[:, None] * p[None, :]).ravel()
                u, inv = np.unique(np.round(nx, 9), return_inverse=True)
                xs, ws = u, np.bincount(inv, weights=nw)
        return xs, ws
    dx = float(os.environ.get("EFT_DX", "2e-4"))
    R = sum(z * abs(J) * np.max(np.abs(vals)) for J, z, s in shells) + 10 * dx
    n = int(2 ** np.ceil(np.log2(2 * R / dx + 1)))
    x0 = -0.5 * n * dx
    spec = np.ones(n // 2 + 1, complex)
    for J, z, sgn in shells:
        p = p_same if sgn > 0 else p_other
        g = np.zeros(n); pos = (-J * vals - x0) / dx; i = np.floor(pos).astype(int); f = pos - i
        np.add.at(g, i, p * (1 - f)); np.add.at(g, i + 1, p * f)
        G = np.fft.rfft(g)
        # g is centred at -x0: remove the offset after each power (shift phase), then raise to z
        k = np.arange(n // 2 + 1); ph = np.exp(2j * np.pi * k * (-x0 / dx) / n)
        spec *= (G * ph) ** z
    w = np.fft.irfft(spec * np.exp(-2j * np.pi * np.arange(n // 2 + 1) * (-x0 / dx) / n), n)
    xs = x0 + dx * np.arange(n); m = np.abs(w) > 1e-15
    return xs[m], w[m]

def Fk(vals, x, beta, D, kmax):
    e = beta * (x[:, None] * vals[None, :] + D * vals[None, :] ** 2)
    e -= e.max(axis=1, keepdims=True); w = np.exp(e); Z = w.sum(axis=1)
    return [(w * vals[None, :] ** k).sum(axis=1) / Z for k in range(1, kmax + 1)]

def eft_map(mom_odd, mom_even, vals, T, D, shells):
    # moments of sublattice A; sublattice B = A with odd moments reversed (collinear AFM)
    n = len(vals); mom = np.zeros(n); mom[0] = 1.0
    for k in range(1, n):
        mom[k] = mom_odd[(k - 1) // 2] if k % 2 else mom_even[k // 2 - 1]
    momB = mom * np.array([(-1) ** k for k in range(n)])
    pA, pB = dist_from_moments(vals, mom), dist_from_moments(vals, momB)
    xs, ws = field_dist(vals, pA, pB, shells)
    F = Fk(vals, xs, 1.0 / (KB * T), D, n - 1)
    out = [float((ws * f).sum()) for f in F]
    return np.array(out[0::2]), np.array(out[1::2])         # odd, even

_EV = {}
def para(vals, T, D, shells, it=200):
    # even moments at the paramagnetic fixed point (plain iteration, warm start from the last call)
    n = len(vals); ne = (n - 1) // 2
    if ne == 0: return np.zeros(0)
    ev = _EV.get(len(vals), np.array([np.mean(vals ** (2 * (k + 1))) for k in range(ne)]))
    for _ in range(it):
        _, ev2 = eft_map(np.zeros(len(vals) // 2), ev, vals, T, D, shells)
        if np.max(np.abs(ev2 - ev)) < 1e-11: ev = ev2; break
        ev = ev2
    _EV[len(vals)] = ev
    return ev

def lam(vals, T, D, shells, h=1e-6):
    no = len(vals) // 2; ev = para(vals, T, D, shells)
    Jm = np.zeros((no, no))
    for j in range(no):
        d = np.zeros(no); d[j] = h
        op, _ = eft_map(d, ev, vals, T, D, shells); om, _ = eft_map(-d, ev, vals, T, D, shells)
        Jm[:, j] = (op - om) / (2 * h)
    return max(abs(np.linalg.eigvals(Jm)))

def tn_eft(J, D, vals=np.arange(-2.0, 3.0), lo=None, hi=None, tol=0.002):
    shells = [(J[k], z, s) for k, (z, s) in SHELL_Z.items() if J.get(k, 0.0) != 0.0]
    f = lambda T: lam(vals, T, D, shells) - 1.0
    tm = tn_mfa(J, D)[1]
    a, b = lo or 0.6 * tm, hi or 1.02 * tm          # EFT lies below single-site MFA
    fa, fb = f(a), f(b)
    if fa < 0 or fb > 0: return None
    while b - a > tol:                                # regula falsi (Illinois) with bisection safeguard
        c = b - fb * (b - a) / (fb - fa)
        if not (a < c < b) or min(c - a, b - c) < 1e-3 * (b - a): c = 0.5 * (a + b)
        fc = f(c)
        if fc > 0: a, fa = c, fc; fb *= 0.5
        else: b, fb = c, fc; fa *= 0.5
    return 0.5 * (a + b)

def tn_mfa(J, D):
    Je = 8 * J["J2"] - 2 * J.get("J1", 0) - 4 * J.get("J3", 0) - 8 * J.get("J4", 0) - 4 * J.get("J7", 0)
    s = np.arange(-2, 3.)
    T = 100.0
    for _ in range(300):
        w = np.exp(D * s ** 2 / (KB * T)); T = Je * (s ** 2 * w).sum() / w.sum() / KB
    return Je, T

def test():
    # spin-1/2 Ising, H = -J sum s s (s = +-1): EFT kTc/J = 3.0898 (square), 5.0733 (sc) [Honmura & Kaneyoshi]
    v = np.array([-1.0, 1.0])
    for name, z in (("square", 4), ("sc", 6)):
        sh = [(-KB, z, +1)]         # ferro on one sublattice: same-sublattice shell, J = -1 (units kB=1)
        f = lambda T: lam(v, T, 0.0, sh) - 1.0
        a, b = 0.5, 20.0
        for _ in range(60):
            c = 0.5 * (a + b)
            if f(c) > 0: a = c
            else: b = c
        print("Ising s=1/2 %-6s z=%d  EFT kTc/J = %.4f" % (name, z, 0.5 * (a + b)))

def sets():
    S = {"neutron": {"J1": -0.030, "J2": 0.46, "J3": 0.0, "J4": 0.005, "J7": 0.011}}
    if os.path.exists("fef2_J.dat"):
        by = {}
        for l in open("fef2_J.dat"):
            if l.startswith("#"): continue
            w = l.split(); J = {k: float(v) / 4.0 for k, v in re.findall(r"(J\d)=(-?[\d.]+)\(", l)}
            if w[0] not in by or w[1].startswith("B"): by[w[0]] = (w[1][0] + w[2], J)
        for U, (sn, J) in sorted(by.items(), key=lambda x: float(x[0])): S["dft_U%s_%s" % (U, sn)] = J
    return S

if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "test": test(); sys.exit(0)
    D = float(os.environ.get("EFT_D", "0.82"))
    print("# Blume-Capel S=2, D = %.3f meV: single-site MFA and EFT (decoupling) Neel temperatures" % D)
    print("%-16s %8s %8s %8s" % ("set", "J_eff", "T_MFA", "T_EFT"))
    for n, J in sets().items():
        Je, tm = tn_mfa(J, D); te = tn_eft(J, D)
        print("%-16s %8.4f %8.2f %8s" % (n, Je, tm, "%.2f" % te if te else "-"), flush=True)
