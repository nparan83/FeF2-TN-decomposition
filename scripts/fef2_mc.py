# fef2_mc.py - classical Monte Carlo for FeF2 (rutile; magnetic Fe on a bct lattice, 2 sites/cell)
# H = sum_bonds J_n S_i.S_j - D sum_i (S_i^z)^2   (J > 0 AFM; same convention as Sears et al. 2026 and fef2_fit.py)
# models: bc     S^z in {-2,-1,0,1,2} (Blume-Capel, S = 2), heat-bath
#         ising1 S^z in {-1,1} (test only), heat-bath
#         heis   classical vector, |S| = SLEN (default 2), Metropolis with adaptive step
# shells (as fef2_gen.py): J1 (0,0,+-1) same sublattice z=2 | J2 body diagonal, other sublattice z=8 |
#   J3 (+-1,0,0),(0,+-1,0) z=4 | J4 (+-1,0,+-1),(0,+-1,+-1) z=8 | J7 (+-1,+-1,0) z=4.  J5, J6 = 0 (not resolved).
# Lattice L x L x L cells (L even), periodic. 16-colour parallel update (sublattice x parity of i,j,k).
# Units: J, D in meV, T in K. Order parameter: staggered z magnetisation per site / SLEN (full vector if D = 0).
# ASCII only. numpy only.
import sys, os, time
import numpy as np
KB = 0.08617333262
SAME = {"J1": [(0, 0, 1), (0, 0, -1)],
        "J3": [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)],
        "J4": [(a, 0, c) for a in (1, -1) for c in (1, -1)] + [(0, b, c) for b in (1, -1) for c in (1, -1)],
        "J7": [(a, b, 0) for a in (1, -1) for b in (1, -1)]}
J2OFF = {0: [(a, b, c) for a in (0, -1) for b in (0, -1) for c in (0, -1)],
         1: [(a, b, c) for a in (0, 1) for b in (0, 1) for c in (0, 1)]}
def shifted(a, o):        # value of the neighbour at cell offset o (last three axes are i, j, k)
    return np.roll(a, shift=(-o[0], -o[1], -o[2]), axis=(-3, -2, -1))
def field(S, J, b):       # local exchange field on sublattice b; S shape (2,L,L,L) or (3,2,L,L,L)
    sb = S[..., b, :, :, :]; so = S[..., 1 - b, :, :, :]
    h = np.zeros_like(sb)
    for k, offs in SAME.items():
        if J.get(k, 0.0) != 0.0:
            for o in offs: h += J[k] * shifted(sb, o)
    for o in J2OFF[b]: h += J["J2"] * shifted(so, o)
    return h
class MC:
    def __init__(s, L, model, J, D, T, seed, slen=2.0):
        s.L, s.model, s.J, s.D, s.beta = L, model, J, D, 1.0 / (KB * T)
        s.rng = np.random.default_rng(seed)
        ijk = np.indices((L, L, L)) % 2
        s.masks = [(ijk[0] == p) & (ijk[1] == q) & (ijk[2] == r) for p in (0, 1) for q in (0, 1) for r in (0, 1)]
        if model in ("bc", "ising1"):
            s.vals = np.array([-2., -1., 0., 1., 2.]) if model == "bc" else np.array([-1., 1.])
            s.slen = s.vals.max()
            s.S = np.empty((2, L, L, L)); s.S[0] = s.slen; s.S[1] = -s.slen      # AFM start
        else:
            s.slen = slen; s.S = np.zeros((3, 2, L, L, L)); s.S[2, 0] = slen; s.S[2, 1] = -slen
            s.step = 0.5 * slen; s.acc = [0, 0]
    def sweep(s, adapt=False):
        for b in (0, 1):
            for m in s.masks:
                h = field(s.S, s.J, b)
                if s.model != "heis":
                    hm = h[m]
                    e = hm[:, None] * s.vals[None, :] - s.D * s.vals[None, :] ** 2
                    w = np.exp(-s.beta * (e - e.min(axis=1, keepdims=True)))
                    c = np.cumsum(w, axis=1); u = s.rng.random(len(hm)) * c[:, -1]
                    s.S[b][m] = s.vals[(c < u[:, None]).sum(axis=1)]
                else:
                    old = s.S[:, b][:, m]; hm = h[:, m]
                    new = old + s.step * s.rng.standard_normal(old.shape)
                    new *= s.slen / np.linalg.norm(new, axis=0)
                    dE = ((new - old) * hm).sum(axis=0) - s.D * (new[2] ** 2 - old[2] ** 2)
                    ok = s.rng.random(dE.shape) < np.exp(-s.beta * np.clip(dE, 0, None))
                    sub = s.S[:, b]; blk = sub[:, m]; blk[:, ok] = new[:, ok]; sub[:, m] = blk
                    s.acc[0] += ok.sum(); s.acc[1] += ok.size
        if adapt and s.model == "heis" and s.acc[1]:
            r = s.acc[0] / s.acc[1]; s.step = float(np.clip(s.step * (1.1 if r > 0.5 else 0.9), 0.02 * s.slen, 2 * s.slen)); s.acc = [0, 0]
    def measure(s):
        N = 2 * s.L ** 3
        if s.model != "heis":
            ms = (s.S[0].sum() - s.S[1].sum()) / (N * s.slen)
            e = 0.5 * sum((s.S[b] * field(s.S, s.J, b)).sum() for b in (0, 1)) - s.D * (s.S ** 2).sum()
        else:
            if s.D == 0.0:   # isotropic: full staggered vector (z-only is not an order parameter without anisotropy)
                ms = np.linalg.norm(s.S[:, 0].sum(axis=(1, 2, 3)) - s.S[:, 1].sum(axis=(1, 2, 3))) / (N * s.slen)
            else:
                ms = (s.S[2, 0].sum() - s.S[2, 1].sum()) / (N * s.slen)
            e = 0.5 * sum((s.S[:, b] * field(s.S, s.J, b)).sum() for b in (0, 1)) - s.D * (s.S[2] ** 2).sum()
        return ms, e / N
def run(L, model, J, D, T, neq, nmeas, seed, slen=2.0, every=2):
    sim = MC(L, model, J, D, T, seed, slen)
    for i in range(neq): sim.sweep(adapt=(i < neq // 2))
    ms, es = [], []
    for i in range(nmeas):
        for _ in range(every): sim.sweep()
        m, e = sim.measure(); ms.append(m); es.append(e)
    ms = np.abs(np.array(ms)); es = np.array(es); N = 2 * L ** 3
    m2, m4 = (ms ** 2).mean(), (ms ** 4).mean()
    out = {"m": ms.mean(), "m2": m2, "U4": 1 - m4 / (3 * m2 ** 2) if m2 > 0 else 0.0, "e": es.mean(),
           "C": N * es.var() / (KB * T ** 2), "chi": N * (m2 - ms.mean() ** 2) / (KB * T),
           "acc": (sim.acc[0] / sim.acc[1] if model == "heis" and sim.acc[1] else 1.0)}
    nb = 10; bl = np.array_split(ms, nb)      # binning error of U4
    u = [1 - (x ** 4).mean() / (3 * (x ** 2).mean() ** 2) for x in bl if len(x)]
    out["dU4"] = np.std(u) / np.sqrt(len(u)); return out
if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "test":
    # bcc nearest-neighbour checks (J2 only, D = 0): Ising +-1 Tc = 6.3515 J/kB; classical Heisenberg |S|=1 Tc = 2.054 J/kB
    for model, Tc in ((("heis", 2.054 / KB),) if "heis" in sys.argv else (("ising1", 6.3515 / KB), ("heis", 2.054 / KB))):
        for T in (0.97 * Tc, 1.0 * Tc, 1.03 * Tc):
            r = [run(L, model, {"J2": 1.0}, 0.0, T, 300, 600, 7, slen=1.0) for L in (6, 10)]
            print("%-6s T/Tc=%.2f  U4(L=6)=%.3f+-%.3f  U4(L=10)=%.3f+-%.3f" % (model, T / Tc, r[0]["U4"], r[0]["dU4"], r[1]["U4"], r[1]["dU4"]), flush=True)
