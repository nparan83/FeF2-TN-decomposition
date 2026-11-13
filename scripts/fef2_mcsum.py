# fef2_mcsum.py - summary of fef2_mc.dat: Binder-cumulant crossings (pairs of L), chi peaks, and the
# single-site mean-field T_N of the same Hamiltonian (Blume-Capel S=2 and classical Heisenberg |S|=2 + D).
# MFA: kT_N = J_eff <(S^z)^2>_0(T_N),  J_eff = 8J2 - 2J1 - 4J3 - 8J4 - 4J7 (AFM order, J > 0 AFM).
# ASCII only. numpy only.
import re, os, sys
import numpy as np
exec(open("fef2_mcrun.py").read().split('\nif __name__ == "__main__":')[0])
KB = 0.08617333262; D = 0.82
def mfa(J, model):
    Je = 8 * J["J2"] - 2 * J.get("J1", 0) - 4 * J.get("J3", 0) - 8 * J.get("J4", 0) - 4 * J.get("J7", 0)
    def sz2(T):
        b = 1 / (KB * T)
        if model == "bc":
            s = np.arange(-2, 3.); w = np.exp(b * D * s ** 2); return (s ** 2 * w).sum() / w.sum()
        c = np.linspace(-1, 1, 4001); w = np.exp(b * D * 4 * c ** 2); return 4 * (c ** 2 * w).sum() / w.sum()
    T = 100.0
    for _ in range(200): T = Je * sz2(T) / KB
    return Je, T
R = read(); S = sets()
print("%-16s %-5s %8s %8s | %s" % ("set", "model", "J_eff", "T_MFA", "Binder crossings (K)   chi-peak T per L"))
for n, J in S.items():
    for m in ("bc", "heis"):
        Ls = sorted(L for (a, b, L) in R if a == n and b == m)
        Je, Tm = mfa(J, m); cr = []
        for L1, L2 in zip(Ls[:-1], Ls[1:]):
            a, b = R[(n, m, L1)], R[(n, m, L2)]
            T = np.intersect1d(np.round(a[:, 0], 3), np.round(b[:, 0], 3))
            if len(T) < 3: continue
            u1 = np.array([a[np.round(a[:, 0], 3) == t][0, 3] for t in T]); u2 = np.array([b[np.round(b[:, 0], 3) == t][0, 3] for t in T])
            d = u2 - u1; k = np.where(np.sign(d[:-1]) != np.sign(d[1:]))[0]
            for i in k[:1]: cr.append("%d/%d:%.1f" % (L1, L2, T[i] + (T[i + 1] - T[i]) * d[i] / (d[i] - d[i + 1])))
        pk = ["L%d:%.0f" % (L, R[(n, m, L)][np.argmax(R[(n, m, L)][:, 7]), 0]) for L in Ls]
        print("%-16s %-5s %8.4f %8.1f | %s   %s" % (n, m, Je, Tm, " ".join(cr) if cr else "-", " ".join(pk)))
