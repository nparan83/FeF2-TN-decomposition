# fef2_tabB.py - Plan-B result tables from the FSS outputs (run inside mcB/ after the fss_*.txt files exist).
# Inputs: fss_bc.txt (production Blume-Capel, all sizes), fss_bc1624.txt and fss_syn1624.txt (L = 16/24 only),
#   fss_syn_all.txt (synthetic sets, L = 12/16/24), fss_d0.txt, fss_val.txt, fss_heis082.txt, fss_heisDC.txt,
#   fef2_eft.txt, fef2_eftB.txt, fef2_eft.py (EFT code); ../atom/v7a_U*/jsub_B.txt, ../v7_113/jsub_113.txt,
#   ../v3_U4.51/jsub_B_sub8.txt (DFT controls). v2 (10.11.2026): synthetic sets from fss_syn_all.txt (two of seven
#   have too few points at L = 16/24 only); sections (7) atomic projectors, (8) J6 estimate added.
# Prints: (1) Table 2 with FSS T_N and ratios; (2) decomposition with MC error and neutron-parameter error
#   (sensitivities from the EFT code, checked against MC for J2 + sigma and D + sigma); (3) EFT/MC and MFA/MC of
#   the synthetic sets vs the production sets (same sizes); (4) isotropic Heisenberg (protocol of Ref. [3]);
#   (5) bcc Ising validation; (6) Heisenberg with D and with the El Mendili-Zhitomirsky D_C.
# Missing inputs are skipped. ASCII only. numpy only.
import os
import numpy as np
TEXP = 78.418                                    # Belanger & Yoshizawa 1987, Table IV
SIG = {"J1": 0.003, "J2": 0.02, "J4": 0.001, "J7": 0.0035, "D": 0.04}   # Sears et al. Table I (J7: (J7a+J7b)/2 -> 0.007/2)
def fss(fn):
    out = {}
    if not os.path.exists(fn): return out
    for l in open(fn):
        if l.startswith("#"): continue
        w = l.split()
        if len(w) >= 7:
            try: out[w[0]] = (float(w[3]), float(w[6]))
            except ValueError: pass
    return out
def table(fn, ci):
    out = {}
    if not os.path.exists(fn): return out
    for l in open(fn):
        if l.startswith("#") or not l.strip(): continue
        w = l.split()
        try: out[w[0]] = [float(w[i]) for i in ci]
        except (ValueError, IndexError): pass
    return out
BC, B16, S16, SA = fss("fss_bc.txt"), fss("fss_bc1624.txt"), fss("fss_syn1624.txt"), fss("fss_syn_all.txt")
EFT = table("fef2_eft.txt", (2, 3)); EB = table("fef2_eftB.txt", (9, 10))
print("== (1) Table 2: Blume-Capel S=2, D = 0.82 meV; MC = FSS fit, L = 16/24/32 (total error)")
print("%-16s %8s %8s %16s %8s %8s" % ("set", "MFA", "EFT", "MC", "MFA/MC", "EFT/MC"))
for s in ("neutron", "dft_U3.0_B14", "dft_U4.0_B14", "dft_U4.51_B14", "dft_U5.0_B14", "dft_U6.0_B14", "dft_U7.0_B14"):
    if s in BC and s in EFT:
        m, e = EFT[s]; t, d = BC[s]; print("%-16s %8.2f %8.2f %9.2f +- %4.2f %8.4f %8.4f" % (s, m, e, t, d, m / t, e / t))
if "neutron" in BC and "dft_U4.51_B14" in BC and "dft_U4.51_B14" in EFT:
    tn, dn = BC["neutron"]; tl, dl = BC["dft_U4.51_B14"]; te = EFT["dft_U4.51_B14"][1]
    fJ, fE, fM = tl / tn, te / tl, tn / TEXP; tot = te / TEXP
    eJ = fJ * np.hypot(dl / tl, dn / tn); eE = fE * dl / tl; eM = fM * dn / tn
    print("\n== (2) Decomposition  T_EFT(J_DFT)/T_exp = %.4f = J %.4f x EFT %.4f x model %.4f" % (tot, fJ, fE, fM))
    lt = np.log(tot); print("   log shares: J %.1f %%  EFT %.1f %%  model %.1f %%" % (100 * np.log(fJ) / lt, 100 * np.log(fE) / lt, 100 * np.log(fM) / lt))
    print("   MC errors only:   J %.4f +- %.4f   EFT %.4f +- %.4f   model %.4f +- %.4f" % (fJ, eJ, fE, eE, fM, eM))
    try:
        exec(open("fef2_eft.py").read().split('\nif __name__ == "__main__":')[0])
        NE = {"J1": -0.030, "J2": 0.46, "J3": 0.0, "J4": 0.005, "J7": 0.011}; D0 = 0.82; t0 = tn_eft(NE, D0, tol=0.0005); q = 0.0
        print("   neutron-parameter sensitivities (EFT code, +-1 sigma, central differences):")
        for k, sg in SIG.items():
            Jp, Jm, Dp, Dm = dict(NE), dict(NE), D0, D0
            if k == "D": Dp, Dm = D0 + sg, D0 - sg
            else: Jp[k] += sg; Jm[k] -= sg
            d = (tn_eft(Jp, Dp, tol=0.0005) - tn_eft(Jm, Dm, tol=0.0005)) / 2 / t0; q += d * d
            print("     %-3s sigma %.4f -> dT/T = %+.4f" % (k, sg, d))
        rn = np.sqrt(q); print("   total relative neutron error of T_N(J_n): %.4f" % rn)
        print("   with neutron error:  J %.3f +- %.3f   model %.3f +- %.3f   (EFT factor unaffected)" % (fJ, fJ * np.hypot(eJ / fJ, rn), fM, fM * np.hypot(eM / fM, rn)))
        sJ = np.log(fJ) / lt; sM = np.log(fM) / lt; dl_ = rn / lt
        print("   log shares with neutron error: J %.0f +- %.0f %%, model %.0f +- %.0f %%" % (100 * sJ, 100 * dl_, 100 * sM, 100 * dl_))
        print("   alternative ordering  [EFT(J_DFT)/EFT(J_n)] x [EFT(J_n)/MC(J_n)] x model: %.4f x %.4f x %.4f" % (te / EFT["neutron"][1], EFT["neutron"][1] / tn, fM))
        if "neu_J2p" in S16 and "neu_Dp" in S16 and "neutron" in B16:
            b = B16["neutron"][0]
            for s, k in (("neu_J2p", "J2"), ("neu_Dp", "D")):
                mc = S16[s][0] / b - 1; ef = EB[s][1] / EFT["neutron"][1] - 1
                print("   MC check (L=16/24): %-8s +1 sigma %s: MC dT/T = %+.4f +- %.4f   EFT dT/T = %+.4f" % (s, k, mc, np.hypot(S16[s][1] / S16[s][0], B16["neutron"][1] / b), ef))
    except Exception as ex: print("   (sensitivity step skipped: %s)" % ex)
if SA:
    print("\n== (3) EFT/MC and MFA/MC: production sets (L = 16/24/32) and synthetic/sensitivity sets (L = 12/16/24)")
    rows = [(s, EFT[s][0], EFT[s][1], BC[s]) for s in BC if s in EFT] + [(s, EB[s][0], EB[s][1], SA[s]) for s in SA if s in EB]
    for s, m, e, (t, d) in rows: print("%-16s MC %8.2f +- %4.2f   MFA/MC %.4f   EFT/MC %.4f" % (s, t, d, m / t, e / t))
    r = [e / t for s, m, e, (t, d) in rows]; q = [m / t for s, m, e, (t, d) in rows]
    print("   all %d sets: EFT/MC %.4f - %.4f (mean %.4f)   MFA/MC %.4f - %.4f" % (len(r), min(r), max(r), np.mean(r), min(q), max(q)))
D0R = fss("fss_d0.txt")
if D0R:
    print("\n== (4) isotropic classical Heisenberg (D = 0, spin length 2), x (S+1)/S = 1.5 as in Ref. [3]")
    for s, (t, d) in D0R.items(): print("%-16s T_MC %7.2f +- %4.2f   x1.5 = %7.2f   vs exp %+6.1f %%" % (s, t, d, 1.5 * t, 100 * (1.5 * t / TEXP - 1)))
V = fss("fss_val.txt")
if "bcc_ising" in V:
    tc0 = 1.0 / (0.1573725 * 0.08617333262); t, d = V["bcc_ising"]
    print("\n== (5) bcc nn Ising validation: Tc(FSS) = %.3f +- %.3f K, series %.3f K -> deviation %+.3f %% (+- %.3f %%)" % (t, d, tc0, 100 * (t / tc0 - 1), 100 * d / tc0))
H1, H2 = fss("fss_heis082.txt"), fss("fss_heisDC.txt")
if H1 or H2:
    print("\n== (6) classical Heisenberg with anisotropy (FSS): D = 0.82 meV and D_C mapping (0.875 D); x1.5")
    for s in sorted(set(H1) | set(H2), key=lambda x: (x != "neutron", x)):
        a = H1.get(s); b = H2.get(s); bc = BC.get(s)
        print("%-16s D=0.82: %s   0.875D: %s   mapped/BC %s" % (s, "%.2f+-%.2f" % a if a else "-", "%.2f+-%.2f (x1.5 %.2f)" % (b[0], b[1], 1.5 * b[0]) if b else "-",
              "%.4f" % (1.5 * b[0] / bc[0]) if (b and bc) else "-"))
import glob, re
for f in sorted(glob.glob("../atom/v7a_U*/jsub_B.txt")):
    t = open(f).read(); J = {k: float(v) for k, v in re.findall(r"(J\d)\s+J\*S\^2 =.*?J\(S=2\) =\s+(-?[\d.]+)", t)}
    U = re.search(r"v7a_U([\d.]+)", f).group(1)
    try:
        exec(open("fef2_eft.py").read().split('\nif __name__ == "__main__":')[0])
        je, tm = tn_mfa(J, 0.82); te = tn_eft(J, 0.82, tol=0.0005)
        rr = np.mean([EFT[s][1] / BC[s][0] for s in BC if s in EFT])
        print("\n== (7) atomic projectors, U = %s eV: J(S=2) %s" % (U, " ".join("%s=%.5f" % kv for kv in sorted(J.items()))))
        print("   J2/J2(neutron) = %.3f   J2 S^2 U = %.2f meV eV   J_eff %.4f   T_MFA %.2f   T_EFT %.2f   T_EFT/<EFT/MC> = %.2f K (%+.1f %% vs exp; not a MC result)" % (
              J["J2"] / 0.46, 4 * J["J2"] * float(U), je, tm, te, te / rr, 100 * (te / rr / TEXP - 1)))
    except Exception as ex: print("(7) skipped:", ex)
try:
    t1 = open("../v7_113/jsub_113.txt").read(); t2 = open("../v3_U4.51/jsub_B_sub8.txt").read()
    g = lambda t, k: float(re.search(k + r"\s+J\*S\^2 =.*?J\(S=2\) =\s+(-?[\d.]+)", t).group(1))
    X, J1, J4, J2, J5 = g(t1, r"X=J1\+4J4\+J6"), g(t2, "J1"), g(t2, "J4"), g(t1, "J2"), g(t1, "J5")
    J6 = X - J1 - 4 * J4; je = EFT and None
    print("\n== (8) 1x1x3 cell: J2 = %.5f, J5 = %.5f meV (J5/J2 = %.4f); X = J1+4J4+J6 = %.5f; with J1, J4 of the 8-config set at U_LR:" % (J2, J5, J5 / J2, X))
    print("   J6 = %.5f meV (k-meshes 6x6x9 and 6x6x8 combined); change of J_eff = -2 J6 = %.4f meV" % (J6, -2 * J6))
except Exception as ex: print("(8) skipped:", ex)
