# fef2_mkjobsB.py - writes the Plan-B Monte Carlo job files from the production parameter sets (fef2_mcrun.py
# sets(): neutron + fef2_J.dat / 4, exact production values) and computes MFA/EFT for the synthetic test sets.
# Window centres T0 are search-window centres only (FSS re-analysis of the stage-3 data, 10.11.2026), not results.
# Writes: mcB_tau.txt mcB_focus.txt mcB_scan5.txt mcB_d0.txt mcB_val.txt mcB_syn.txt  and  fef2_eftB.txt
# Run inside mcB/ (needs fef2_mc.py fef2_mcrun.py fef2_eft.py fef2_J.dat). ASCII only. numpy only.
import numpy as np
exec(open("fef2_mcrun.py").read().split('\nif __name__ == "__main__":')[0])      # sets(), D (= 0.82)
E = open("fef2_eft.py").read().split('\nif __name__ == "__main__":')[0]; exec(E)  # tn_eft, tn_mfa, KB
S = sets(); D0 = 0.82
assert abs(D - D0) < 1e-12, "fef2_mcrun.py D must be 0.82 (production copy)"
def line(name, model, slen, Dv, J, T0, win, nT, Ls, nrep, rep0):
    return "%-16s %-6s %.1f %.4f %.7f %.7f %.7f %.7f %.7f %.4f %.4f %d %s %d %d\n" % (
        name, model, slen, Dv, J.get("J1", 0), J.get("J2", 0), J.get("J3", 0), J.get("J4", 0), J.get("J7", 0), T0, win, nT, Ls, nrep, rep0)
TC = {"neutron": 82.36, "dft_U3.0_B14": 149.72, "dft_U4.0_B14": 122.21, "dft_U4.51_B14": 110.48,
      "dft_U5.0_B14": 100.47, "dft_U6.0_B14": 82.91, "dft_U7.0_B14": 68.59}
HEIS = {"neutron": 54.59}
HD = "# name model slen D J1 J2 J3 J4 J7 T0 win nT Ls nrep rep0\n"
with open("mcB_tau.txt", "w") as g:      # T, Ls, nrep for fef2_tau.py
    g.write("# name model slen D J1 J2 J3 J4 J7 T Ls nrep\n")
    for n, m, Dv, T in (("neutron", "bc", D0, TC["neutron"]), ("dft_U4.51_B14", "bc", D0, TC["dft_U4.51_B14"]), ("neutron", "heis", D0, HEIS["neutron"])):
        J = S[n]; g.write("%-16s %-6s 2.0 %.4f %.7f %.7f %.7f %.7f %.7f %.4f 16,24,32 4\n" % (n, m, Dv, J["J1"], J["J2"], J.get("J3", 0), J["J4"], J["J7"], T))
with open("mcB_focus.txt", "w") as g:    # decomposition sets: focused grid, replicas 12-17 (stage 3 used 0-11)
    g.write(HD)
    for n in ("neutron", "dft_U4.51_B14"): g.write(line(n, "bc", 2.0, D0, S[n], TC[n], 0.008, 7, "16,24,32", 6, 12))
with open("mcB_scan5.txt", "w") as g:    # optional: U scan, focused grid, replicas 4-7 (stage 3 used 0-3)
    g.write(HD)
    for n in ("dft_U3.0_B14", "dft_U4.0_B14", "dft_U5.0_B14", "dft_U6.0_B14", "dft_U7.0_B14"): g.write(line(n, "bc", 2.0, D0, S[n], TC[n], 0.008, 7, "16,24,32", 4, 4))
with open("mcB_d0.txt", "w") as g:       # isotropic classical Heisenberg (protocol of Ref. [3]); centres from a local L=8/12 scan
    g.write(HD)
    g.write(line("neutron_D0", "heis", 2.0, 0.0, S["neutron"], 45.0, 0.04, 11, "16,24,32", 4, 0))
    g.write(line("dft_U4.51_D0", "heis", 2.0, 0.0, S["dft_U4.51_B14"], 62.0, 0.04, 11, "16,24,32", 4, 0))
with open("mcB_val.txt", "w") as g:      # bcc nearest-neighbour Ising (+-1): beta_c J = 0.1573725 (Butera & Comi) -> Tc = J/(0.1573725 kB)
    g.write(HD); g.write(line("bcc_ising", "ising1", 1.0, 0.0, {"J2": 1.0}, 1.0 / (0.1573725 * KB), 0.008, 9, "16,24,32", 6, 0))
LR = S["dft_U4.51_B14"]; NE = S["neutron"]
SYN = [("syn_J2only", {"J2": LR["J2"]}, D0), ("syn_J1neg", dict(LR, J1=-LR["J1"]), D0),
       ("syn_J47x2", dict(LR, J4=2 * LR["J4"], J7=2 * LR["J7"]), D0), ("syn_Dlow", LR, 0.20), ("syn_Dhigh", LR, 2.70),
       ("neu_J2p", dict(NE, J2=NE["J2"] + 0.02), D0), ("neu_Dp", NE, D0 + 0.04)]
with open("fef2_eftB.txt", "w") as h, open("mcB_syn.txt", "w") as g:
    h.write("# Blume-Capel S=2 synthetic/sensitivity sets: single-site MFA and EFT (fef2_eft.py, same code as fef2_eft.txt)\n")
    h.write("# set              D     J1       J2       J3       J4       J7      D/J2   J_eff    T_MFA    T_EFT\n")
    g.write(HD)
    for n, J, Dv in SYN:
        je, tm = tn_mfa(J, Dv); te = tn_eft(J, Dv, tol=0.0005)
        h.write("%-16s %5.3f %8.5f %8.5f %8.5f %8.5f %8.5f %6.3f %7.4f %8.2f %8.2f\n" % (n, Dv, J.get("J1", 0), J["J2"], J.get("J3", 0), J.get("J4", 0), J.get("J7", 0), Dv / J["J2"], je, tm, te))
        g.write(line(n, "bc", 2.0, Dv, J, te / 1.08, 0.025, 11, "12,16,24", 4, 0))
print(open("fef2_eftB.txt").read())
for f in ("mcB_tau.txt", "mcB_focus.txt", "mcB_scan5.txt", "mcB_d0.txt", "mcB_val.txt", "mcB_syn.txt"): print("==", f); print(open(f).read())
