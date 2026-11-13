# fef2_conv.py - convergence summary (arg: dir, default S0; CV3 = target orbital, v3 geometry): relaxed u, cutoff/k convergence of E(FM)-E(AFM), commensurate-k check. ASCII only.
import re, os, sys
RY = 13605.693122994
def E(fn):
    if not os.path.exists(fn): return None
    t = open(fn).read(); e = re.findall(r"^!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", t, re.M)
    return float(e[-1]) if (e and "JOB DONE" in t and "convergence has been achieved" in t) else None
d = sys.argv[1] if len(sys.argv) > 1 else "S0"
uf = "u_relaxed_v3.txt" if d != "S0" else "u_relaxed.txt"
print("dir:", d, " u(F):", open(uf).read().strip() if os.path.exists(uf) else "-")
print("%-10s %16s %16s %14s" % ("set", "E_AFM (Ry)", "E_FM (Ry)", "FM-AFM meV/Fe"))
for s in ("e70", "e90", "e110", "k8812"):
    a, f = E("%s/cv_%s_afm.out" % (d, s)), E("%s/cv_%s_fm.out" % (d, s))
    if a is None or f is None: print("%-10s  eksik" % s); continue
    print("%-10s %16.8f %16.8f %14.4f" % (s, a, f, (f - a) / 2 * RY))
a1, f1 = E("%s/cv_e90_afm.out" % d), E("%s/cv_e90_fm.out" % d)
a8, f8 = E("%s/chk212_afm.out" % d), E("%s/chk212_fm.out" % d)
if None not in (a1, f1, a8, f8):
    print("commensurate-k check (per Fe, meV): AFM %.4f  FM %.4f  (2x1x2 vs 1x1x1, should be ~0)" %
          ((a8 / 8 - a1 / 2) * RY, (f8 / 8 - f1 / 2) * RY))
else:
    print("commensurate-k check: eksik")
a5, a4 = E("%s/cv_e90_afm.out" % d), E("%s/cv_orb4_afm.out" % d)
if None not in (a5, a4):
    print("orbital check: E(alt orbital) - E(standard) = %.2f meV/Fe  (must be > 0)" % ((a4 - a5) / 2 * RY))
else:
    print("orbital check: eksik")
