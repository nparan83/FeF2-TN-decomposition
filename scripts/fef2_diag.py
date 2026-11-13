# fef2_diag.py - orbital-state diagnosis of FeF2 DFT+U outputs. ASCII only. Needs numpy.
# Minority-spin d orbital of every Fe -> SIGNED eigenvector, QE real-harmonic order (z2, xz, yz, x2-y2, xy).
# Body-centre Fe is mapped onto the corner one by the mirror y -> -y (yz, xy change sign); classes compared by |dot|.
# Usage: python3 fef2_diag.py U3.0 U4.0 ...   (also reads S0/ if present)
import sys, os, re, glob
import numpy as np
RY = 13605.693122994
src = open("fef2_fit.py").read()          # reuse parse() and sublat() from the fit script (same directory)
exec(src[:src.index("rows = []")])
classes = []
def label(v):
    for k, c in enumerate(classes):
        if abs(np.dot(v, c)) > 0.95: return "ABCDEFGHIJ"[k]
    classes.append(v); return "ABCDEFGHIJ"[len(classes) - 1]
def canon(v, s):
    v = v.copy()
    if s == 1: v[2], v[4] = -v[2], -v[4]
    return v
def summ(fo, fi):
    r = parse(fo); sl = sublat(fi)
    lab = "".join(label(canon(a["minvec"], s)) for a, s in zip(r["atoms"], sl)) if len(r["atoms"]) == len(sl) else "?"
    return r, lab
for d in sys.argv[1:] + (["S0"] if os.path.isdir("S0") else []):
    print("== %s" % d)
    cf = os.path.join(d, "configs.dat")
    keep = [l.split()[0] for l in open(cf) if not l.startswith("#")] if os.path.exists(cf) else None
    outs = sorted(glob.glob(os.path.join(d, "*.out")))
    if keep: outs = [os.path.join(d, n + ".out") for n in keep if os.path.exists(os.path.join(d, n + ".out"))]
    if keep: print("  configs: %d, outputs: %d" % (len(keep), len(outs)))
    rows = []
    for fo in outs:
        n = os.path.basename(fo)[:-4]; fi = os.path.join(d, n + ".in")
        if not os.path.exists(fi) or "_david_fail" in n: continue
        r, lab = summ(fo, fi)
        nfe = len(sublat(fi)); e = r["E"] / nfe * RY if r["E"] is not None else None
        rows.append((n, r["conv"] and r["done"], e, lab, r["TS"]))
    ft = os.path.join(d, "orb_target.dat")
    if os.path.exists(ft): print("".join("  " + l for l in open(ft)), end="")
    if d == "S0":
        for n, ok, e, lab, ts in rows: print("  %-14s conv=%d  labels=%s" % (n, ok, lab))
        continue
    cnt = {}
    for n, ok, e, lab, ts in rows:
        k = "".join(sorted(set(lab))) if ok else "notconv"; cnt[k] = cnt.get(k, 0) + 1
    print("  label sets over configs:", cnt)
    e0 = [e for n, ok, e, lab, ts in rows if n == "c111_ud" and e]
    e0 = e0[0] if e0 else None
    for n, ok, e, lab, ts in rows:
        if n.startswith("c111") or n.startswith("c112") or n.startswith("c211"):
            print("  %-14s conv=%d  E/Fe-E(c111_ud)=%s  labels=%s  TS=%s" % (n, ok, "%.2f" % (e - e0) if e and e0 else "-", lab, ts))
print("== orbital classes (corner-site frame, order z2 xz yz x2-y2 xy)")
for k, c in enumerate(classes): print("  %s  %s" % ("ABCDEFGHIJ"[k], np.round(c, 2)))
