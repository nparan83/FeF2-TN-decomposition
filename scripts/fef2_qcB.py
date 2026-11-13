# fef2_qcB.py - Supplementary tables from the production DFT outputs (referee M6, M7): run in the FeF2 root.
# For every U directory and every configuration used in the fit (fef2_qc.dat OK lines):
#   total energy (Ry), energy per Fe relative to the AFM unit cell (meV/Fe), residual of the production set-B fit
#   (fef2_J.dat), smearing contribution -TS (meV per cell), smallest per-Fe orbital overlap, range of |m|.
# Writes qcB_energies.txt (one row per run) and prints a per-U summary.  Usage: python3 fef2_qcB.py
# ASCII only. numpy only.
import os, re
import numpy as np
RY = 13605.693122994
src = open("fef2_fit.py").read(); exec(src[src.index("def parse"):src.index("rows = []")])
J = {}
for l in open("fef2_J.dat"):
    if l.startswith("#") or "B_J1" not in l: continue
    w = l.split(); J[w[0]] = {k: float(v) for k, v in re.findall(r"(J\d)=(-?[\d.]+)\(", l)}
ok = {}
for l in open("fef2_qc.dat"):
    w = l.split()
    if len(w) >= 3 and w[2] == "OK": ok.setdefault(w[0], []).append(w[1])
out = open("qcB_energies.txt", "w")
out.write("# U name nFe E_total(Ry) dE/Fe_vs_AFM(meV) resid_setB(meV/Fe) -TS(meV/cell) min_overlap |m|min |m|max\n")
print("# U    N  min_overlap  |m| range      max|-TS| (meV/cell)  TS(FM)-TS(AFM) (meV/Fe)  rms resid (meV/Fe)")
for U in sorted(ok, key=float):
    d = "v3_U%s" % U
    if not os.path.isdir(d): d = "v3_U%s" % U.rstrip("0").rstrip(".")
    cf = {l.split()[0]: (int(l.split()[1]), np.array(l.split()[2:9], float)) for l in open(os.path.join(d, "configs.dat")) if not l.startswith("#")}
    ref = {int(l.split()[1]): np.array(l.split()[2:7], float) for l in open(os.path.join(d, "orb_target.dat")) if l.startswith("sub")}
    rows = []
    for n in ok[U]:
        r = parse(os.path.join(d, n + ".out")); sl = sublat(os.path.join(d, n + ".in"))
        ov = [abs(float(np.dot(a["minvec"], ref[s]))) for a, s in zip(r["atoms"], sl)]
        mg = [abs(a["mag"]) for a in r["atoms"]]
        rows.append((n, cf[n][0], cf[n][1], r["E"], (r["TS"] or 0.0) * RY, min(ov), min(mg), max(mg)))
    epf = {n: E / nf * RY for n, nf, c, E, *_ in rows}; e0 = epf.get("c111_ud", np.mean(list(epf.values())))
    jj = J.get(U) or J.get(U.rstrip("0").rstrip("."))
    jv = np.array([jj["J1"], jj["J2"], jj["J3"], jj["J4"], jj["J7"]]); X = np.array([[c[0], c[1], c[2], c[3], c[6]] for n, nf, c, *_ in rows])
    y = np.array([epf[r[0]] for r in rows]); res = y - X @ jv; res -= res.mean()
    for (n, nf, c, E, ts, mo, m1, m2), rr in zip(rows, res):
        out.write("%-5s %-24s %2d %18.8f %10.4f %8.4f %10.4f %6.3f %5.2f %5.2f\n" % (U, n, nf, E, epf[n] - e0, rr, ts, mo, m1, m2))
    tsd = (dict((r[0], r[4] / r[1]) for r in rows).get("c111_uu", np.nan) - dict((r[0], r[4] / r[1]) for r in rows).get("c111_ud", np.nan))
    print("%-5s %2d  %.4f      %.3f-%.3f     %.3e            %+.2e                %.4f" % (U, len(rows), min(r[5] for r in rows), min(r[6] for r in rows),
          max(r[7] for r in rows), max(abs(r[4]) for r in rows), tsd, np.sqrt((res ** 2).mean())))
out.close(); print("-> qcB_energies.txt")
