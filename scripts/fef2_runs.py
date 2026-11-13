# fef2_runs.py - per-run table for one or more v3 directories: wall time, SCF iterations, convergence,
# retry, orbital match with orb_target.dat. Totals per directory (node-hours). ASCII only. Needs numpy.
# Usage: python3 fef2_runs.py v3_U5.0 [v3_U6.0 ...] [-v]     (-v: list every run, else summary + failures)
import sys, os, re, glob
import numpy as np
src = open("fef2_fit.py").read(); exec(src[src.index("def parse"):src.index("rows = []")])
def wall(t):
    m = re.findall(r"PWSCF\s*:.*CPU\s+(.*?)\s*WALL", t)
    if not m: return None
    s = 0.0
    for v, u in re.findall(r"([\d.]+)\s*([dhms])", m[-1]): s += float(v) * {"d": 86400, "h": 3600, "m": 60, "s": 1}[u]
    return s
verbose = "-v" in sys.argv
for d in [a for a in sys.argv[1:] if a != "-v"]:
    ref = {int(l.split()[1]): np.array(l.split()[2:7], float) for l in open(os.path.join(d, "orb_target.dat")) if l.startswith("sub")}
    rows = []
    for fo in sorted(glob.glob(os.path.join(d, "*.out"))):
        n = os.path.basename(fo)[:-4]
        b = re.sub(r"(_[sm]\d)?(_david_fail)?$", "", n)
        fi = os.path.join(d, b + ".in")
        if not os.path.exists(fi) or n == b and os.path.exists(os.path.join(d, b + "_v1.out")) and False: continue
        t = open(fo).read(); r = parse(fo); sl = sublat(fi)
        it = re.findall(r"convergence has been achieved in\s+(\d+) iterations", t)
        ov = [abs(float(np.dot(a["minvec"], ref[s]))) for a, s in zip(r["atoms"], sl)] if len(r["atoms"]) == len(sl) else [0.0]
        ok = min(ov) > 0.95
        worst = r["atoms"][int(np.argmin(ov))]["minvec"] if len(r["atoms"]) == len(sl) else None
        e = r["E"] / len(sl) * 13605.693122994 if r["E"] is not None else None
        kind = "fail" if "_david_fail" in n else ("scr" if re.search(r"_s\d$", n) else ("full" if re.search(r"_m\d$", n) else "copy"))
        rows.append((n, kind, len(sl), wall(t), int(it[-1]) if it else None, r["conv"] and r["done"], ok, min(ov), worst, e))
    run = [x for x in rows if x[1] != "copy"]
    tot = sum(x[3] or 0 for x in run)
    print("== %s : %d runs, %.1f node-h total" % (d, len(run), tot / 3600))
    for k in ("scr", "full", "fail"):
        for nfe in (2, 4, 8):
            g = [x for x in run if x[1] == k and x[2] == nfe]
            if not g: continue
            w = [x[3] for x in g if x[3]]; its = [x[4] for x in g if x[4]]
            print("  %-4s nFe=%d  n=%3d  wall/run %6.1f min (max %6.1f)  iters %s  conv %d/%d  target %d/%d" % (
                k, nfe, len(g), np.mean(w) / 60 if w else 0, max(w) / 60 if w else 0,
                "%.0f" % np.mean(its) if its else "-", sum(x[5] for x in g), len(g), sum(x[6] for x in g), len(g)))
    for x in run:
        if verbose or (x[1] in ("full", "fail") and not (x[5] and x[6])):
            e0 = [y[9] for y in rows if y[0] == "c111_ud" and y[9]]
            de = "%9.2f" % (x[9] - e0[0]) if e0 and x[9] else "        -"
            print("   %-24s %-4s nFe=%d wall=%7.1f min it=%-4s conv=%d target=%d minov=%.2f E/Fe-E(c111_ud)=%s worst=%s" % (
                x[0], x[1], x[2], (x[3] or 0) / 60, x[4], x[5], x[6], x[7], de,
                " ".join("%.2f" % v for v in x[8]) if x[8] is not None else "-"))
