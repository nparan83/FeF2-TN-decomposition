# fef2_orbfix.py v3 - orbital-state control for FeF2 configs. ASCII only. Needs numpy. Run inside a U dir or S1/.
#   plan              : c111_ud_m{1..5}.in (minority electron started in eigenvector m)
#   target            : lowest-energy converged c111_ud_m* -> orb_target.dat, orb_morder.txt; c111_ud.out := it
#   check NAME        : exit 0 if NAME.out converged and every Fe minority orbital matches orb_target.dat
#   mkin NAME M [s]   : NAME_mM.in (full settings) or NAME_sM.in (screening: 50/600 Ry, conv 1e-6)
#   accept NAME OUT   : keep original as NAME_v1.out (once), copy OUT.out -> NAME.out
#   order             : m try-order (most successful m in this dir first, then orb_morder.txt)
#   reorder N         : configs.dat := first N configs in priority order (full list kept in configs_all.dat)
#   relaxpick         : S1/: lowest-energy finished relax_afm_m* -> ../u_relaxed_v3.txt, orb_relax.dat
import sys, os, re, shutil
import numpy as np
RY = 13605.693122994
src = open("../fef2_fit.py").read(); exec(src[src.index("def parse"):src.index("rows = []")])
PRIO = ("c111_uu c111_ud c112_uuud c112_uudd c211_uuud c211_uudd c221_uuuuuuud c212_uudddduu c221_uudddduu "
        "c212_uuududuu c212_uuududdu c212_uuuddudd c221_uuuuuddd c212_uudduddu").split()
def base(n): return re.sub(r"_[sm]\d$", "", n)
def mkin(name, m, screen=False):
    new = "%s_%s%d" % (name, "s" if screen else "m", m)
    t = open(name + ".in").read()
    t = re.sub(r"starting_ns_eigenvalue\((\d),(\d),(\d)\) = \S+",
               lambda g: "starting_ns_eigenvalue(%s,%s,%s) = %s" % (g.group(1), g.group(2), g.group(3), "1.0" if int(g.group(1)) == m else "0.0"), t)
    t = t.replace("prefix = '%s'" % name, "prefix = '%s'" % new).replace("tmp_%s_U" % name, "tmp_%s_U" % new)
    if screen:
        w, r_ = os.environ.get("FEF2_SCREENCUT", "50,600").split(",")
        t = re.sub(r"ecutwfc = \S+", "ecutwfc = " + w, t); t = re.sub(r"ecutrho = \S+", "ecutrho = " + r_, t)
        t = re.sub(r"\n  nr[123]s? = \d+", "", t)
        t = re.sub(r"conv_thr = \S+", "conv_thr = 1.0d-6", t); t = t.replace("electron_maxstep = 300", "electron_maxstep = 120")
    open(new + ".in", "w").write(t); return new
def target():
    lines = []
    for m in range(1, 6):
        n = "c111_ud_m%d" % m
        if not os.path.exists(n + ".out"): continue
        r = parse(n + ".out")
        if r["conv"] and r["done"] and len(r["atoms"]) == 2: lines.append((r["E"], m, r))
    if not lines: print("### TARGET YOK"); sys.exit(1)
    lines.sort(); e0, m0, r0 = lines[0]; sl = sublat("c111_ud_m%d.in" % m0)
    ref = {s: a["minvec"] for a, s in zip(r0["atoms"], sl)}
    hit = lambda r: all(abs(np.dot(a["minvec"], ref[s])) > 0.95 for a, s in zip(r["atoms"], sl))
    same = [m for e, m, r in lines if hit(r)]
    order = same + [m for e, m, r in lines if m not in same] + [m for m in range(1, 6) if m not in [x[1] for x in lines]]
    with open("orb_target.dat", "w") as g:
        g.write("# target = lowest-energy c111_ud (m=%d); try order: %s\n" % (m0, " ".join(map(str, order))))
        for s in sorted(ref): g.write("sub %d " % s + " ".join("%.4f" % x for x in ref[s]) + "\n")
        for e, m, r in lines:
            g.write("# m=%d  E-Emin = %8.2f meV/Fe  mag %.3f  %s  vec %s\n" % (m, (e - e0) / 2 * RY, abs(r["atoms"][0]["mag"]),
                    "target" if hit(r) else "other", " ".join("%.2f" % x for x in r["atoms"][0]["minvec"])))
    open("orb_morder.txt", "w").write(" ".join(map(str, order)) + "\n")
    accept("c111_ud", "c111_ud_m%d" % m0)
    print(open("orb_target.dat").read())
def check(name):
    if not os.path.exists(name + ".out"): return False
    ref = {int(l.split()[1]): np.array(l.split()[2:7], float) for l in open("orb_target.dat") if l.startswith("sub")}
    r = parse(name + ".out"); sl = sublat(base(name) + ".in")
    if not (r["conv"] and r["done"]) or len(r["atoms"]) != len(sl): return False
    return all(abs(np.dot(a["minvec"], ref[s])) > 0.95 for a, s in zip(r["atoms"], sl))
def accept(name, out):
    if out == name: return
    if os.path.exists(name + ".out") and not os.path.exists(name + "_v1.out"): shutil.move(name + ".out", name + "_v1.out")
    shutil.copy(out + ".out", name + ".out")
    mm = re.search(r"_m(\d)$", out)
    if mm and name != "c111_ud": open("orb_success.txt", "a").write(mm.group(1) + "\n")
def relaxpick():
    rows = []
    for m in range(1, 6):
        fo = "relax_afm_m%d.out" % m
        if not os.path.exists(fo): continue
        t = open(fo).read()
        if "End of BFGS Geometry Optimization" not in t: continue
        r = parse(fo); blk = t[t.index("Begin final coordinates"):]
        u = float(re.search(r"^F\s+(\S+)", blk, re.M).group(1))
        rows.append((r["E"], m, u, r))
    if not rows: print("### RELAX YOK"); sys.exit(1)
    rows.sort(); e0, m0, u0, r0 = rows[0]
    with open("orb_relax.dat", "w") as g:
        g.write("# relaxed AFM unit cell, U=4: m  E-Emin(meV/Fe)  u(F)  |mag|  minority orbital (z2 xz yz x2-y2 xy)\n")
        for e, m, u, r in rows:
            g.write("m=%d %9.2f %.6f %.3f  %s\n" % (m, (e - e0) / 2 * RY, u, abs(r["atoms"][0]["mag"]), " ".join("%.2f" % x for x in r["atoms"][0]["minvec"])))
    open("../u_relaxed_v3.txt", "w").write("%.10f\n" % u0)
    print(open("orb_relax.dat").read()); print("u(F) v3 =", u0, "(m=%d)" % m0)
c = sys.argv[1]
if c == "plan":
    for m in range(1, 6): mkin("c111_ud", m)
elif c == "target": target()
elif c == "check": sys.exit(0 if check(sys.argv[2]) else 1)
elif c == "mkin": print(mkin(sys.argv[2], int(sys.argv[3]), len(sys.argv) > 4))
elif c == "accept": accept(sys.argv[2], sys.argv[3])
elif c == "relaxpick": relaxpick()
elif c == "order":
    b = open("orb_morder.txt").read().split()
    hits = open("orb_success.txt").read().split() if os.path.exists("orb_success.txt") else []
    print(" ".join(sorted(b, key=lambda m: (-hits.count(m), b.index(m)))))
elif c == "reorder":
    if not os.path.exists("configs_all.dat"): shutil.copy("configs.dat", "configs_all.dat")
    L = [l for l in open("configs_all.dat")]; hd = [l for l in L if l.startswith("#")]; rows = [l for l in L if not l.startswith("#")]
    rows.sort(key=lambda l: PRIO.index(l.split()[0]) if l.split()[0] in PRIO else 99)
    open("configs.dat", "w").write("".join(hd + rows[:int(sys.argv[2])]))
