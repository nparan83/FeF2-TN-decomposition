# fef2_pa.py - per-Fe orbital control for FeF2 configs that fail with one m per spin species.
# Every Fe becomes its own species (Fe1..FeN, F = N+1), so each Fe gets its own starting eigenvector m.
# Run inside a v3_U* dir (needs orb_target.dat, orb_morder.txt, NAME.in from fef2_gen.py). ASCII only. numpy.
#   mk NAME CODE [s] : write NAME_pCODE.in (CODE = one digit m per Fe, input order; s = screening 50/600 Ry)
#   chk NAME CODE [s]: per-Fe overlap with the target orbital; exit 0 if converged and all on target
#   val NAME         : run NAME (a passing config) in per-Fe form with probe m; energy must equal NAME.out (exit 1 if not)
#   search NAME      : probe (1 SCF iteration) -> m per Fe -> full run; if that misses, greedy per-Fe m search
#                      (4 Fe: full settings; 8 Fe: screening, then one full run);
#                      on success NAME.out := found run (original kept as NAME_v1.out). Prints "PA_OK"/"PA_FAIL".
import sys, os, re, shutil, subprocess
import numpy as np
RY = 13605.693122994
src = open("../fef2_fit.py").read(); exec(src[src.index("def parse"):src.index("rows = []")])
BUDGET = {4: int(os.environ.get("PA_BUDGET4", "14")), 8: int(os.environ.get("PA_BUDGET8", "6"))}

def order():
    b = open("orb_morder.txt").read().split()
    hits = open("orb_success.txt").read().split() if os.path.exists("orb_success.txt") else []
    return [int(m) for m in sorted(b, key=lambda m: (-hits.count(m), b.index(m)))]

def tag(code, screen): return "_%s%s" % ("q" if screen else "p", "".join(map(str, code)))

def mk(name, code, screen=False):
    new = name + tag(code, screen)
    t = open(name + ".in").read()
    head, pos = t.split("ATOMIC_POSITIONS crystal\n")
    pl = pos.split("K_POINTS")[0].rstrip("\n").split("\n"); tail = "K_POINTS" + pos.split("K_POINTS", 1)[1]
    fe = [l for l in pl if re.match(r"Fe[12]\s", l)]; nfe = len(fe); assert len(code) == nfe
    spin = [1 if l.startswith("Fe1") else -1 for l in fe]
    U = re.search(r"U Fe1-3d (\S+)", t) or re.search(r"U Fe2-3d (\S+)", t); U = U.group(1)
    sysl = [l for l in head.split("\n") if not re.match(r"\s*(starting_magnetization|starting_ns_eigenvalue)", l)]
    head = "\n".join(sysl).replace("ntyp = 3", "ntyp = %d" % (nfe + 1))
    add = []
    for i, s in enumerate(spin):
        add.append("  starting_magnetization(%d) = %s" % (i + 1, "0.8" if s > 0 else "-0.8"))
        for m in range(1, 6):
            add.append("  starting_ns_eigenvalue(%d,%d,%d) = %s" % (m, 2 if s > 0 else 1, i + 1, "1.0" if m == code[i] else "0.0"))
    add.append("  starting_magnetization(%d) = 0.0" % (nfe + 1))
    head = head.replace("  occupations = 'smearing'", "\n".join(add) + "\n  occupations = 'smearing'")
    psf = re.search(r"Fe1 55.845 (\S+)", head).group(1); psff = re.search(r"\nF 18.998 (\S+)", head).group(1)
    sp = head.split("ATOMIC_SPECIES\n"); rest = sp[1][sp[1].index("CELL_PARAMETERS"):]
    head = sp[0] + "ATOMIC_SPECIES\n" + "".join("Fe%d 55.845 %s\n" % (i + 1, psf) for i in range(nfe)) + "F 18.998 %s\n" % psff + rest
    k = 0; npl = []
    for l in pl:
        if re.match(r"Fe[12]\s", l): k += 1; l = ("Fe%d" % k).ljust(3) + l[3:]
        npl.append(l)
    tail = re.sub(r"HUBBARD ortho-atomic\n.*", "HUBBARD ortho-atomic\n" + "\n".join("U Fe%d-3d %s" % (i + 1, U) for i in range(nfe)) + "\n", tail, flags=re.S)
    t = head + "ATOMIC_POSITIONS crystal\n" + "\n".join(npl) + "\n" + tail
    t = t.replace("prefix = '%s'" % name, "prefix = '%s'" % new).replace("tmp_%s_U" % name, "tmp_%s_U" % new)
    if screen:
        w, r_ = os.environ.get("FEF2_SCREENCUT", "50,600").split(",")
        t = re.sub(r"ecutwfc = \S+", "ecutwfc = " + w, t); t = re.sub(r"ecutrho = \S+", "ecutrho = " + r_, t)
        t = re.sub(r"\n  nr[123]s? = \d+", "", t)
        t = re.sub(r"conv_thr = \S+", "conv_thr = 1.0d-6", t); t = t.replace("electron_maxstep = 300", "electron_maxstep = 120")
    open(new + ".in", "w").write(t); return new

def ovl(name, run):
    ref = {int(l.split()[1]): np.array(l.split()[2:7], float) for l in open("orb_target.dat") if l.startswith("sub")}
    sl = sublat(name + ".in")
    if not os.path.exists(run + ".out"): return None, []
    r = parse(run + ".out")
    if len(r["atoms"]) != len(sl): return r, [0.0] * len(sl)
    return r, [abs(float(np.dot(a["minvec"], ref[s]))) for a, s in zip(r["atoms"], sl)]

def good(r, ov): return 0 if (r is None or not (r["conv"] and r["done"])) else sum(o > 0.95 for o in ov)

def run(name, code, screen, log):
    n = mk(name, code, screen)
    subprocess.run(["bash", "-c", ". ../fef2_run.sh && run_one " + n], check=False)
    r, ov = ovl(name, n); g = good(r, ov)
    e = r["E"] / len(ov) * RY if (r is not None and r["E"] is not None) else None
    log.write("%-28s %s conv=%d ontarget=%d/%d ov=%s E/Fe=%s\n" % (n, "scr " if screen else "full", int(bool(r and r["conv"] and r["done"])),
              g, len(ov), " ".join("%.2f" % o for o in ov), "%.3f" % e if e else "-")); log.flush()
    return n, g, len(ov)

def probe(name):
    # One SCF iteration at full settings; the ns printed after iteration 1 (before "Modifying starting ns")
    # holds the eigenvectors that starting_ns_eigenvalue(m,...) indexes. For each Fe pick the minority
    # eigenvector closest to the target orbital -> m per Fe. Deterministic (startingwfc='atomic').
    nfe = len(sublat(name + ".in")); n = mk(name, [1] * nfe); p = name + "_probe"
    t = open(n + ".in").read().replace("electron_maxstep = 300", "electron_maxstep = 1")
    t = t.replace("prefix = '%s'" % n, "prefix = '%s'" % p).replace("tmp_%s_U" % n, "tmp_%s_U" % p)
    open(p + ".in", "w").write(t); os.remove(n + ".in")
    if not (os.path.exists(p + ".out") and "Modifying starting ns" in open(p + ".out").read()):
        subprocess.run(["bash", "-c", ". ../fef2_run.sh && $MPIRUN pw.x -nk $NK -in %s.in > %s.out 2>&1; rm -rf tmp_%s_U*" % (p, p, p)], check=False)
    o = open(p + ".out").read()
    i = o.find("iteration #  1 "); j = o.find("Modifying starting ns", i)
    if i < 0 or j < 0: return None, []
    ref = {int(l.split()[1]): np.array(l.split()[2:7], float) for l in open("orb_target.dat") if l.startswith("sub")}
    code, best = [], []
    for a, s in zip(re.split(r"-+ ATOM\s+\d+ -+", o[i:j])[1:], sublat(name + ".in")):
        tr = re.search(r"\(up, down, total\) =\s+(\S+)\s+(\S+)", a); up, dn = float(tr.group(1)), float(tr.group(2))
        blk = re.split(r"SPIN\s+\d", a)[1:3][1 if up > dn else 0]
        ev = blk.split("eigenvalues:")[1].split("eigenvectors")[0].split()
        vec = np.array(blk.split("(columns):")[1].split("occupation")[0].split(), float).reshape(len(ev), len(ev))
        ov = np.abs(ref[s] @ vec); code.append(int(np.argmax(ov)) + 1); best.append(float(ov.max()))
    return code, best

def search(name):
    nfe = len(sublat(name + ".in")); screen = False; od = order(); log = open("pa_%s.log" % name, "a")
    code, pov = probe(name)
    if code and len(code) == nfe:
        log.write("probe: code=%s max overlap per Fe=%s\n" % ("".join(map(str, code)), " ".join("%.2f" % x for x in pov)))
        n, g, N = run(name, code, False, log)
        if g == N:
            if os.path.exists(name + ".out") and not os.path.exists(name + "_v1.out"): shutil.move(name + ".out", name + "_v1.out")
            shutil.copy(n + ".out", name + ".out"); print("PA_OK %s code=%s (probe)" % (name, "".join(map(str, code))))
            log.write("PA_OK %s\n" % n); return 0
        log.write("probe code missed (%d/%d) -> greedy fallback\n" % (g, N))
    screen = nfe > 4
    code = [od[0]] * nfe; tried = {}; best = None
    for m in od:                       # start from the best existing one-m-per-species run (no rerun)
        nm = "%s_%s%d" % (name, "s" if screen else "m", m)
        if not os.path.exists(nm + ".out"): continue
        r, ov = ovl(name, nm); g = good(r, ov)
        key = (g, -(r["E"] if r is not None and r["E"] is not None and r["conv"] else 1e9))
        if best is None or key > best[0]: best = (key, m, nm, g)
    if best is not None:
        code = [best[1]] * nfe; n, g, N = best[2], best[3], nfe; tried[tuple(code)] = g; nrun = 0
        log.write("start from %s: %d/%d on target\n" % (n, g, N))
    else:
        n, g, N = run(name, code, screen, log); tried[tuple(code)] = g; nrun = 1
    while g < N and nrun < BUDGET[8 if screen else 4]:
        r, ov = ovl(name, n); bad = [i for i, o in enumerate(ov) if o <= 0.95] or list(range(N))
        improved = False
        for i in bad:
            for m in od:
                if m == code[i]: continue
                c2 = list(code); c2[i] = m
                if tuple(c2) in tried: continue
                n2, g2, _ = run(name, c2, screen, log); tried[tuple(c2)] = g2; nrun += 1
                if g2 > g: code, n, g, improved = c2, n2, g2, True; break
                if nrun >= BUDGET[8 if screen else 4]: break
            if improved or nrun >= BUDGET[8 if screen else 4]: break
        if not improved: break
    if g == N and screen:
        n, g, N = run(name, code, False, log)
    if g == N:
        if os.path.exists(name + ".out") and not os.path.exists(name + "_v1.out"): shutil.move(name + ".out", name + "_v1.out")
        shutil.copy(n + ".out", name + ".out"); print("PA_OK %s code=%s runs=%d" % (name, "".join(map(str, code)), nrun))
        log.write("PA_OK %s\n" % n); return 0
    print("PA_FAIL %s best=%s %d/%d on target, runs=%d" % (name, "".join(map(str, code)), g, N, nrun)); log.write("PA_FAIL\n"); return 1

def val(name):
    nfe = len(sublat(name + ".in")); log = open("pa_%s.log" % name, "a")
    code, pov = probe(name)
    n, g, N = run(name, code if code and len(code) == nfe else [order()[0]] * nfe, False, log)
    r0, ov0 = ovl(name, name); r1, ov1 = ovl(name, n)
    if r0 is None or r1 is None or r0["E"] is None or r1["E"] is None: print("PA_VAL %s eksik" % name); return 1
    d = (r1["E"] - r0["E"]) / nfe * RY
    print("PA_VAL %s dE(per-Fe species - original) = %.4f meV/Fe  on target: original %d/%d, per-Fe %d/%d" % (
        name, d, good(r0, ov0), nfe, good(r1, ov1), nfe))
    return 0 if abs(d) < 0.01 and good(r1, ov1) == nfe else 1

if __name__ == "__main__":
    c = sys.argv[1]
    if c == "mk": print(mk(sys.argv[2], [int(x) for x in sys.argv[3]], len(sys.argv) > 4))
    elif c == "chk":
        r, ov = ovl(sys.argv[2], sys.argv[2] + tag([int(x) for x in sys.argv[3]], len(sys.argv) > 4))
        print(" ".join("%.2f" % o for o in ov)); sys.exit(0 if good(r, ov) == len(ov) and ov else 1)
    elif c == "val": sys.exit(val(sys.argv[2]))
    elif c == "search": sys.exit(search(sys.argv[2]))
