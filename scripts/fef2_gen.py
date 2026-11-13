# fef2_gen.py - FeF2 (rutile) collinear spin configurations for DFT+U energy mapping
# Writes QE 7.2 scf inputs for selected supercells, one U value per call.
# v3: startingwfc='atomic' (no random component -> reproducible first-iteration ns for orbital control).
# Exchange convention (same as Sears et al. arXiv:2601.04303): E = E0 + sum_bonds J_n S_i.S_j, J>0 = AFM.
# Usage: python3 fef2_gen.py U u_F [ecutwfc ecutrho] [outdir]   (u_F = F x-parameter of 4f site)
# ASCII only. Needs numpy.
import sys, os, itertools
import numpy as np
A, C = 4.6974, 3.3082           # exp. lattice (290 K), FeF2 rutile P4_2/mnm
KBASE = tuple(int(x) for x in os.environ.get("FEF2_KBASE", "6,6,8").split(","))  # Gamma-centred mesh of the 6-atom cell; supercells use KBASE/n (commensurate)
CELLS = [(1,1,1), (1,1,2), (2,1,1), (2,2,1), (2,1,2)]
DSH = {1:3.308, 2:3.711, 3:4.697, 4:5.745, 5:5.971, 6:6.616, 7:6.643}
PSF = {"Fe":"Fe.pbe-spn-kjpaw_psl.0.2.1.UPF", "F":"F.pbe-n-kjpaw_psl.0.1.UPF"}

def shells():
    sh = {n: [] for n in DSH}
    for i in range(-2, 3):
        for j in range(-2, 3):
            for k in range(-3, 4):
                for s in (0, 1):
                    v = np.array([i, j, k], float) + (0.5 if s else 0.0)
                    d = np.sqrt((v[0]*A)**2 + (v[1]*A)**2 + (v[2]*C)**2)
                    for n, dd in DSH.items():
                        if abs(d - dd) < 0.02: sh[n].append((s, v))
    return sh
SH = shells()

def sites(N):
    return [(i, j, k, b) for i in range(N[0]) for j in range(N[1]) for k in range(N[2]) for b in (0, 1)]

def coeffs(N, sig):
    st = sites(N); idx = {s: t for t, s in enumerate(st)}; out = []
    for n in sorted(DSH):
        tot = 0.0
        for t, (i, j, k, b) in enumerate(st):
            for (s, v) in SH[n]:
                p = np.array([i, j, k], float) + (0.5 if b else 0.0) + v
                pb = 0 if abs(p[0] - np.floor(p[0])) < 1e-6 else 1
                q = np.floor(p - (0.5 if pb else 0.0) + 1e-9).astype(int)
                tot += sig[t] * sig[idx[(q[0] % N[0], q[1] % N[1], q[2] % N[2], pb)]]
        out.append(tot / 2.0 / len(st))          # bonds per Fe
    return out

def configs():
    seen = {}; res = []
    for N in CELLS:
        nfe = 2 * N[0] * N[1] * N[2]
        for bits in itertools.product((1, -1), repeat=nfe - 1):
            sig = (1,) + bits
            cf = tuple(round(x, 6) for x in coeffs(N, sig))
            if cf in seen: continue
            seen[cf] = 1
            name = "c%d%d%d_%s" % (N + ("".join("u" if x > 0 else "d" for x in sig),))
            res.append((name, N, sig, cf))
    return res

def fft_lines(N, ecw, ecr):
    # commensurate FFT grids: supercell grid = unit-cell grid x N (file written by step 0 from cv_e90_afm.out)
    for fn in ("fft_base.txt", "../fft_base.txt"):
        if os.path.exists(fn):
            w = open(fn).read().split()
            if float(w[0]) != float(ecw) or float(w[1]) != float(ecr): return []
            b = [int(x) for x in w[2:8]]
            d = [b[i % 3] * N[i % 3] for i in range(3)] + [b[3 + i] * N[i] for i in range(3)]
            return ["  nr1 = %d" % d[0], "  nr2 = %d" % d[1], "  nr3 = %d" % d[2],
                    "  nr1s = %d" % d[3], "  nr2s = %d" % d[4], "  nr3s = %d" % d[5]]
    return []

def omc_lines(occ):
    # occupation-matrix control: minority channel (spin 2 for Fe1=up, spin 1 for Fe2=down) starts with the
    # single minority electron in eigenvector 'occ' of the first-iteration ns (5 = crystal-field favoured).
    # Without it a flipped Fe can land in a metastable orbital (local test: +240 meV per cell).
    L = []
    for sp, isp in ((1, 2), (2, 1)):
        for m in range(1, 6):
            L.append("  starting_ns_eigenvalue(%d,%d,%d) = %s" % (m, isp, sp, "1.0" if m == occ else "0.0"))
    return L

def write_input(name, N, sig, U, u, ecw, ecr, fn, occ=5):
    if os.environ.get("FEF2_TESTCUT"): ecw, ecr = os.environ["FEF2_TESTCUT"].split(",")   # local test only
    pos = []
    for (i, j, k, b), sg in zip(sites(N), sig):
        f0 = np.array([i, j, k], float) + (0.5 if b else 0.0)
        pos.append(("Fe1" if sg > 0 else "Fe2", f0))
    for i in range(N[0]):
        for j in range(N[1]):
            for k in range(N[2]):
                for f in ((u, u, 0), (-u, -u, 0), (0.5+u, 0.5-u, 0.5), (0.5-u, 0.5+u, 0.5)):
                    pos.append(("F", np.array([i, j, k], float) + np.array(f)))
    kk = [KBASE[d] // N[d] for d in range(3)]
    L = []
    L += ["&CONTROL", "  calculation = 'scf'", "  prefix = '%s'" % name, "  outdir = './tmp_%s_U%s'" % (name, U),
          "  pseudo_dir = '../pseudo'", "  verbosity = 'high'", "/"]
    L += ["&SYSTEM", "  ibrav = 0", "  nat = %d" % len(pos), "  ntyp = 3",
          "  ecutwfc = %s" % ecw, "  ecutrho = %s" % ecr, "  nspin = 2"] + fft_lines(N, ecw, ecr) + [
          "  starting_magnetization(1) = 0.8", "  starting_magnetization(2) = -0.8",
          "  starting_magnetization(3) = 0.0",
          "  occupations = 'smearing'", "  smearing = 'mv'", "  degauss = 0.005"] + omc_lines(occ) + ["/"]
    L += ["&ELECTRONS", "  startingwfc = 'atomic'", "  conv_thr = %s" % os.environ.get("FEF2_CONV", "1.0d-10"), "  mixing_beta = 0.3", "  electron_maxstep = 300", "/"]
    L += ["ATOMIC_SPECIES", "Fe1 55.845 %s" % PSF["Fe"], "Fe2 55.845 %s" % PSF["Fe"], "F 18.998 %s" % PSF["F"]]
    L += ["CELL_PARAMETERS angstrom", "%14.8f %14.8f %14.8f" % (A*N[0], 0, 0),
          "%14.8f %14.8f %14.8f" % (0, A*N[1], 0), "%14.8f %14.8f %14.8f" % (0, 0, C*N[2])]
    L += ["ATOMIC_POSITIONS crystal"]
    for sp, f in pos:
        f = (f / np.array(N, float)) % 1.0
        L.append("%-3s %14.10f %14.10f %14.10f" % (sp, f[0], f[1], f[2]))
    L += ["K_POINTS automatic", "%d %d %d 0 0 0" % tuple(kk)]
    L += ["HUBBARD ortho-atomic", "U Fe1-3d %s" % U, "U Fe2-3d %s" % U]
    open(fn, "w").write("\n".join(L) + "\n")

def special(U, u, od):
    # relax (AFM unit cell), convergence set, commensurate-k check (AFM/FM in 2x1x2)
    os.makedirs(od, exist_ok=True)
    afm = (1, -1); fm = (1, 1)
    write_input("relax_afm", (1,1,1), afm, U, u, "90", "1080", os.path.join(od, "relax_afm.in"))
    t = open(os.path.join(od, "relax_afm.in")).read()
    t = t.replace("&ELECTRONS", "&ELECTRONS\n  diagonalization = 'cg'")
    t = t.replace("calculation = 'scf'", "calculation = 'relax'\n  forc_conv_thr = 1.0d-4\n  etot_conv_thr = 1.0d-6")
    t = t.replace("&ELECTRONS", "&ELECTRONS").replace("  electron_maxstep = 300\n/", "  electron_maxstep = 300\n/\n&IONS\n  ion_dynamics = 'bfgs'\n/")
    open(os.path.join(od, "relax_afm.in"), "w").write(t)
    for ecw, ecr in (("70", "840"), ("90", "1080"), ("110", "1320")):
        for nm, sg in (("afm", afm), ("fm", fm)):
            write_input("cv_e%s_%s" % (ecw, nm), (1,1,1), sg, U, u, ecw, ecr, os.path.join(od, "cv_e%s_%s.in" % (ecw, nm)))
    for nm, sg in (("afm", afm), ("fm", fm)):
        fn = os.path.join(od, "cv_k8812_%s.in" % nm)
        write_input("cv_k8812_%s" % nm, (1,1,1), sg, U, u, "90", "1080", fn)
        t = open(fn).read(); t = t[:t.index("K_POINTS")] + "K_POINTS automatic\n8 8 12 0 0 0\n" + t[t.index("HUBBARD"):]; open(fn, "w").write(t)
    write_input("cv_orb4_afm", (1,1,1), afm, U, u, "90", "1080", os.path.join(od, "cv_orb4_afm.in"), occ=4)
    write_input("chk212_afm", (2,1,2), (1,-1)*4, U, u, "90", "1080", os.path.join(od, "chk212_afm.in"))
    write_input("chk212_fm", (2,1,2), (1,1)*4, U, u, "90", "1080", os.path.join(od, "chk212_fm.in"))
    print("special inputs written:", od)

if __name__ == "__main__":
    if sys.argv[1] == "special":
        special(sys.argv[2], float(sys.argv[3]), sys.argv[4]); sys.exit(0)
    U = sys.argv[1]; u = float(sys.argv[2])
    ecw = sys.argv[3] if len(sys.argv) > 3 else "90"; ecr = sys.argv[4] if len(sys.argv) > 4 else "1080"
    od = sys.argv[5] if len(sys.argv) > 5 else "U%s" % U
    os.makedirs(od, exist_ok=True)
    cf = configs()
    with open(os.path.join(od, "configs.dat"), "w") as g:
        g.write("# name nFe c1 c2 c3 c4 c5 c6 c7  (bonds per Fe, sum sigma_i sigma_j)\n")
        for name, N, sig, c in cf:
            write_input(name, N, sig, U, u, ecw, ecr, os.path.join(od, name + ".in"))
            g.write("%-28s %3d " % (name, len(sig)) + " ".join("%8.4f" % x for x in c) + "\n")
    print("configs:", len(cf), "dir:", od)
