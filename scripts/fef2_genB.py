# fef2_genB.py - Plan-B input generator = fef2_gen.py (same physics and settings) with environment overrides:
#   FEF2_A, FEF2_C  lattice constants (Angstrom) for the positions; neighbour shells are still identified with
#                   the production lattice (topology unchanged), so +-1 % in c does not lose a shell
#   FEF2_CELLS      supercells, default "111,112,211,221,212" (as fef2_gen.py); "113" for the J5 test
#   FEF2_PROJ       Hubbard projector: ortho-atomic (default) | atomic
#   FEF2_KBASE      Gamma-centred k-mesh of the 6-atom cell (read by fef2_gen.py; "6,6,9" for 1x1x3 cells)
#   FEF2_FFTBASE    file "ecutwfc ecutrho nr1 nr2 nr3 nr1s nr2s nr3s" of the unit cell (supercell grid = base x N);
#                   "none" = let pw.x choose; unset = fef2_gen.py default (fft_base.txt in . or ..)
# Usage: python3 fef2_genB.py U u_F ecutwfc ecutrho outdir [name1 name2 ...]   (names: write only these)
# Writes outdir/configs.dat (bond counts per Fe, c1..c7) and outdir/genB_info.txt. ASCII only. numpy.
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, "fef2_gen.py")).read().split('\nif __name__ == "__main__":')[0])   # SH from production A, C
A0, C0 = A, C
if os.environ.get("FEF2_A"): A = float(os.environ["FEF2_A"])
if os.environ.get("FEF2_C"): C = float(os.environ["FEF2_C"])
CELLS = [tuple(int(ch) for ch in x) for x in os.environ.get("FEF2_CELLS", "111,112,211,221,212").split(",")]
PROJ = os.environ.get("FEF2_PROJ", "ortho-atomic")
FB = os.environ.get("FEF2_FFTBASE")
if FB is not None:
    def fft_lines(N, ecw, ecr):
        if FB == "none" or not os.path.exists(FB): return []
        w = open(FB).read().split()
        if float(w[0]) != float(ecw) or float(w[1]) != float(ecr): return []
        b = [int(x) for x in w[2:8]]
        d = [b[i] * N[i] for i in range(3)] + [b[3 + i] * N[i] for i in range(3)]
        return ["  nr1 = %d" % d[0], "  nr2 = %d" % d[1], "  nr3 = %d" % d[2], "  nr1s = %d" % d[3], "  nr2s = %d" % d[4], "  nr3s = %d" % d[5]]
if __name__ == "__main__":
    U, u, ecw, ecr, od = sys.argv[1], float(sys.argv[2]), sys.argv[3], sys.argv[4], sys.argv[5]
    only = sys.argv[6:]; os.makedirs(od, exist_ok=True); n = 0
    with open(os.path.join(od, "configs.dat"), "w") as g:
        g.write("# name nFe c1 c2 c3 c4 c5 c6 c7  (bonds per Fe, sum sigma_i sigma_j)\n")
        for name, N, sig, c in configs():
            if only and name not in only: continue
            fn = os.path.join(od, name + ".in"); write_input(name, N, sig, U, u, ecw, ecr, fn)
            if PROJ != "ortho-atomic":
                t = open(fn).read().replace("HUBBARD ortho-atomic", "HUBBARD " + PROJ); open(fn, "w").write(t)
            g.write("%-28s %3d " % (name, len(sig)) + " ".join("%8.4f" % x for x in c) + "\n"); n += 1
    open(os.path.join(od, "genB_info.txt"), "w").write("A %.6f C %.6f (shells from %.4f %.4f)  cells %s  proj %s  kbase %s  fftbase %s  U %s  u %.10f  ecut %s %s\n" % (
        A, C, A0, C0, ",".join("%d%d%d" % x for x in CELLS), PROJ, ",".join(map(str, KBASE)), FB, U, u, ecw, ecr))
    print("genB: %d configs -> %s (%s)" % (n, od, open(os.path.join(od, "genB_info.txt")).read().strip()))
