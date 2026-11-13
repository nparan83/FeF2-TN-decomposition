# fef2_fit.py - parse FeF2 DFT+U scf outputs, quality checks, least-squares exchange fit
# Convention (Sears et al. arXiv:2601.04303): E = E0 + sum_bonds J_n S_i.S_j, J>0 AFM.
# Collinear mapping gives (J_n * S^2) in meV; J_n for S = 2 is also printed (divide by 4).
# Usage: python3 fef2_fit.py U4.0 [U5.0 ...]   (each dir: configs.dat + <name>.out)
# Orbital QC: signed minority eigenvector vs orb_target.dat (per sublattice), |overlap| >= 0.95.
# Writes fef2_J.dat (one line per U and fit set) and fef2_qc.dat. ASCII only. Needs numpy.
import sys, os, re
import numpy as np
RY = 13605.693122994   # meV
SETS = {"B_J1-J4+J7": [0, 1, 2, 3, 6], "A_J1-J4": [0, 1, 2, 3]}
# NOTE: in all supercells used c5 == c2 (J5 aliases onto J2) and c6 is constant (J6 -> E0):
#       the fitted "J2" is J2 + J5. Set B (with J7) is primary; omitting J7 biases J1, J3 (synthetic test).

def parse(fn):
    t = open(fn).read()
    r = {"done": "JOB DONE" in t, "conv": "convergence has been achieved" in t}
    e = re.findall(r"^!\s+total energy\s+=\s+(-?\d+\.\d+)\s+Ry", t, re.M)
    r["E"] = float(e[-1]) if e else None
    ts = re.findall(r"smearing contrib\. \(-TS\)\s+=\s+(-?\d+\.\d+)\s+Ry", t)
    r["TS"] = float(ts[-1]) if ts else None
    # last HUBBARD OCCUPATIONS block
    i = t.rfind("HUBBARD OCCUPATIONS")
    atoms = []
    if i > 0:
        blk = t[i:]
        for a in re.split(r"-+ ATOM\s+\d+ -+", blk)[1:]:
            m = re.search(r"Atomic magnetic moment for atom\s+\d+\s+=\s+(-?\d+\.\d+)", a)
            tr = re.search(r"\(up, down, total\) =\s+(\S+)\s+(\S+)", a)
            sp = []
            for s in re.split(r"SPIN\s+\d", a)[1:3]:
                ev = np.array(s.split("eigenvalues:")[1].split("eigenvectors")[0].split(), float)
                vv = s.split("(columns):")[1].split("occupation")[0].split()
                vec = np.array(vv, float).reshape(len(ev), len(ev))
                sp.append((ev, vec))
            if m and tr and len(sp) == 2:
                up, dn = float(tr.group(1)), float(tr.group(2))
                ev, vec = sp[1] if up > dn else sp[0]          # minority spin channel
                k = int(np.argmax(ev))
                atoms.append({"mag": float(m.group(1)), "minocc": ev[k], "minvec": vec[:, k]})
    r["atoms"] = atoms
    return r

def sublat(fn_in):
    # sublattice of each Fe: 0 = cell corner, 1 = body centre (from crystal coords x supercell size)
    t = open(fn_in).read()
    cell = re.search(r"CELL_PARAMETERS angstrom\n(.*)\n(.*)\n(.*)\n", t)
    nx = float(cell.group(1).split()[0]) / 4.6974
    return [0 if abs(float(x) * nx - round(float(x) * nx)) < 1e-3 else 1
            for x in re.findall(r"^Fe[12]\s+(\S+)", t[t.index("ATOMIC_POSITIONS"):], re.M)]

rows = []; qc = []
for d in sys.argv[1:]:
    U = d.rstrip("/").split("U")[-1]
    cfg = [l.split() for l in open(os.path.join(d, "configs.dat")) if not l.startswith("#")]
    P = {}
    for w in cfg:
        fo = os.path.join(d, w[0] + ".out")
        if os.path.exists(fo): P[w[0]] = (parse(fo), sublat(os.path.join(d, w[0] + ".in")))
    ref = {}
    ft = os.path.join(d, "orb_target.dat")
    if os.path.exists(ft):
        for l in open(ft):
            if l.startswith("sub"): w = l.split(); ref[int(w[1])] = np.array(w[2:7], float)
    elif "c111_ud" in P:
        for a, s in zip(P["c111_ud"][0]["atoms"], P["c111_ud"][1]): ref[s] = a["minvec"]
    data = []
    for w in cfg:
        name, nfe, c = w[0], int(w[1]), np.array(w[2:], float)
        if name not in P:
            qc.append((U, name, "MISSING")); continue
        r, sl = P[name]; flags = []
        if not (r["done"] and r["conv"] and r["E"] is not None): flags.append("NOTCONV")
        if r["TS"] is not None and abs(r["TS"]) > 1e-5: flags.append("TS=%.1e" % r["TS"])
        if len(r["atoms"]) != nfe: flags.append("NOOCC")
        else:
            for a in r["atoms"]:
                if abs(abs(a["mag"]) - 3.7) > 0.3: flags.append("MOM%.2f" % a["mag"])
            ov = [abs(float(np.dot(a["minvec"], ref[s]))) for a, s in zip(r["atoms"], sl) if s in ref]
            if ov and min(ov) < 0.95: flags.append("ORB%.2f" % min(ov))
        qc.append((U, name, " ".join(flags) if flags else "OK"))
        if not flags or all(f.startswith("TS") for f in flags):
            data.append((name, nfe, c, r["E"] / nfe * RY))
    if len(data) < 6:
        print("U=%s: only %d usable configs" % (U, len(data))); continue
    E = np.array([x[3] for x in data]); Cm = np.array([x[2] for x in data])
    for sname, sel in SETS.items():
        M = np.hstack([np.ones((len(E), 1)), Cm[:, sel]])
        if np.linalg.matrix_rank(M) < M.shape[1]:
            print("U=%s %s: rank deficient" % (U, sname)); continue
        p, *_ = np.linalg.lstsq(M, E, rcond=None)
        res = E - M @ p; rms = np.sqrt(np.mean(res ** 2))
        loo = []
        for k in range(len(E)):
            mk = np.ones(len(E), bool); mk[k] = False
            if np.linalg.matrix_rank(M[mk]) == M.shape[1]:
                loo.append(np.linalg.lstsq(M[mk], E[mk], rcond=None)[0][1:])
        sd = np.std(np.array(loo), axis=0) if loo else np.zeros(len(sel))
        rows.append((U, sname, len(E), rms, np.max(np.abs(res)), sel, p[1:], sd))
        print("U=%s %s  N=%d  rms=%.4f max=%.4f meV/Fe" % (U, sname, len(E), rms, np.max(np.abs(res))))
        for j, s, v, e in zip(range(len(sel)), sel, p[1:], sd):
            print("   J%d*S^2 = %9.4f +- %.4f meV   J%d(S=2) = %8.4f meV" % (s + 1, v, e, s + 1, v / 4))
with open("fef2_J.dat", "w") as g:
    g.write("# U set Nconf rms_meV_perFe maxres  then J_n*S^2 (meV) and LOO std, convention E=E0+sum J S.S, J>0 AFM\n")
    for U, s, n, rms, mx, sel, p, sd in rows:
        g.write("%s %s %d %.5f %.5f " % (U, s, n, rms, mx) + " ".join("J%d=%.5f(%.5f)" % (k + 1, v, e) for k, v, e in zip(sel, p, sd)) + "\n")
with open("fef2_qc.dat", "w") as g:
    for q in qc: g.write("%s %-28s %s\n" % q)
bad = [q for q in qc if q[2] != "OK"]
print("QC: %d outputs, %d flagged (see fef2_qc.dat)" % (len(qc), len(bad)))
for q in bad[:15]: print("  ", *q)
