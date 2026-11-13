# FeF2 DFT+U -> J -> MFA / EFT / MC: Figs. 2, 3, 4 and S1 of the paper, plus the decomposition printout.
# Run from anywhere:  python3 figures/figs2_4_S1.py   (reads data/ and mcB/ of this repository; writes figures/output/)
# Every plotted number comes from a file in data/ or mcB/; the only hard-coded values are the neutron parameters of
# Sears et al. (same values as scripts/fef2_mcrun.py and scripts/fef2_eft.py) and T_N = 78.418 K (Belanger & Yoshizawa,
# Phys. Rev. B 35, 4823 (1987)). Same code as the Colab cell used for the paper; only the file locations differ.
import os, sys
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))          # repository root (this file is in figures/)
OUT = os.path.join(REPO, 'figures', 'output'); os.makedirs(OUT, exist_ok=True)
P = lambda f: os.path.join(REPO, f if f.startswith('mcB/') else os.path.join('data', f))
import re, glob, numpy as np, matplotlib
matplotlib.use(os.environ.get('MPLBACKEND', 'Agg'))
import matplotlib.pyplot as plt
TN_EXP = 78.418                                                   # K, Belanger & Yoshizawa, PRB 35, 4823 (1987), Table IV
NEU = {'J1': -0.030, 'J2': 0.46, 'J4': 0.005, 'J7': 0.011}        # meV, S=2 (J7 = (J7a+J7b)/2)
hp = [l.split() for l in open(P('hp_scf.dat')) if l.strip()]
U_HP = float([w for w in hp if w[1] == 'q2'][-1][4].split('=')[1])   # U_out of the q-converged hp.x step
# ---- exchange constants (set B) ----
J = {}
for l in open(P('fef2_J.dat')):
    if l.startswith('#'): continue
    w = l.split()
    if w[1].startswith('B'):
        J[float(w[0])] = {k: (float(v) / 4, float(e) / 4) for k, v, e in re.findall(r'(J\d)=(-?[\d.]+)\(([\d.]+)\)', l)}
Us = np.array(sorted(J))
U_LR = float(Us[np.argmin(abs(Us - U_HP))])                       # production U closest to the hp.x value
# ---- MFA / EFT ----
eft = {}
for l in open(P('fef2_eft.txt')):
    w = l.split()
    if len(w) == 4 and not l.startswith(('#', 'set')): eft[w[0]] = (float(w[2]), float(w[3]))
# ---- MC: precision crossings if present (v6 summary: local crossing, jackknife error), else coarse ----
fine = {}
if os.path.exists(P('fef2_mcfine_sum.txt')):
    for l in open(P('fef2_mcfine_sum.txt')):
        m = re.match(r'(\S+)\s+(bc|heis)\s+(\d+)/(\d+)\s+(\S+) \+- (\S+) \((\d+)\)', l)
        if m and 'EDGE' not in l and m.group(5) != 'nan':     # skip pairs without a crossing or with the crossing at the window edge
            fine.setdefault((m.group(1), m.group(2)), {})[(int(m.group(3)), int(m.group(4)))] = (float(m.group(5)), float(m.group(6)))
def mc_fine(name, model):
    d = fine.get((name, model))
    if not d: return None
    pairs = sorted(d); t2, s2 = d[pairs[-1]]; t1 = d[pairs[0]][0]
    return t2, max(s2 if np.isfinite(s2) else 0, abs(t1 - t2))          # largest-L pair; error covers the pair drift
coarse = {}
for l in open(P('fef2_mcsum.txt')):
    w = l.split()
    if len(w) > 3 and w[1] in ('bc', 'heis'):
        c = [float(x) for x in re.findall(r'\d+/\d+:(\d+\.\d+)', l)]
        if c: coarse[(w[0], w[1])] = (np.mean(c), abs(c[-1] - c[0]) / 2 if len(c) > 1 else 0)
fss = {}                                                          # Plan B: finite-size-scaling fits (mcB/fss_bc.txt)
if os.path.exists(P('mcB/fss_bc.txt')):
    for l in open(P('mcB/fss_bc.txt')):
        w = l.split()
        if not l.startswith('#') and len(w) > 6: fss[w[0]] = (float(w[3]), float(w[6]))
name = lambda U: [k for k in eft if k.startswith('dft_U%s_' % (('%g' % U) if U != int(U) else '%.1f' % U))][0]
def mc(nm, model):
    if model == 'bc' and nm in fss: return (fss[nm][0], fss[nm][1], 'FSS')
    f = mc_fine(nm, model)
    return (f[0], f[1], 'fine') if f else (coarse[(nm, model)][0], coarse[(nm, model)][1], 'coarse')
plt.rcParams.update({'font.size': 10, 'axes.linewidth': 0.8, 'savefig.dpi': 300, 'figure.dpi': 110})
C = {'MFA': '#888888', 'EFT': '#1f77b4', 'MC': '#d62728', 'exp': 'k'}
# ---- Fig. 1: J_n(U) ----
fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.9))
ax[0].errorbar(Us, [J[u]['J2'][0] for u in Us], yerr=[J[u]['J2'][1] for u in Us], fmt='o-', color='k', ms=4, label='DFT+U')
ax[0].axhline(NEU['J2'], ls='--', color='C3', label='neutron')
ax[0].set_ylabel(r'$J_2$ (meV)')
for k, mk, ls in (('J1', 'o', '-'), ('J4', 's', '--'), ('J7', '^', ':'), ('J3', 'v', '-')):
    ax[1].errorbar(Us, [J[u][k][0] for u in Us], yerr=[J[u][k][1] for u in Us], fmt=mk, ls=ls, ms=4 if k != 'J4' else 5, mfc='white' if k == 'J4' else None, label=r'$J_%s$' % k[1])
    if k in NEU: ax[1].axhline(NEU[k], ls='--', lw=0.8, color=ax[1].lines[-1].get_color())
ax[1].axhline(0, color='0.7', lw=0.6); ax[1].set_ylabel(r'$J_n$ (meV)')
for a in ax:
    a.axvline(U_LR, color='0.5', ls=':', lw=1); a.set_xlabel(r'$U$ (eV)')
ax[0].text(U_LR + 0.05, ax[0].get_ylim()[1] * 0.95, r'$U_{\rm LR}$', va='top', fontsize=9, color='0.4')
ax[0].legend(frameon=False, fontsize=8); ax[1].legend(frameon=False, fontsize=8, ncol=2)
ax[0].set_title('(a)', loc='left', fontsize=10); ax[1].set_title('(b)', loc='left', fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'fig2_J_vs_U.png')); fig.savefig(os.path.join(OUT, 'fig2_J_vs_U.pdf'))
# ---- Fig. 2: T_N(U) by method ----
fig, ax = plt.subplots(figsize=(4.2, 3.2))
mfa = [eft[name(u)][0] for u in Us]; ef = [eft[name(u)][1] for u in Us]
mcs = [mc(name(u), 'bc') for u in Us]
ax.plot(Us, mfa, 's-', color=C['MFA'], ms=4, label='MFA')
ax.plot(Us, ef, 'o-', color=C['EFT'], ms=4, label='EFT')
ax.errorbar(Us, [m[0] for m in mcs], yerr=[m[1] for m in mcs], fmt='D-', color=C['MC'], ms=4, capsize=2, label='MC')
for u, m in zip(Us, mcs):
    if m[2] == 'coarse': ax.plot(u, m[0], 'D', mfc='white', mec=C['MC'], ms=4, zorder=5)
ax.axhline(TN_EXP, color='k', ls='--', lw=1, label=r'exp. $T_N$')
ax.axvline(U_LR, color='0.5', ls=':', lw=1)
ax.set_xlabel(r'$U$ (eV)'); ax.set_ylabel(r'$T_N$ (K)'); ax.legend(frameon=False, fontsize=8)
if any(m[2] == 'coarse' for m in mcs): ax.set_title('(open: coarse MC)', fontsize=8, loc='left')
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'fig3_TN_vs_U.png')); fig.savefig(os.path.join(OUT, 'fig3_TN_vs_U.pdf'))
# ---- Fig. 3: method hierarchy for neutron and U_LR sets ----
nl = name(U_LR); rows = []
for lab, nm in (('neutron $J$', 'neutron'), (r'DFT+$U_{\rm LR}$ $J$', nl)):
    m = mc(nm, 'bc'); rows.append((lab, eft[nm][0], eft[nm][1], m[0], m[1]))
fig, ax = plt.subplots(figsize=(4.2, 3.2)); x = np.arange(2); wd = 0.26
for i, (k, col) in enumerate((('MFA', C['MFA']), ('EFT', C['EFT']), ('MC', C['MC']))):
    v = [r[1 + i] for r in rows]; e = None if k != 'MC' else [r[4] for r in rows]
    b = ax.bar(x + (i - 1) * wd, v, wd, color=col, yerr=e, capsize=2, label=k)
    for xi, vi in zip(x + (i - 1) * wd, v): ax.text(xi, vi + 2, '%+.0f%%' % (100 * (vi / TN_EXP - 1)), ha='center', fontsize=7)
ax.axhline(TN_EXP, color='k', ls='--', lw=1, label='exp.')
ax.set_xticks(x); ax.set_xticklabels([r[0] for r in rows]); ax.set_ylabel(r'$T_N$ (K)'); ax.legend(frameon=False, fontsize=8, ncol=2)
ax.set_ylim(0, max(r[1] for r in rows) * 1.15)
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'fig4_hierarchy.png')); fig.savefig(os.path.join(OUT, 'fig4_hierarchy.pdf'))
# ---- Fig. S1: Binder cumulants (precision MC, replica mean +- s.e.) ----
if os.path.exists(P('fef2_mcfine.dat')):
    D = {}
    for l in open(P('fef2_mcfine.dat')):
        if l.startswith('#'): continue
        w = l.split(); D.setdefault((w[0], w[1], int(w[2])), {}).setdefault(round(float(w[4]), 4), []).append(float(w[7]))
    sets_ = sorted({k[0] for k in D}); fig, ax = plt.subplots(len(sets_), 2, figsize=(7.0, 2.6 * len(sets_)), squeeze=False)
    for i, s in enumerate(sets_):
        for j, mod in enumerate(('bc', 'heis')):
            for L in sorted({k[2] for k in D if k[0] == s and k[1] == mod}):
                d = D[(s, mod, L)]; T = np.array(sorted(d)); u = np.array([np.mean(d[t]) for t in T])
                e = np.array([np.std(d[t], ddof=1) / np.sqrt(len(d[t])) if len(d[t]) > 1 else 0 for t in T])
                ax[i, j].errorbar(T, u, yerr=e, fmt='o-', ms=3, capsize=1.5, label='L=%d' % L)
            ax[i, j].set_title('%s, %s' % (s, mod), fontsize=9, loc='left'); ax[i, j].set_xlabel('T (K)'); ax[i, j].set_ylabel(r'$U_4$')
            ax[i, j].legend(frameon=False, fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, 'figS1_binder.png')); fig.savefig(os.path.join(OUT, 'figS1_binder.pdf'))
# ---- table: decomposition of the DFT+U_LR + EFT error ----
mn, mne = mc('neutron', 'bc')[:2]; ml, mle = mc(nl, 'bc')[:2]; el = eft[nl][1]
fJ, fE, fM = ml / mn, el / ml, mn / TN_EXP; tot = el / TN_EXP
print('\nhp.x U (q-converged) = %.4f eV; production set U = %.2f eV' % (U_HP, U_LR))
print('%-22s %8s %8s %12s' % ('set', 'T_MFA', 'T_EFT', 'T_MC (BC)'))
for lab, nm in (('neutron', 'neutron'), ('DFT+U_LR', nl)):
    m = mc(nm, 'bc'); print('%-22s %8.2f %8.2f %7.2f+-%.2f (%s)' % (lab, eft[nm][0], eft[nm][1], m[0], m[1], m[2]))
print('\nDFT+U_LR+EFT vs experiment: x%.3f = J x%.3f (+-%.3f)  *  EFT/MC x%.3f  *  MC(neutron)/exp x%.3f' % (
      tot, fJ, fJ * np.hypot(mne / mn, mle / ml), fE, fM))
lt = np.log(tot); print('log shares: J %.0f%%, EFT %.0f%%, model %.0f%%' % (100 * np.log(fJ) / lt, 100 * np.log(fE) / lt, 100 * np.log(fM) / lt))
print('figures written to', OUT)
