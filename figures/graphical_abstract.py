# FeF2 graphical abstract (13 x 5 cm, 300 dpi). Run: python3 figures/graphical_abstract.py (reads data/fef2_eft.txt
# and mcB/fss_bc.txt of this repository; writes figures/output/). Same code as the Colab cell used for the paper.
import os, sys
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))          # repository root (this file is in figures/)
OUT = os.path.join(REPO, 'figures', 'output'); os.makedirs(OUT, exist_ok=True)
P = lambda f: os.path.join(REPO, f if f.startswith('mcB/') else os.path.join('data', f))
import numpy as np, matplotlib
matplotlib.use(os.environ.get('MPLBACKEND', 'Agg'))
import matplotlib.pyplot as plt
TN_EXP = 78.418                                   # K, Belanger & Yoshizawa, PRB 35, 4823 (1987), Table IV
eft = {w[0]: (float(w[2]), float(w[3])) for w in (l.split() for l in open(P('fef2_eft.txt'))) if len(w) == 4 and w[0] != 'set'}
mc = {w[0]: float(w[3]) for w in (l.split() for l in open(P('mcB/fss_bc.txt')) if not l.startswith('#')) if len(w) > 6}
n, d = 'neutron', 'dft_U4.51_B14'
fJ, fE, fM = mc[d] / mc[n], eft[d][1] / mc[d], mc[n] / TN_EXP; tot = eft[d][1] / TN_EXP
sh = np.array([np.log(fJ), np.log(fE), np.log(fM)]) / np.log(tot)
plt.rcParams.update({'font.size': 8, 'axes.linewidth': 0.7})
fig = plt.figure(figsize=(13 / 2.54, 5 / 2.54), dpi=300)
ax = fig.add_axes([0.095, 0.17, 0.39, 0.78]); C = ('#888888', '#1f77b4', '#d62728'); x = np.arange(2); wd = 0.27
for i, (lab, vals) in enumerate((('MFA', [eft[n][0], eft[d][0]]), ('EFT', [eft[n][1], eft[d][1]]), ('MC', [mc[n], mc[d]]))):
    ax.bar(x + (i - 1) * wd, vals, wd, color=C[i], label=lab)
    for xi, v in zip(x + (i - 1) * wd, vals): ax.text(xi, v + 2, '%+.0f%%' % (100 * (v / TN_EXP - 1)), ha='center', fontsize=5.5)
ax.axhline(TN_EXP, color='k', ls='--', lw=0.8); ax.text(1.62, TN_EXP + 2, 'exp.', fontsize=6, ha='right')
ax.set_xticks(x); ax.set_xticklabels(['neutron $J$', r'DFT+$U_{\rm LR}$ $J$']); ax.set_ylabel(r'$T_N$ (K)'); ax.set_ylim(0, 172)
ax.legend(frameon=False, fontsize=6, ncol=3, loc='upper left')
bx = fig.add_axes([0.53, 0.30, 0.44, 0.22]); left = 0
for s, c, lab, f in zip(sh, ('#7b3294', C[1], '#4d9221'), ('exchange', 'EFT', 'model'), (fJ, fE, fM)):
    bx.barh(0, s, left=left, color=c, height=0.6); bx.text(left + s / 2, 0, '%s\n%.0f%%' % (lab, 100 * s), ha='center', va='center', color='w', fontsize=6)
    left += s
bx.set_xlim(0, 1); bx.axis('off')
fig.text(0.75, 0.83, r'$T_N^{\rm EFT}(J_{\rm DFT})\,/\,T_N^{\rm exp} = %.2f$' % tot, ha='center', fontsize=8)
fig.text(0.75, 0.66, r'$= %.2f_{\rm exch}\ \times\ %.3f_{\rm EFT}\ \times\ %.2f_{\rm model}$' % (fJ, fE, fM), ha='center', fontsize=8)
fig.text(0.75, 0.13, 'FeF$_2$, $S=2$ Blume-Capel, $D$ from neutrons; shares on a log scale', ha='center', fontsize=5.5, color='0.3')
out = OUT
fig.savefig(os.path.join(out, 'graphical_abstract.png'), dpi=300); fig.savefig(os.path.join(out, 'graphical_abstract.pdf'))
print('GA: %.4f = %.4f x %.4f x %.4f; shares %s -> %s' % (tot, fJ, fE, fM, np.round(100 * sh, 1), out))
