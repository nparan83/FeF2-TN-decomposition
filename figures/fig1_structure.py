# FeF2 rutile unit cell with the exchange paths of the energy model (Fig. 1; no data needed).
# Geometry = production values of the paper: a = 4.6974, c = 3.3082 A, u = 0.3008929724 (u_relaxed_v3.txt, Section 2.1).
# Bond lengths and coordination numbers in the legend are computed here from that geometry, not typed in.
import os, itertools, numpy as np, matplotlib
matplotlib.use(os.environ.get('MPLBACKEND', 'Agg'))
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(REPO, 'figures', 'output'); os.makedirs(OUT, exist_ok=True)
A, C, U = 4.6974, 3.3082, 0.3008929724
S = np.array([A, A, C])
corners = [np.array(p, float) for p in itertools.product((0, 1), repeat=3)]
centre = np.array([0.5, 0.5, 0.5])
F = [np.array(p) for p in ((U, U, 0), (1 - U, 1 - U, 0), (U, U, 1), (1 - U, 1 - U, 1), (0.5 + U, 0.5 - U, 0.5), (0.5 - U, 0.5 + U, 0.5))]
# ---- shells from the geometry (Fe lattice: corner + body centre), distance and coordination number ----
fe_all = [np.array(t, float) + s for t in itertools.product(range(-3, 4), repeat=3) for s in (np.zeros(3), centre)]
d = sorted(round(float(np.linalg.norm(p * S)), 3) for p in fe_all if np.linalg.norm(p) > 1e-9)
shell = sorted(set(d))[:7]; z = {r: d.count(r) for r in shell}
J = {1: shell[0], 2: shell[1], 3: shell[2], 4: shell[3], 7: shell[6]}   # J5 = shell[4], J6 = shell[5] not drawn
# one representative bond per J, placed on different edges/faces to keep labels apart (fractional coordinates)
B = {1: ((1, 1, 0), (1, 1, 1)), 2: ((0.5, 0.5, 0.5), (1, 0, 1)), 3: ((0, 0, 0), (1, 0, 0)),
     4: ((0, 0, 0), (0, 1, 1)), 7: ((0, 1, 1), (1, 0, 1))}
for n, (p, q) in B.items():                                       # check that every drawn bond has its shell length
    assert abs(np.linalg.norm((np.array(q) - np.array(p)) * S) - J[n]) < 2e-3, n
COL = {1: '#E69F00', 2: '#0072B2', 3: '#009E73', 4: '#D55E00', 7: '#CC79A7'}    # Okabe-Ito
plt.rcParams.update({'font.size': 8, 'font.family': 'DejaVu Sans'})
fig = plt.figure(figsize=(3.5, 4.0), dpi=150)
ax = fig.add_subplot(projection='3d'); ax.set_proj_type('ortho')
ax.set_box_aspect((A, A, C)); ax.view_init(elev=14, azim=-66)
edges = [(p * S, q * S) for p, q in itertools.combinations(corners, 2) if np.isclose(np.abs(p - q).sum(), 1)]
ax.add_collection3d(Line3DCollection(edges, colors='0.55', linewidths=0.6, linestyles=':'))
# FeF6 octahedron of the body-centre Fe (Fe-F < 2.2 A)
oct_ = [f * S for f in F if np.linalg.norm((f - centre) * S) < 2.2]
octe = [(p, q) for p, q in itertools.combinations(oct_, 2) if np.linalg.norm(p - q) < 3.1]
ax.add_collection3d(Line3DCollection(octe, colors='0.75', linewidths=0.45))
for n, (p, q) in B.items():
    p, q = np.array(p) * S, np.array(q) * S
    ax.plot(*zip(p, q), color=COL[n], lw=2.4, solid_capstyle='round', zorder=5,
            label=r'$J_{%d}$: %.2f $\AA$, $z$ = %d' % (n, J[n], z[J[n]]))
X = np.array([f * S for f in F]); ax.scatter(*X.T, s=18, c='#9ccc65', edgecolors='0.3', linewidths=0.4, depthshade=False, zorder=6)
for p, up in [(c, True) for c in corners] + [(centre, False)]:
    r = p * S; ax.scatter(*r, s=60, c='#4a6fa5' if up else '#c0504d', edgecolors='k', linewidths=0.5, depthshade=False, zorder=7)
    ax.quiver(*(r - [0, 0, 0.4 if up else -0.4]), 0, 0, 0.8 if up else -0.8, color='k', lw=0.9, arrow_length_ratio=0.3, zorder=8)
ax.set_xlim(0, A); ax.set_ylim(0, A); ax.set_zlim(0, C); ax.set_axis_off()
for xyz, lab in (((0.5 * A, -0.45, -0.35), r'$a$'), ((-0.45, 0.5 * A, -0.35), r'$b$'), ((-0.35, -0.35, 0.5 * C), r'$c$')):
    ax.text(*xyz, lab, color='0.35', fontsize=9, ha='center', va='center')      # cell-edge labels
from matplotlib.lines import Line2D
h, l = ax.get_legend_handles_labels()
ha = [Line2D([], [], ls='', marker='o', ms=7, mfc='#4a6fa5', mec='k', mew=0.5), Line2D([], [], ls='', marker='o', ms=7, mfc='#c0504d', mec='k', mew=0.5),
      Line2D([], [], ls='', marker='o', ms=4.5, mfc='#9ccc65', mec='0.3', mew=0.4)]
fig.legend(ha, [r'Fe (corner, $S^z>0$)', r'Fe (body centre, $S^z<0$)', 'F'], loc='lower center', bbox_to_anchor=(0.5, 0.22), ncol=3,
           frameon=False, fontsize=6.5, handletextpad=0.3, columnspacing=0.8)
fig.legend(h, l, loc='lower center', bbox_to_anchor=(0.5, 0.06), ncol=2, frameon=False, fontsize=7, handlelength=1.6, columnspacing=1.0)
fig.subplots_adjust(left=0, right=1, top=1.02, bottom=0.18)
for e in ('png', 'pdf'): fig.savefig(os.path.join(OUT, 'fig1_structure.' + e), dpi=600)
print('shells (A):', shell, ' z:', [z[r] for r in shell]); print('saved', os.path.join(OUT, 'fig1_structure.png'))
