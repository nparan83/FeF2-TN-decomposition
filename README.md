# FeF2-TN-decomposition

Data and code for

> N. Paran, *Separating the errors of a first-principles Néel temperature: DFT+U exchange, effective-field theory and Monte Carlo for FeF₂*, submitted to Journal of Physics and Chemistry of Solids.

Every number in the tables and figures of the article is produced from the files in this repository by the scripts listed below. The DFT calculations were run with Quantum ESPRESSO 7.2 on the TRUBA cluster (TÜBİTAK ULAKBİM); the Monte Carlo and EFT codes need only Python 3 with NumPy, and the figure scripts also need Matplotlib.

## Layout

| Folder | Contents |
|---|---|
| `calc/` | One compressed archive per calculation directory, with every Quantum ESPRESSO input (`*.in`) and output (`*.out`) and the small text files written next to them (orbital targets, configuration lists, fit results). Wavefunctions, charge densities and scratch directories are not included. Extract all archives in the repository root with `for f in calc/*.tar.gz; do tar xzf "$f"; done`. |
| `scripts/` | DFT input generation, orbital control, fitting, quality control, Monte Carlo (stage-3 runs) and EFT codes, and the SLURM job scripts. |
| `data/` | Derived data: fitted exchange constants, fit logs, EFT/MFA table, stage-3 Monte Carlo data and summaries, linear-response U, convergence and quality-control tables. |
| `mcB/` | Final Monte Carlo analysis: job files, raw Binder-cumulant data, finite-size-scaling (FSS) fits, autocorrelation times, validation, synthetic test sets and the script that prints Tables 2–3 and Section 3.4–3.5 numbers. |
| `mcDC/` | Heisenberg runs with the anisotropy mapped as D → 0.875 D (Section 3.5). |
| `figures/` | Scripts for Figs. 1–4, Fig. S1 and the graphical abstract (output in `figures/output/`). |
| `MANIFEST.txt` | MD5 checksums of all files. |

Pseudopotentials are not included: `Fe.pbe-spn-kjpaw_psl.0.2.1.UPF` and `F.pbe-n-kjpaw_psl.0.1.UPF` from PSlibrary (`scripts/fef2_setup.sh` downloads them).

## Calculation directories (`calc/`)

| Archive | Article | Contents |
|---|---|---|
| `S1` | Sections 2.1, 3.1 | AFM unit cell relaxed at U = 4 eV from the five minority-orbital starting points; `S1/orb_relax.dat` lists energy, orbital and u of each solution. The lowest one gives the common geometry, `data/u_relaxed_v3.txt` (u = 0.3008929724). |
| `CV3` | Section 2.1 | Cut-off and k-mesh convergence and supercell consistency in the target orbital state; summary `data/conv_CV3.txt`. |
| `geomchk` | Section 3.4 | Relaxation of u at U_LR; with `geom_dj` outputs in `top_level` (`fef2_gdj_*.out`): change of FM−AFM energy (−0.19%). |
| `v3_U3.0` … `v3_U7.0`, `v3_U4.51` | Sections 2.2–2.3, 3.1–3.2, Tables 1, S1, Fig. 2 | Production energy mapping: 14 collinear configurations per U in 1×1×1, 1×1×2, 2×1×1, 2×1×2 and 2×2×1 cells (`configs.dat`), the c111_ud orbital scans (`*_m1` … `*_m5`), the per-Fe orbital probes and `orb_target.dat`. |
| `hp_U5.00`, `hp_U4.44`, `hp_U4.53` | Section 2.4 | Self-consistent DFPT (hp.x) iterations for U with ortho-atomic projectors; summary `data/hp_scf.dat` (U_LR = 4.51 eV). |
| `atom` | Section 2.4, 3.2, Table S2 | Same with non-orthogonalized atomic projectors (`atom/hp_scf.dat`, U = 2.91 eV) and the eight-configuration fit at that U (`atom/v7a_U2.91/jsub_B.txt`). |
| `v7_113` | Section 2.3, Table S2 | 1×1×3 cell: J₂ and J₅ separately, J₁ + 4J₄ + J₆ (`v7_113/jsub_113.txt`). |
| `v7_c099`, `v7_c101` (+ `_fft`) | Section 3.2, Table S2 | Eight-configuration fits with c scaled by 0.99 and 1.01 (`jsub_B.txt`); `_fft` = FFT-grid probes. |
| `v7_lat_SR`, `v7_lat_BK` | Section 2.1, Table S2 | FM−AFM with the lattice constants of Stout & Reed and of Baur & Khan (`jsub_FA.txt`). |
| `sub_test` | Section 2.3 | Configuration-subset (selection-bias) test of the U = 3 eV fit (`sub_test/sub_fit.log`). |
| `top_level` | — | Job logs and calculations run in the top directory. |

The superseded first production rounds (directories S0 and U3.0–U7.0, run before per-Fe orbital control in a metastable orbital state) are not part of the article and are not included.

## Where each result comes from

| Article item | Data file(s) | Produced by |
|---|---|---|
| Fig. 1 | — (geometry of Section 2.1) | `figures/fig1_structure.py` |
| Exchange constants, Table 1, Table S1, Fig. 2 | `data/fit_v4b.log` (U = 3–7 eV), `data/fit_uhp.log` (U_LR), collected in `data/fef2_J.dat` (J·S², set B = J₁–J₄ + J₇) | `scripts/fef2_fit.py` on the outputs in `calc/v3_U*` |
| Orbital control and quality checks (Section 2.2: overlaps, moments, −TS) | `data/qcB_summary.txt`, `data/qcB_energies.txt`, `data/runs_v3*.log` | `scripts/fef2_orbfix.py`, `scripts/fef2_pa.py`, `scripts/fef2_qcB.py` |
| Control calculations (Section 2.3, Table S2) | `jsub_*.txt` in `calc/v7_*`, `calc/atom`, `v3_U4.51/jsub_B_sub8.txt` | `scripts/fef2_genB.py`, `scripts/fef2_jsub.py`, `scripts/fef2_dftB.sh` |
| U_LR (Section 2.4) | `data/hp_scf.dat`, `atom/hp_scf.dat` (in `calc/atom`) | `scripts/fef2_hp.sh` |
| MFA and EFT T_N (Tables 2–3) | `data/fef2_eft.txt`, `mcB/fef2_eftB.txt` | `scripts/fef2_eft.py`, `mcB/fef2_mkjobsB.py` |
| Monte Carlo raw data | `data/fef2_mcfine.dat` (stage 3), `mcB/mcB_*.dat` (final runs), `mcDC/fef2_mcfine.dat` | `scripts/fef2_mc.py`, `scripts/fef2_mcfine.py`, `mcB/fef2_mcB.py` |
| Monte Carlo T_N by finite-size scaling (Tables 2–3, Fig. 3) | `mcB/fss_bc.txt`, `mcB/fss_syn_all.txt`, `mcB/fss_heis082.txt`, `mcB/fss_heisDC.txt`, `mcB/fss_d0.txt`, `mcB/fss_val.txt` | `mcB/fef2_fss.py` |
| Autocorrelation times (Section 2.6) | `mcB/tau.dat`, `mcB/tau_sum.txt` | `mcB/fef2_tau.py` |
| Tables 2–3, decomposition and neutron-parameter errors (Section 3.4), Section 3.5, bcc Ising validation | `mcB/tabB.txt` | `mcB/fef2_tabB.py` |
| Figs. 3, 4, S1 and the decomposition printout | `data/`, `mcB/fss_bc.txt` | `figures/figs2_4_S1.py` |
| Graphical abstract | `data/fef2_eft.txt`, `mcB/fss_bc.txt` | `figures/graphical_abstract.py` |

Older fit logs in `data/` (`fit_v3*.log`, `fit_v4.log`) are intermediate rounds before all configurations were in the target orbital state; the article uses `fit_v4b.log` and `fit_uhp.log` only. The neutron parameters (Sears et al., Phys. Rev. Lett. 136, 226701 (2026)) and the experimental T_N = 78.418 K (Belanger & Yoshizawa, Phys. Rev. B 35, 4823 (1987)) are the only numbers entered by hand; they appear in `scripts/fef2_mcrun.py`, `scripts/fef2_eft.py` and the figure scripts.

The out-of-sample test of Section 2.3 uses configurations computed for a separate study that is not yet published; their energies are available from the author on request.

## Reproducing

```bash
python3 figures/fig1_structure.py
python3 figures/figs2_4_S1.py          # prints MFA/EFT/MC T_N and the factorization 1.520 = 1.342 x 1.079 x 1.050
python3 figures/graphical_abstract.py
for f in calc/*.tar.gz; do tar xzf "$f"; done
cd mcB && python3 fef2_tabB.py          # reproduces mcB/tabB.txt
```

`fef2_tabB.py` needs the extracted `calc/` archives only for its sections 7–8 (control calculations); sections 1–6 run from `mcB/` alone.

## License

MIT (see `LICENSE`).
