"""Paper figures 1-3 from the extracted spectra.
Usage: python fit_lines.py ... (A1 and NW first), then
       python make_figures.py fits_A1.json fits_NW.json"""
import sys, json
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common import load_spectrum, to_kelvin, continuum_rms, LINES, C_KMS

fA1, fNW = sys.argv[1], sys.argv[2]
V_A1, V_NW = 112.2, 109.3
fa, sa = load_spectrum('../data/spectrum_A1.txt'); ca, _ = continuum_rms(sa)
fn, sn = load_spectrum('../data/spectrum_NW.txt'); cn, _ = continuum_rms(sn)

# Figure 1: labelled A1 spectrum
labels = [(220594.4231, 'CH$_3$CN 12$_6$', 'k'), (220539.3235, 'CH$_3$CN 12$_7$', 'k'),
          (220475.8072, 'CH$_3$CN 12$_8$', 'k'), (220403.9, 'CH$_3$CN 12$_9$?', 'gray'),
          (220585.2, 'HNCO 10$_{1,9}$', 'C3'), (220561.391, 'CH$_2$CHCN', 'C2'),
          (220599.9797, 'CH$_3^{13}$CN 12$_3$', 'C0'), (220570.3728, 'CH$_3^{13}$CN 12$_4$', 'C0'),
          (220532.3273, 'CH$_3^{13}$CN 12$_5$', 'C0'), (220485.8587, 'CH$_3^{13}$CN 12$_6$', 'C0'),
          (220430.9858, 'CH$_3^{13}$CN 12$_7$', 'C0'), (220367.7305, 'CH$_3^{13}$CN 12$_8$ (blend)', 'C0'),
          (220398.6765, '$^{13}$CO 2–1', 'C4')]
f = fa / 1e3; T = to_kelvin(sa, fa); Tc = to_kelvin(ca, 220400.)
fig, ax = plt.subplots(figsize=(11, 4.6))
ax.step(f, T, where='mid', color='k', lw=0.8); ax.axhline(Tc, color='gray', ls='--', lw=0.7)
ymax = T.max()
for i, (nr, lab, col) in enumerate(labels):
    fo = nr * (1 - V_A1 / C_KMS) / 1e3; yy = ymax * 1.03 + (i % 3) * 10
    ax.plot([fo, fo], [ymax * 0.98, yy], color=col, lw=0.7)
    ax.text(fo, yy + 1, lab, rotation=90, ha='center', va='bottom', fontsize=7.5, color=col)
for a_, b_ in [(220.3590, 220.3725), (220.3150, 220.3235)]:
    ax.axvspan(a_, b_, color='C4', alpha=0.12)
ax.set_xlim(f.min() - 0.003, f.max() + 0.003); ax.set_ylim(-5, ymax * 1.03 + 75)
ax.set_xlabel('Observed frequency (LSRK) [GHz]'); ax.set_ylabel('T$_B$ [K]')
ax2 = ax.secondary_xaxis('top', functions=(lambda x: x * (1 + V_A1 / C_KMS), lambda x: x / (1 + V_A1 / C_KMS)))
ax2.set_xlabel('Rest frequency (v$_{LSR}$ = 112.2 km s$^{-1}$) [GHz]')
plt.tight_layout(); plt.savefig('fig1_spectrum_A1.png', dpi=200); plt.close()

# Figure 2: rotation diagrams of CH3-13CN
k, h, c = 1.380649e-16, 6.62607e-27, 2.99792458e10
MU, J = 3.92197, 12

def upper(L, n):
    nu, elo, gu, sp, K = LINES[n]
    A = 1.16395e-20 * nu**3 * ((J**2 - K**2) / J) * MU**2 / (2 * J + 1)
    Eu = 1.43877 * (elo + nu / 29979.2458)
    Nu = 8 * np.pi * k * (nu * 1e6) ** 2 * L[n]['W'] * 1e5 / (h * c**3 * A)
    return Eu, np.log(Nu / gu), L[n]['eW'] / L[n]['W']

fig, ax = plt.subplots(figsize=(5.5, 4.2))
for path, col, lab, excl in [(fA1, 'C0', 'A1', []), (fNW, 'C3', 'NW', ['CH3-13CN 12_5'])]:
    L = json.load(open(path))['lines']
    names = [n for n in L if n.startswith('CH3-13CN')]
    pts = {n: upper(L, n) for n in names}
    use = [n for n in names if n not in excl]
    X, Y, E = map(np.array, zip(*[pts[n] for n in use]))
    pf = np.polyfit(X, Y, 1, w=1 / E)
    ax.errorbar(X, Y, yerr=E, fmt='o', color=col, capsize=3, label=f'{lab}: T$_{{rot}}$ ≈ {-1/pf[0]:.0f} K')
    for n in excl:
        if n in pts:
            ax.errorbar(pts[n][0], pts[n][1], yerr=pts[n][2], fmt='o', mfc='none', color=col, capsize=3)
    xx = np.linspace(150, 450); ax.plot(xx, np.polyval(pf, xx), color=col, ls='--')
ax.set_xlabel('E$_u$/k [K]'); ax.set_ylabel('ln(N$_u$/g$_u$)  [cm$^{-2}$]'); ax.legend(fontsize=9)
ax.set_title('CH$_3^{13}$CN rotation diagrams', fontsize=10)
plt.tight_layout(); plt.savefig('fig2_rotation.png', dpi=200); plt.close()

# Figure 3: continuum-subtracted A1 and NW spectra in the rest frame
fig, axs = plt.subplots(2, 1, figsize=(11, 6), sharex=True)
for ax, (fr, s, cont, v, lab, col, y) in zip(axs, [
        (fa, sa, ca, V_A1, 'A1  (v$_{LSR}$ = 112.2 km s$^{-1}$, T$_{cont}$ ≈ 35 K)', 'k', 0.52),
        (fn, sn, cn, V_NW, 'NW  (v$_{LSR}$ = 109.3 km s$^{-1}$, T$_{cont}$ ≈ 11 K)', 'C3', 0.72)]):
    ax.step(fr / 1e3 * (1 + v / C_KMS), to_kelvin(s - cont, fr), where='mid', color=col, lw=0.8)
    ax.axhline(0, color='gray', lw=0.5); ax.set_ylabel('T$_B$ − continuum [K]')
    ax.text(0.40, y, lab, transform=ax.transAxes, fontsize=9.5, color=col)
short = {'CH3CN 12_6': 'CH$_3$CN 12$_6$', 'HNCO 10_1,9': 'HNCO', 'CH3-13CN 12_4': '$^{13}$C 12$_4$',
         'CH2CHCN 24_1,24': 'CH$_2$CHCN', 'CH3CN 12_7': 'CH$_3$CN 12$_7$', 'CH3-13CN 12_5': '$^{13}$C 12$_5$',
         'CH3-13CN 12_6': '$^{13}$C 12$_6$', 'CH3CN 12_8': 'CH$_3$CN 12$_8$', 'CH3-13CN 12_7': '$^{13}$C 12$_7$'}
for n, lab in short.items():
    nr = LINES[n][0] / 1e3
    for ax in axs:
        ax.axvline(nr, color='C0', lw=0.4, ls=':')
    axs[0].text(nr, 95, lab, rotation=90, fontsize=7.5, ha='center', va='bottom', color='C0',
                bbox=dict(fc='w', ec='none', pad=0.5))
axs[0].set_ylim(-45, 160); axs[1].set_ylim(-20, 115); axs[1].set_xlim(220.363, 220.612)
axs[1].set_xlabel('Rest frequency [GHz]   ($^{13}$C = CH$_3^{13}$CN)')
plt.tight_layout(); plt.savefig('fig3_A1_vs_NW.png', dpi=200)
print('Saved fig1_spectrum_A1.png, fig2_rotation.png, fig3_A1_vs_NW.png')
