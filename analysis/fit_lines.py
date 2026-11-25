"""Gaussian fits of identified lines -> JSON (v_LSR, FWHM, T_peak, W).
Usage: python fit_lines.py ../data/spectrum_A1.txt --vlsr 112.2 --out fits_A1.json"""
import argparse, json
import numpy as np
from scipy.optimize import curve_fit
from common import load_spectrum, to_kelvin, continuum_rms, LINES, GROUPS, C_KMS

p = argparse.ArgumentParser(); p.add_argument('spectrum')
p.add_argument('--vlsr', type=float, required=True); p.add_argument('--out', required=True)
a = p.parse_args()

f, s = load_spectrum(a.spectrum)
cont, rms = continuum_rms(s)
T = to_kelvin(s - cont, f); Trms = to_kelvin(rms, 220400.)
shift = 1 - a.vlsr / C_KMS                       # rest -> observed

def multi(x, *p):
    y = np.zeros_like(x)
    for i in range(0, len(p), 3):
        amp, x0, sig = p[i:i+3]; y += amp * np.exp(-0.5 * ((x - x0) / sig) ** 2)
    return y

res = {}
for names, (lo, hi) in GROUPS:
    m = (f >= lo * shift) & (f <= hi * shift)
    p0 = []
    for n in names:
        p0 += [max(T[m].max() * 0.7, 3 * Trms), LINES[n][0] * shift, 1.5]
    try:
        p, cov = curve_fit(multi, f[m], T[m], p0=p0, sigma=np.full(m.sum(), Trms),
                           absolute_sigma=True, maxfev=20000)
    except RuntimeError:
        print(f'{", ".join(names)}: fit did not converge (blended), skipped')
        continue
    e = np.sqrt(np.diag(cov))
    for i, n in enumerate(names):
        amp, x0, sig = p[3*i:3*i+3]; ea, ex, es = e[3*i:3*i+3]
        fw = 2.3548 * abs(sig) / x0 * C_KMS; efw = 2.3548 * es / x0 * C_KMS
        W = 1.0645 * amp * fw; eW = W * np.hypot(ea / amp, efw / fw)
        res[n] = dict(v=C_KMS * (1 - x0 / LINES[n][0]), ev=ex / x0 * C_KMS,
                      Tpk=amp, eT=ea, fwhm=fw, efw=efw, W=W, eW=eW)
        print(f'{n:17s} v={res[n]["v"]:6.1f}±{res[n]["ev"]:.1f}  Tpk={amp:5.1f}±{ea:.1f} K  '
              f'FWHM={fw:4.1f}±{efw:.1f}  W={W:6.1f}±{eW:.1f} K km/s')
json.dump(dict(cont_K=float(to_kelvin(cont, 220400.)), rms_K=float(Trms), lines=res),
          open(a.out, 'w'), indent=1)
