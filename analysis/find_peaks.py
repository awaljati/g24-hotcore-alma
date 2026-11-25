"""Continuum, rms and >3 sigma emission/absorption features of a spectrum.
Usage: python find_peaks.py ../data/spectrum_A1.txt --vlsr 112.2"""
import argparse
from scipy.signal import find_peaks
from common import load_spectrum, continuum_rms, C_KMS

p = argparse.ArgumentParser(); p.add_argument('spectrum'); p.add_argument('--vlsr', type=float, default=112.2)
a = p.parse_args()
f, s = load_spectrum(a.spectrum)
cont, rms = continuum_rms(s)
print(f'continuum {cont*1e3:.1f} mJy, rms {rms*1e3:.1f} mJy')
fac = 1 + a.vlsr / C_KMS
for sign, lab in [(1, 'EM'), (-1, 'AB')]:
    idx, _ = find_peaks(sign * (s - cont), height=3 * rms, distance=3)
    for i in idx:
        print(f'{lab} ch{i:3d}  f_obs {f[i]/1e3:.4f}  f_rest {f[i]*fac/1e3:.4f} GHz  S/N {(s[i]-cont)/rms:+.1f}')
