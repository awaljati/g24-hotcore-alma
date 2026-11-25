"""Rotation diagram of CH3-13CN, optical depth of CH3CN 12_6, column densities.
Usage: python rotation_diagram.py fits_A1.json [--exclude "CH3-13CN 12_5"]"""
import argparse, json
import numpy as np
from scipy.optimize import brentq
from common import LINES

k, h, c = 1.380649e-16, 6.62607e-27, 2.99792458e10   # cgs
MU, J, Q300, R1213 = 3.92197, 12, 14683.6, 40.4

p = argparse.ArgumentParser(); p.add_argument('fits'); p.add_argument('--exclude', nargs='*', default=[])
a = p.parse_args()
L = json.load(open(a.fits))['lines']

def upper(n):
    nu, elo, gu, sp, K = LINES[n]
    A = 1.16395e-20 * nu**3 * ((J**2 - K**2) / J) * MU**2 / (2 * J + 1)
    Eu = 1.43877 * (elo + nu / 29979.2458)
    Nu = 8 * np.pi * k * (nu * 1e6) ** 2 * L[n]['W'] * 1e5 / (h * c**3 * A)
    return Eu, np.log(Nu / gu), L[n]['eW'] / L[n]['W']

names = [n for n in L if n.startswith('CH3-13CN') and n not in a.exclude]
print('lines used:', ', '.join(names))
X, Y, E = map(np.array, zip(*[upper(n) for n in names]))
pfit, cov = np.polyfit(X, Y, 1, w=1 / E, cov='unscaled')
Trot = -1 / pfit[0]; eT = np.sqrt(cov[0, 0]) / pfit[0] ** 2
N13 = np.exp(pfit[1]) * Q300 * (Trot / 300) ** 1.5
ratio = L['CH3CN 12_6']['W'] / L['CH3-13CN 12_6']['W']
tau = brentq(lambda t: (1 - np.exp(-t)) / (1 - np.exp(-t / R1213)) - ratio, 0.01, 200)
print(f'T_rot = {Trot:.0f} ± {eT:.0f} K   N(CH3-13CN) = {N13:.2e} cm^-2')
print(f'tau(CH3CN 12_6) = {tau:.1f}   N(CH3CN) ≈ {R1213*N13:.1e} cm^-2 (12C/13C = {R1213})')
