"""Shared constants: beam, unit conversion and the line list (JPL/CDMS)."""
import numpy as np

C_KMS = 299792.458
BEAM_MAJ, BEAM_MIN = 0.18828, 0.15437      # arcsec
BEAM_PIX, BOX_PIX = 36.59132, 25            # pixels per beam, pixels in 5x5 box

# name: (rest freq MHz, E_low cm^-1, g_up, species, K or None)
LINES = {
    'CH3CN 12_6':     (220594.4231, 219.1546, 100, 'CH3CN', 6),
    'CH3CN 12_7':     (220539.3235, 283.6085,  50, 'CH3CN', 7),
    'CH3CN 12_8':     (220475.8072, 357.9389,  50, 'CH3CN', 8),
    'HNCO 10_1,9':    (220585.2000,  63.1900,  63, 'HNCO', None),
    'CH2CHCN 24_1,24':(220561.3910,  86.4030, 147, 'CH2CHCN', None),
    'CH3-13CN 12_4':  (220570.3728, 119.9152,  50, 'CH3C13N', 4),
    'CH3-13CN 12_5':  (220532.3273, 164.5757,  50, 'CH3C13N', 5),
    'CH3-13CN 12_6':  (220485.8587, 219.1402, 100, 'CH3C13N', 6),
    'CH3-13CN 12_7':  (220430.9858, 283.5961,  50, 'CH3C13N', 7),
}
# Lines fitted together, with fitting windows in REST frequency (MHz)
GROUPS = [
    (['CH3CN 12_6', 'HNCO 10_1,9'],             (220577.3, 220599.2)),
    (['CH2CHCN 24_1,24', 'CH3-13CN 12_4'],      (220556.4, 220576.4)),
    (['CH3CN 12_7'],                            (220531.3, 220546.3)),
    (['CH3-13CN 12_5'],                         (220527.3, 220536.3)),
    (['CH3CN 12_8', 'CH3-13CN 12_6'],           (220470.4, 220492.4)),
    (['CH3-13CN 12_7'],                         (220426.4, 220435.4)),
]


def load_spectrum(path):
    """specflux log -> frequency (MHz), flux density (Jy in the 5x5 box)."""
    d = np.loadtxt(path, comments='#')
    return d[:, 2] * 1e3, d[:, 4]


def to_kelvin(s_jy, nu_mhz):
    """Box flux density -> mean brightness temperature (K)."""
    i_beam = s_jy * BEAM_PIX / BOX_PIX
    return 1.222e6 * i_beam / ((nu_mhz / 1e3) ** 2 * BEAM_MAJ * BEAM_MIN)


def continuum_rms(s, nsig=2.5, niter=15):
    """Iterative sigma clipping -> continuum level and rms (same units as s)."""
    m = np.ones_like(s, bool)
    for _ in range(niter):
        c = np.median(s[m]); r = 1.4826 * np.median(np.abs(s[m] - c))
        m = np.abs(s - c) < nsig * r
    return c, r
