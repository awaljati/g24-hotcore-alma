# Run inside CASA after 01: plots the mean-intensity map around A1 and NW.
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

ia.open('g24_mean.im'); img = ia.getchunk(dropdeg=True); ia.close()
x0, y0, r, pix = 565, 489, 150, 0.03        # A1 peak pixel, half-size, arcsec/pixel
sub = img[x0-r:x0+r, y0-r:y0+r].T * 1e3
plt.figure(figsize=(6, 5))
plt.imshow(sub, origin='lower', cmap='inferno', extent=[r*pix, -r*pix, -r*pix, r*pix])
plt.colorbar(label='mJy/beam')
plt.plot(0, 0, 'c+', ms=12)
plt.xlabel('ΔRA (arcsec)'); plt.ylabel('ΔDec (arcsec)')
plt.title('G24.78+0.08 — mean intensity, spw 2')
plt.savefig('mean_map.png', dpi=200)
