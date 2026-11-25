# Run inside CASA after 01: continuum-subtracted moment-0 maps of individual lines.
import os, shutil
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

def clean(*names):
    for d in names:
        if os.path.exists(d):
            shutil.rmtree(d)

# 0. Cut a 141 x 141 pixel (4.2") sub-cube around A1 and NW (keeps memory low)
clean('g24_sub.im')
imsubimage(imagename='g24_spw2.im', outfile='g24_sub.im', box='515,430,655,570')

# 1. Continuum subtraction: zeroth-order fit to line-free channels
clean('g24_line.im', 'g24_cont.im')
imcontsub(imagename='g24_sub.im', linefile='g24_line.im', contfile='g24_cont.im',
          fitorder=0, chans='120~175,295~320,365~385')

# 2. Moment 0 per line (channel ranges cover both A1 and NW)
lines = {'HNCO': '18~32', 'CH3CN 12_7': '96~110',
         'CH2CHCN': '57~70', 'CH3-13CN 12_4': '43~54'}
for nm, ch in lines.items():
    out = 'mom0_' + nm.replace(' ', '_') + '.im'
    clean(out)
    immoments(imagename='g24_line.im', moments=[0], chans=ch, outfile=out)

# 3. Figure
def load(name):
    ia.open(name); a = ia.getchunk(dropdeg=True); ia.close(); return a
cont = load('g24_cont.im')
if cont.ndim == 3:
    cont = cont[:, :, 0]
x0, y0, pix = 70, 70, 0.03
ext = [x0*pix, -x0*pix, -y0*pix, y0*pix]
pos = {'A1': (565-515, 489-430), 'NW': (603-515, 511-430)}
lev = [f*cont.max() for f in (0.1, 0.3, 0.5, 0.7, 0.9)]

panels = [('Continuum [mJy/beam]', cont.T*1e3)]
for nm in lines:
    panels.append((nm + ' [Jy/beam km/s]', load('mom0_' + nm.replace(' ', '_') + '.im').T))

fig, axs = plt.subplots(2, 3, figsize=(13, 8.5))
for ax, (title, data) in zip(axs.flat, panels):
    im = ax.imshow(data, origin='lower', cmap='inferno', extent=ext)
    ax.contour(cont.T, levels=lev, colors='cyan', linewidths=0.6, extent=ext, origin='lower')
    for lab, (px, py) in pos.items():
        ax.text(-(px-x0)*pix, (py-y0)*pix, lab, color='w', fontsize=9, ha='center', va='center')
    ax.set_title(title, fontsize=10)
    ax.set_xlabel('ΔRA (")'); ax.set_ylabel('ΔDec (")')
    plt.colorbar(im, ax=ax, fraction=0.046)
axs.flat[-1].axis('off')
plt.tight_layout()
plt.savefig('moment0_maps.png', dpi=200)
