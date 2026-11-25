# Run inside CASA:  casa --nogui -c 01_extract_spectrum.py
# Imports the ALMA cube, writes the header, finds the brightest pixel of the
# mean-intensity map (continuum peak of A1) and extracts its spectrum.
import os, shutil

fits = 'G24.78+0.08_BB_2_SPW_2.image.pbcor.fits'   # ALMA product (not in repo)
img = 'g24_spw2.im'

importfits(fitsimage=fits, imagename=img, overwrite=True)

hdr = imhead(img, mode='list')
with open('header.txt', 'w') as f:
    for k, v in hdr.items():
        f.write(f'{k}: {v}\n')

# Mean over all channels -> brightest pixel
if os.path.exists('g24_mean.im'):
    shutil.rmtree('g24_mean.im')
immoments(imagename=img, moments=[-1], outfile='g24_mean.im')
st = imstat('g24_mean.im')
x, y = int(st['maxpos'][0]), int(st['maxpos'][1])
print('Peak pixel:', x, y, '|', st['maxposf'])

# 5x5-pixel box (~1 beam) spectrum at the peak
specflux(imagename=img, region=f'centerbox[[{x}pix,{y}pix],[5pix,5pix]]',
         unit='GHz', logfile='spectrum_A1.txt', overwrite=True)
