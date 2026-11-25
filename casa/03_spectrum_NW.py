# Run inside CASA after 01: spectrum of the NW core (cube maximum pixel).
specflux(imagename='g24_spw2.im', region='centerbox[[603pix,511pix],[5pix,5pix]]',
         unit='GHz', logfile='spectrum_NW.txt', overwrite=True)
