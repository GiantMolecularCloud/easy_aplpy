import sys
sys.path.append('/Users/krieger/modules')
import easy_aplpy

from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib.pyplot as plt

plt.ion()


####################################################################################################

easy_aplpy.plot.map('map.fits',
    out = {'filename': 'map.simple.png'}
    )

easy_aplpy.plot.map('map.fits',
    out  = {'filename': 'map.complex.png', 'dpi': 300, 'transparent': False},
    cmap = 'magma',
    vmin = 100,
    vmax = 8000,
    stretch  = 'log',
    recenter = [SkyCoord('00h47m33.07s -25d17m20.0s'), 40.0*u.arcsec, 32.0*u.arcsec],
    contour  = [['map.fits', [2e3,4e3,6e3], 'black'],
                ['map.fits', [1e3,3e3,5e3], 'white']],
    clabel   = {'fmt': '%i'},
    legend   = True,
    colorbar = ['top', 'awesomeness scale'],
    scalebar = [1.0*u.arcsec, 'a few parsec', 'bottom'],
    beam     = 'bottom left',
    circles  = [[SkyCoord('00h47m33.07s -25d17m20.0s'), 10.0*u.arcsec, {'linewidth': 1.0, 'edgecolor':'red'}],
                [SkyCoord('00h47m33.07s -25d17m20.0s'), 20.0*u.arcsec, {'linewidth': 1.0, 'edgecolor':'red'}]]
    )

####################################################################################################

easy_aplpy.plot.map('cube.fits',
    out = {'filename': 'cube.simple.chan_int.png'},
    channel = 250
    )

easy_aplpy.plot.map('cube.fits',
    out = {'filename': 'cube.simple.chan_unit.png'},
    channel = 250*u.km/u.s
    )

easy_aplpy.plot.map('cube.fits',
    out      = {'filename': 'cube.complex.png', 'dpi': 300, 'transparent': False},
    channel  = 250*u.km/u.s,
    cmap     = 'viridis',
    vmin     = 1,
    vmax     = 60,
    stretch  = 'log',
    recenter = [SkyCoord('00h47m33.07s -25d17m20.0s'), 40.0*u.arcsec, 32.0*u.arcsec],
    contour  = [['map.fits', [2e3,4e3,6e3], 'black'],
                ['map.fits', [1e3,3e3,5e3], 'white']],
    clabel   = {'fmt': '%i'},
    legend   = True,
    colorbar = ['top', 'brightness temperature [K]'],
    scalebar = [1.0*u.arcsec, 'a few parsec', 'bottom'],
    beam     = 'bottom left',
    circles  = [[SkyCoord('00h47m33.07s -25d17m20.0s'), 10.0*u.arcsec, {'linewidth': 1.0, 'edgecolor':'red'}],
                [SkyCoord('00h47m33.07s -25d17m20.0s'), 20.0*u.arcsec, {'linewidth': 1.0, 'edgecolor':'red'}]]
    )

####################################################################################################

easy_aplpy.plot.map('pv.fits',
    out = {'filename': 'pv.simple.png'}
    )

easy_aplpy.plot.map('pv.fits',
    out     = {'filename': 'pv.complex.png'},
    cmap    = 'jet',
    vmin    = 0.1,
    vmax    = 17,
    stretch = 'log',
    labels  = ['offset ["]','velocity [km\,s$^{-1}$]'],
    recenter = [0.*u.arcsec, 250.*u.km/u.s, 40.*u.arcsec, 500.*u.km/u.s],
    contour  = [['pv.fits', [5,10,15], 'black']],
    clabel   = {'fmt': '%i'},
    colorbar = ['right', 'brightness temperature [K]']
    )

####################################################################################################

easy_aplpy.plot.grid('cube.fits', [2,3], [150,200,250,300,350,400])

easy_aplpy.plot.grid('cube.fits', [2,3], [150,200,250,300,350,400],
    out      = {'filename': 'cube.channelmap.complex.png', 'dpi': 300, 'transparent': True},
    vmin     = 0,
    vmax     = 60,
    stretch  = 'linear',
    contours = [[['cube.fits', 150, [2.5,5,10,20,40], 'black']],
                [['cube.fits', 200, [2.5,5,10,20,40], 'black']],
                [['cube.fits', 250, [2.5,5,10,20,40], 'black']],
                [['cube.fits', 300, [2.5,5,10,20,40], 'black']],
                [['cube.fits', 350, [2.5,5,10,20,40], 'black']],
                [['cube.fits', 400, [2.5,5,10,20,40], 'black']]],
    )

####################################################################################################
