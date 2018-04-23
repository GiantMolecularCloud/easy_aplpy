#####################################################################
#                          APLPY PLOTTING                           #
#####################################################################
# These functions will produce plots of channel maps, moment maps,  #
# pV diagrams, ... in a quality that (hopefully) allows publishing. #
#####################################################################

__all__ = ['map']


###################################################################################################

# load helper functions
#######################

from ._helpers import *


###################################################################################################

# plot a map
############

def map(fitsfile, **kwargs):
    """
    easy_aplpy.plot.map
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Provides a simple interface to APLpy plots. A single filename is enough to get a basic plot.

    Mandatory unnamed arguments:
        fitsfile    Path and file name of the fits image to be plotted


    Optional arguments:
        out         Path and file name of the created plot. If not specified, the plot will be
                    saved as png where the input image is located.
                    OR dictionary specifying matplotlib.figure.save parameters.

        figsize     Fiugre size as tupel in inches. Defaults to A4 size (8.267 * 11.692 inches)

        cmap        Colormap to plot the image. Defaults to viridis. Every named matplotlib
                    colormap can be used or any matplotlib colormap object.

        vmin        Minimum value for colormap normalization. Autoscales if not provided.
        vmax        Maximum value for colormap normalization. Autoscales if not provided.

        stretch     Colormap strech. APLpy provides 'linear', 'log', 'sqrt', 'arcsin', 'power'.

        recenter    Center the image on this location and set image width/height. You either have
                    to specify radius or width and height. Must be a list array containing an
                    astropy.SkyCoord object plus one or two angular distances:
                    [SkyCoord(...), 1*u.arcmin] or [SkyCoord(...), 1*u.arcmin, 2*u.arcmin]

        contour     List of contour elements. Each contour element must be a list of 'image file',
                    list of contour levels and list of colors for each contour. If only one color
                    is given, it will be used for all contours. When specifying a fits cube, you
                    have to provide the slice to be used as second list element. Valid mpl.contour
                    kwargs to fromat the contour can be given as last list element.
        clabel      Draw labels on the contours? Needs contours to be present. Either set to True
                    to use the matplotlib defaults or give a dictionary with allowed mpl.axes.clabel
                    arguments. Tip: Use clabel={'fmt': '%i'} to get labels without trailing zeros.
        legend      Add a legend that lists the contour objects if legend=True or provide a list of
                    strings to label the contours.

        colorbar    Show a colorbar. Defaults to True. Turn off with colorbar=None. Specify
                    colorbar position and label by a list of ['position','label'].

        scalebar    Show a scalebar. Must be a list of length, label and position given as
                    astropy.units, string and string. E.g. ['1.0*u.arcsec', '10pc','bottom left']

        beam        Plot a beam ellipse as given in the fits header at this position in the plot.
                    Defaults to 'bottom right'. Turn off by beam=None.

        circles, markers, polygons, arrows, ellipses, lines, rectangles (overlays in general)
                    List of overlay elements. Each overlay element needs to provide the necessary
                    information, e.g. center, radius and kwargs for a circle.
                    See http://aplpy.readthedocs.io/en/stable/api/aplpy.FITSFigure.html to check
                    what type of overlay needs which information. Positions and lengths need to be
                    given as SkyCoord and astropy.units. The last element needs to be a kwarg dict
                    which can also be empty if you do not want to provide e.g. linewidth or style.

        execute_code    This option allows to pass arbitrary code that is executed just before
                        saving the figure and can be used to access aplpy functionality that is
                        not mapped by easy_aplpy.plot. The code must be given in a list of strings.
                        Note this code is executed in the local namespace of easy_aplpy.plot.map.
                        Be careful with this argument! It executes anything you pass to it.

    General style settings
        Settings that do not have to be changed for each plot but maybe once per script or once per
        project. Often used ones are tick_label_xformat, ticks_xspacing and the corresponding
        settings for y. These can be accessed through easy_aplpy.settings
        import easy_aplpy
        easy_aplpy.settings.some_setting = ...
        Check the auto-completion of easy_aplpy.settings or look at the file settings.py to see
        what can be changed.


    examples:

    easy_aplpy.plot.map('test_data/map.fits')

    easy_aplpy.plot.map('test_data/map.fits',
            out  = {'filename': 'abc.png', 'dpi': 300, 'transparent': False},
            cmap = 'magma',
            vmin = 100,
            vmax = 8000,
            stretch  = 'log',
            recenter = [SkyCoord('00h47m33.07s -25d17m20.0s'), 40.0*u.arcsec, 32.0*u.arcsec],
            contour  = [['test_data/map.fits', [2e3,4e3,6e3], 'black'],
                        ['test_data/map.fits', [1e3,3e3,5e3], 'white']],
            clabel   = {'fmt': '%i'},
            legend   = True,
            colorbar = ['top', 'awesomeness scale'],
            scalebar = [1.0*u.arcsec, 'a few parsec', 'bottom'],
            beam     = 'bottom left',
            circles  = [[SkyCoord('00h47m33.07s -25d17m20.0s'), 10.0*u.arcsec, {'linewidth': 1.0}],
                        [SkyCoord('00h47m33.07s -25d17m20.0s'), 20.0*u.arcsec, {'linewidth': 1.0}]]
            )
    """

    fig = _set_up_figure(fitsfile, kwargs)
    _show_map(fig, kwargs)
    _recenter_plot(fig, kwargs)
    _show_contours(fig, kwargs)
    _overplot_regions(fig, kwargs)
    _show_colorbar(fitsfile, fig, kwargs)
    _show_scalebar(fig, kwargs)
    _show_beam(fig, kwargs)
    _show_label(fig, kwargs)
    _show_overlays(fig, kwargs)
    _show_ticksNlabels(fig, kwargs)
    _show_legend(fig, kwargs)
    _execute_code(fig, kwargs)
    _save_figure(fitsfile, fig, kwargs)


###################################################################################################
