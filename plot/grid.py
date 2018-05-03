#####################################################################
#                          APLPY PLOTTING                           #
#####################################################################
# These functions will produce plots of channel maps, moment maps,  #
# pV diagrams, ... in a quality that (hopefully) allows publishing. #
#####################################################################

__all__ = ['grid']


###################################################################################################

# load helper functions
#######################

import numpy as np
from astropy.io import fits
from ._helpers import *


###################################################################################################

# plot a grid of maps
#####################

def grid(fitsfile, shape, channels, **kwargs):
    """
    easy_aplpy.plot.grid
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    Provides a simple interface to APLpy plots. This function generates a single image output which
    is saved to disk by default. A single filename is enough to get a basic plot but many more
    options are available for customisation. More details are given below in the description of the
    many optional arguments.

    Returned objects:
    matplotlib.figure instance
                    The figure instance is returned for further manual customisation. The individual
                    panels can be accessed through fig.get_axes() which returns a list of axes
                    isinstances of the panels.

    Mandatory unnamed arguments:
        fitsfile    Path and file name of the fits image to be plotted. Can be a 2D image (map),
                    3D image (cube) or a 2D position velocity diagram.
        shape       The number of columns and rows as a list, e.g. [2,3]
        channels    A list specifying the channels of the input fits file that should be plotted.


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
                    For position velocity diagram a list of astropy.units objects must be given as
                    [offset center, velocity center, width, height].

        contour     List of list of contour elements. Each contour element must be a list of
                    'image file', list of contour levels and list of colors for each contour. If
                    only one color is given, it will be used for all contours. When specifying a
                    fits cube, you have to provide the slice to be used as second list element.
                    Valid mpl.contour kwargs to format the contour can be given as last list element.
                    You can draw an arbitrary list of contours in each panel. If no contours should
                    be plotted in a particular panel just give an empty list. Structure:
                    [[[panel1/contour1],[panel1/contour2]], [[panel2/contour1],[panel2/contour2]], ... ]
                    Also see example below.

        clabel      Draw labels on the contours? Needs contours to be present. Either set to True
                    to use the matplotlib defaults or give a dictionary with allowed mpl.axes.clabel
                    arguments. Tip: Use clabel={'fmt': '%i'} to get labels without trailing zeros.
        legend      Draw a legend of the contours. Not implemented yet.

        colorbar    Show a colorbar. Defaults to True. Turn off with colorbar=None. Specify
                    colorbar position and label by a list of ['position','label']. For the position
                    "last panel" and "right" are implemented which draw a horizontal or vertical
                    colorbar as the last panel or to the right of the last panel, respectively.

        scalebar    Show a scalebar. Must be a list of length, label and position given as
                    astropy.units, string and string. E.g. ['1.0*u.arcsec', '10pc','bottom left'].

        beam        Plot a beam ellipse as given in the fits header at this position in the plot.
                    Defaults to 'bottom right'. Turn off by beam=None.

        circles, markers, polygons, arrows, ellipses, lines, rectangles (overlays in general)
                    List of lists of overlay elements. Each overlay element needs to provide the
                    necessary information, e.g. center, radius and kwargs for a circle.
                    See http://aplpy.readthedocs.io/en/stable/api/aplpy.FITSFigure.html to check
                    what type of overlay needs which information. Positions and lengths need to be
                    given as SkyCoord and astropy.units. The last element needs to be a kwarg dict
                    which can also be empty if you do not want to provide e.g. linewidth or style.
                    You can draw an arbitrary list of overlays in each panel. If no overlays should
                    be plotted in a particular panel just give an empty list. Structure:
                    [[[panel1/olay1],[panel1/olay]], [[panel2/olay],[panel2/olay2]], ... ]
                    Also see example below.

    General style settings
        Settings that do not have to be changed for each plot but maybe once per script or once per
        project. Often used ones are tick_label_xformat, ticks_xspacing and the corresponding
        settings for y. These can be accessed through easy_aplpy.settings like so:
        import easy_aplpy
        easy_aplpy.settings.some_setting = ...
        Check the auto-completion of easy_aplpy.settings or look at the file settings.py to see
        what can be changed.


    examples:

    easy_aplpy.plot.grid('cube.fits', [2,3], [150,200,250,300,350,400])

    easy_aplpy.plot.grid('cube.fits', [2,3], [150,200,250,300,350,400],
        out      = {'filename': 'cube.channelmap.complex.right.png', 'dpi': 300, 'transparent': True},
        vmin     = 0,
        vmax     = 60,
        stretch  = 'linear',
        contours = [[['cube.fits', 150, [2.5,5,10,20,40], 'black']],
                    [['cube.fits', 200, [2.5,5,10,20,40], 'black']],
                    [['cube.fits', 250, [2.5,5,10,20,40], 'black']],
                    [['cube.fits', 300, [2.5,5,10,20,40], 'black']],
                    [['cube.fits', 350, [2.5,5,10,20,40], 'black']],
                    [['cube.fits', 400, [2.5,5,10,20,40], 'black']]],
        circles  = [[],
                    [],
                    [[SkyCoord('00h47m33.07s -25d17m20.0s'), 10.0*u.arcsec, {'linewidth': 1.0, 'edgecolor':'red'}],[SkyCoord('00h47m33.07s -25d17m20.0s'), 20.0*u.arcsec, {'linewidth': 1.0, 'edgecolor':'red'}]],
                    [[SkyCoord('00h47m33.07s -25d17m20.0s'), 10.0*u.arcsec, {'linewidth': 1.0, 'edgecolor':'red'}],[SkyCoord('00h47m33.07s -25d17m20.0s'), 20.0*u.arcsec, {'linewidth': 1.0, 'edgecolor':'red'}]],
                    [],
                    []],
        scalebar = [5.0*u.arcsec, 'some distance', 'bottom'],
        colorbar = ['right', 'some units']
        )
    """

    kwargs = _check_image_type(fitsfile, kwargs)
    main_fig = _set_up_grid(fitsfile, shape, kwargs)
    panels = _grid_panels(fitsfile, shape, channels, kwargs)

    for panel in panels:
        if ( panel['type'] == 'map' ):

            fig = _set_up_panel_figure(main_fig, panel, kwargs)
            _show_map(panel['file'], fig, kwargs)
            _recenter_plot(panel['file'], fig, kwargs)
            _show_beam(panel['file'], fig, kwargs)
            _format_grid_ticksNlabels(panel, fig, kwargs)
            _show_channel_label(panel, fig, kwargs)
            _show_contours(panel['file'], fig, kwargs, panel)
            _overplot_regions(panel['file'], fig, kwargs, panel)
            _show_overlays(panel['file'], fig, kwargs, panel)
            _show_scalebar(panel['file'], fig, kwargs, panel)
            _execute_code(panel['file'], fig, kwargs, panel)

    _show_grid_colorbar(fitsfile, main_fig, fig, panels, kwargs)
    #_show_grid_legend(panel['file'], fig, kwargs)
    _save_figure(fitsfile, main_fig, kwargs)

    return main_fig


###################################################################################################
