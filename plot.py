#####################################################################
#                          APLPY PLOTTING                           #
#####################################################################
# These functions will produce plots of channel maps, moment maps,  #
# pV diagrams, ... in a quality that (hopefully) allows publishing. #
#####################################################################

__all__ = ['map','map_grid']


###################################################################################################

# common import
###############

import os
import aplpy
import numpy
from astropy.coordinates import SkyCoord as SkyCoord
from astropy import units as u
from astropy.coordinates import Angle as Angle
from astropy.io import fits
from matplotlib import rc as rc
rc('text',usetex=True)

from .plot_helpers import *
from .settings import *


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

        stretch     Colormap strech. APLpy provides ‘linear’, ‘log’, ‘sqrt’, ‘arcsinh’, ‘power’.

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
            vmin = 0,
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

    print("\x1b[0;34;40m[easy_aplpy]\x1b[0m plotting map "+fitsfile)

    # set up figure
    figsize  = kwargs.get('figsize', (8.267,11.692))     # A4 in inches
    fig = aplpy.FITSFigure(fitsfile, figsize=figsize)

    # plot map
    cmap     = kwargs.get('cmap', 'viridis')             # the recommended cmap
    stretch  = kwargs.get('stretch', 'linear')
    vmin     = kwargs.get('vmin')                        # no default, aplpy scales automatically
    vmax     = kwargs.get('vmax')                        # no default, aplpy scales automatically
    fig.show_colorscale(cmap=cmap, vmin=vmin, vmax=vmax, stretch=stretch)

    # recenter image
    recenter = kwargs.get('recenter')
    if ( not recenter is None ) and test_recenter_format(recenter):
        if (len(recenter) == 2):
            fig.recenter(recenter[0].ra.degree, recenter[0].dec.degree, radius=recenter[1].to(u.degree).value)
        elif (len(recenter) == 3):
            fig.recenter(recenter[0].ra.degree, recenter[0].dec.degree, width=recenter[1].to(u.degree).value, height=recenter[2].to(u.degree).value)


    # plot contours
    contours = kwargs.get('contours')
    clabel   = kwargs.get('clabel')
    legend   = kwargs.get('legend')
    if contours:

        contournum = 0                          # conting variable for # of contours
        for idx,cont in enumerate(contours):
            if len(cont) == 3:
                fig.show_contour(data=cont[0], levels=cont[1], colors=cont[2])
            elif len(cont) == 4:
                # two options when four arguments are given: slice argument (int) as second or kwargs (dict) as last element
                if type(cont[1]) is int:
                    fig.show_contour(data=cont[0], slices=[cont[1]], levels=cont[2], colors=cont[3])
                elif type(cont[3]) is dict:
                    fig.show_contour(data=cont[0], levels=cont[1], colors=cont[2], **cont[3])
                else:
                    raise TypeError("Contour: could not interpret contour list.")
            elif len(cont) == 5:
                fig.show_contour(data=cont[0], slices=cont[1], levels=cont[2], colors=cont[3], **cont[4])
            else:
                raise TypeError("Contour: wrong number or format of contour parameters in contour "+str(cont)+".")

            if 'legend' in kwargs:
                if legend is True:
                    fig._ax1.collections[contournum].set_label(cont[0].replace('_','$\_$'))
                elif ( isinstance(legend, (list,tuple)) ):
                    fig._ax1.collections[contournum].set_label(legend[idx])
                else:
                    raise TypeError("Legend: either True or list of names for each contour.")
                contournum += len(cont[1])     # count up plotted contours

            if 'clabel' in kwargs:
                if clabel == True:
                    fig._layers['contour_set_'+str(idx+1)].clabel()
                if isinstance(clabel,dict):
                    fig._layers['contour_set_'+str(idx+1)].clabel(**clabel)

    # overplot regions
    regions = kwargs.get('regions')
    if isinstance(regions,(list,tuple)):
        for region in regions:
            fig.show_regions(region)

    # show colorbar
    colorbar = kwargs.get('colorbar', ['right',fits.open(fitsfile)[0].header['bunit']])
    if not colorbar is None:
        fig.add_colorbar()
        fig.colorbar.show()
        fig.colorbar.set_location(colorbar[0])
        fig.colorbar.set_axis_label_text(colorbar[1])
        if ( stretch == 'log' ):
            log_ticks = [float('{:.2f}'.format(round(x,int(-1*np.log10(kwargs['vmin']))))) for x in np.logspace(np.log10(kwargs['vmin']),np.log10(kwargs['vmax']),num=10, endpoint=True)]
            fig.colorbar.set_ticks(log_ticks)
        fig.colorbar.set_font(size=colorbar_fontsize)
        fig.colorbar.set_axis_label_font(size=colorbar_fontsize)
        fig.colorbar.set_frame_color(frame_color)

    # show scale bar
    scalebar = kwargs.get('scalebar')
    if isinstance(scalebar,list) and ( len(scalebar) == 3 ):
        fig.add_scalebar(length=scalebar[0].to(u.degree).value, label=scalebar[1], corner=scalebar[2], frame=scalebar_frame)
        fig.scalebar.set_font(size=scalebar_fontsize)
        fig.scalebar.set_linestyle(scalebar_linestyle)
        fig.scalebar.set_linewidth(scalebar_linewidth)
        fig.scalebar.set_color(scalebar_color)

    # show beam
    beam = kwargs.get('beam', 'bottom right')
    if not beam is None:
        fig.add_beam()
        fig.beam.show()
        fig.beam.set_corner(beam)
        fig.beam.set_frame(beam_frame)
        fig.beam.set_color(beam_color)

    # plot set overlay
    label = kwargs.get('label')
    if label:
        if isinstance(label, str):
            lbl = [[0.5,0.9],label, {}]
        elif isinstance(label,(list,tuple)):
            if ( len (label) == 2 ):
                lbl = [label[0], label[1], {}]
        fig.add_label(label[0][0], label[0][1], label[1].replace('_','$\_$'), color='black', relative=True, size=velo_fontsize, **label[2])


    # plot overlays
    circles = kwargs.get('circles')
    if circles:
        if all(isinstance(x,(list,tuple)) for x in circles):
            for circle in circles:
                fig.show_circles(xw=circle[0].ra.degree, yw=circle[0].dec.degree, radius=circle[1].to(u.degree).value, **circle[2])
        else:
            raise TypeError("Overlays: Must be list of lists. I.e. a single overlay needs double brackets [[]].")

    markers = kwargs.get('markers')
    if markers:
        if all(isinstance(x,(list,tuple)) for x in markers):
            for marker in markers:
                fig.show_markers(xw=marker[0].ra.degree, yw=marker[0].dec.degree, **marker[1])
        else:
            raise TypeError("Overlays: Must be list of lists. I.e. a single overlay needs double brackets [[]].")

    polygons = kwargs.get('polygons')
    if polygons:
        if all(isinstance(x,(list,tuple)) for x in polygons):
            for polygon in polygons:
                fig.show_polygons(polygon_list=polygon[0] **polygon[1])
        else:
            raise TypeError("Overlays: Must be list of lists. I.e. a single overlay needs double brackets [[]].")

    arrows = kwargs.get('arrows')
    if arrows:
        raise NotImplementedError("Overplotting arrows is not supported yet.")

    ellipses = kwargs.get('ellipses')
    if ellipses:
        raise NotImplementedError("Overplotting ellipses is not supported yet.")

    lines = kwargs.get('lines')
    if lines:
        raise NotImplementedError("Overplotting lines is not supported yet.")

    rectangles = kwargs.get('rectangles')
    if rectangles:
        raise NotImplementedError("Overplotting rectangles is not supported yet.")

    # add ticks + labels
    fig.tick_labels.show()
    fig.tick_labels.set_xformat(tick_label_xformat)
    fig.tick_labels.set_yformat(tick_label_yformat)
    fig.tick_labels.set_font(size=tick_label_fontsize)
    fig.ticks.show()
    fig.ticks.set_xspacing(ticks_xspacing.to(u.degree).value)
    fig.ticks.set_yspacing(ticks_yspacing.to(u.degree).value)
    fig.ticks.set_minor_frequency(ticks_minor_frequency)
    fig.ticks.set_color(ticks_color)
    fig.frame.set_color(frame_color)
    fig.axis_labels.set_font(size=tick_label_fontsize)

    # add legend
    legend = kwargs.get('legend')
    if legend:
        fig._ax1.legend(loc=0, fontsize=colorbar_fontsize)
        if isinstance(legend,dict):
            fig._ax1.legend(fontsize=colorbar_fontsize, **legend)

    # execute additional code passed by the user
    execute_code = kwargs.get('execute_code')
    if execute_code:
        if isinstance(execute_code, (list,tuple)):
            for codes in execute_code:
                exec(codes) in locals()
        else:
            raise TypeError("Execute code: Code to execute must be given in a list of strings")

    # save figure to disk
    out = kwargs.get('out',os.path.splitext(fitsfile)[0]+'.png')
    if isinstance(out,str):
        fig.save(out, dpi=300, transparent=True, adjust_bbox=True)
    if isinstance(out,dict):
        fig.save(**out)
    print("\x1b[0;34;40m[easy_aplpy]\x1b[0m saved plot as "+out)


###################################################################################################
