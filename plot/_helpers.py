#####################################################################
#                          APLPY PLOTTING                           #
#####################################################################
# These functions will produce plots of channel maps, moment maps,  #
# pV diagrams, ... in a quality that (hopefully) allows publishing. #
#####################################################################

__all__ = ['_check_image_type','_set_up_figure','_set_up_grid','_set_up_panel_figure','_grid_panels','_show_map','_recenter_plot','_test_recenter_format','_show_contours','_overplot_regions','_show_colorbar','_show_grid_colorbar','_show_scalebar','_show_beam','_show_label','_show_overlays','_format_grid_ticksNlabels','_show_ticksNlabels','_show_legend','_show_channel_label','_execute_code','_save_figure']


###################################################################################################

# common imports
################

import os
import aplpy
import numpy as np
import warnings
from astropy.coordinates import SkyCoord as SkyCoord
from astropy import units as u
from astropy.coordinates import Angle as Angle
from astropy.io import fits
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc as rc
rc('text',usetex=True)

import easy_aplpy


###################################################################################################
# plotting subfunctions
#
# TODO:
# check naxis to be 2 or 3
#
###################################################################################################

def _check_image_type(fitsfile, kwargs):
    header = fits.open(fitsfile)[0].header
    ctypes = [header['ctype'+str(i)] for i in np.arange(1,1+header['naxis'])]
    if ( [i for i in ctypes if i in ['OFFSET','offset','POSITION','position']] ):
        kwargs['imtype'] = 'pv'
    else:
        kwargs['imtype'] = 'pp'
    return kwargs


###################################################################################################

def _set_up_figure(fitsfile, kwargs):
    print("\x1b[0;34;40m[easy_aplpy]\x1b[0m plotting map "+fitsfile)
    figsize = kwargs.get('figsize', None)     # A4 in inches: (8.267,11.692)``
    channel = kwargs.get('channel', None)
    if ( channel == None ):
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            fig = aplpy.FITSFigure(fitsfile, figsize=figsize)
    else:
        channel,physical = _channel_physical(fitsfile, channel)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore")
            fig = aplpy.FITSFigure(fitsfile, slices=[channel], figsize=figsize)
    return fig


###################################################################################################

def _channel_physical(fitsfile, user_channel):
    header = fits.open(fitsfile)[0].header
    naxis = header['naxis']
    #TODO automatically convert velocity/frequency if necessary: http://docs.astropy.org/en/stable/units/equivalencies.html#spectral-doppler-equivalencies
    for i in np.arange(1,naxis+1):
        if header['cunit'+str(i)] in ['m/s','km/s','Hz','kHz','MHz','GHz']:
            freqax = str(i)
    if not freqax:
        raise TypeError("Could not find a velocity/frequency axis in the input image.")
    crval = u.Quantity(str(header['crval'+freqax])+header['cunit'+freqax])
    cdelt = u.Quantity(str(header['cdelt'+freqax])+header['cunit'+freqax])
    crpix = int(header['crpix'+freqax])

    if isinstance(user_channel, int):
        channel  = user_channel                                                # channel number already given by user
        physical = (user_channel-crpix)*cdelt+crval                            # physical value calculated
    elif isinstance(user_channel, u.quantity.Quantity):
        channel  = int(((user_channel-crval)/cdelt)+crpix)                     # channel number calculated
        physical = user_channel                                                # physical value already given by user
    else:
        raise TypeError("channel needs to be channel number (int) or velocity/frequency (astropy.units object.)")
    if ( 'Hz' in str(physical) ):
        physical = physical.to(u.GHz)
    if ( 'm / s' in str(physical) ):
        physical = physical.to(u.km/u.s)
    return channel,physical


###################################################################################################

def _set_up_grid(fitsfile, shape, kwargs):
    print("\x1b[0;34;40m[easy_aplpy]\x1b[0m plotting a "+str(shape[0])+"x"+str(shape[1])+" grid")
    figsize = kwargs.get('figsize', (8.267,11.692))     # A4 in inches
    main_fig = plt.figure(figsize=figsize)
    return main_fig


###################################################################################################

def _set_up_panel_figure(main_fig, panel, kwargs):
    if ( panel['type'] == 'map' ):
        print("\x1b[0;34;40m[easy_aplpy]\x1b[0m plotting panel "+str(panel['num']+1)+" of "+str(panel['npanels'])+", file: "+panel['file'])
    elif ( panel['type'] == 'colorbar' ):
        print("\x1b[0;34;40m[easy_aplpy]\x1b[0m plotting panel "+str(panel['num']+1)+" of "+str(panel['npanels'])+", colorbar")
    figsize = kwargs.get('figsize', (8.267,11.692))     # A4 in inches
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        fig = aplpy.FITSFigure(panel['file'], figure=main_fig, subplot=[panel['x'],panel['y'],panel['width'],panel['height']], dimensions=[0,1], slices=[panel['channel']])
    return fig


###################################################################################################

def _grid_panels(fitsfile, shape, channels, kwargs):
    ncols = float(shape[0])                                                    # convert to nrows/ncols to float for python 2 compatibility
    nrows = float(shape[1])                                                    # convert to nrows/ncols to float for python 2 compatibility
    head = fits.open(fitsfile)[0].header
    if ( 'bunit' in head ):
        bunit = head['bunit']
    else:
        print("Header keyword 'BUNIT' not present in "+fitsfile)
        bunit = 'unknown quantity'
    colorbar = kwargs.get('colorbar', ['right',bunit])        # add the colorbar panel

    # all panels to plot
    margins = easy_aplpy.settings.margins
    panels_width  = (1.-margins[0]-margins[1])
    panels_height = (1.-margins[2]-margins[3])
    panels = []
    for idx, channel in enumerate(channels):
        pos = ''
        if ( idx%ncols == 0 ):
            pos += 'left'
        if ( idx >= (nrows-1)*ncols ):
            pos += 'bottom'
        channel,physical = _channel_physical(fitsfile,channel)
        panels.append({'num': idx,
                       'npanels': shape[0]*shape[1],
                       'position': pos,
                       'type': 'map',
                       'x': margins[0]+(idx%ncols)*panels_width/ncols,                         # lower left corner
                       'y': (1.-margins[2])-np.ceil((idx+1.)/ncols)*panels_height/nrows,       # lower left corner
                       'width': panels_width/ncols,
                       'height': panels_height/nrows,
                       'file': fitsfile,
                       'channel': channel,
                       'physical': physical
                       })

    if not colorbar is None:
        idx = idx+1
        if ( colorbar[0] == 'last panel' ):
            panels.append({'num': idx+1,
                           'npanels': shape[0]*shape[1],
                           'position': None,
                           'type': 'colorbar',
                           'x': margins[0]+((idx%ncols)+0.05)*panels_width/ncols,                  # lower left corner
                           'y': (1.-margins[2])-(np.ceil((idx+1.)/ncols)-0.5)*panels_height/nrows, # lower left corner
                           'width': panels_width/ncols*0.9,                                        # cbar fills 90% of panel width
                           'height': panels_height/nrows*easy_aplpy.settings.colorbar_width,       # cbar height is specified fraction of panel height
                           'file': fitsfile,
                           'channel': None,
                           'physical': None
                           })
        if ( colorbar[0] == 'right' ):
            panels.append({'num': idx+1,
                           'npanels': shape[0]*shape[1],
                           'position': None,
                           'type': 'colorbar',
                           'x': margins[0]+((idx%ncols)+ncols)*panels_width/ncols,                 # lower left corner
                           'y': (1.-margins[2])-(np.ceil((idx+1.)/ncols)-1.)*panels_height/nrows,  # lower left corner
                           'width': panels_width/ncols*easy_aplpy.settings.colorbar_width,         # cbar width is specified fraction of panel width
                           'height': panels_height/nrows,                                          # cbar is as high as panels
                           'file': fitsfile,
                           'channel': None,
                           'physical': None
                           })
    return panels


###################################################################################################

def _show_map(fitsfile, fig, kwargs):
    cmap     = kwargs.get('cmap', 'viridis')                                   # the recommended cmap
    stretch  = kwargs.get('stretch', 'linear')
    vmin     = kwargs.get('vmin')                                              # no default, aplpy scales automatically
    vmax     = kwargs.get('vmax')                                              # no default, aplpy scales automatically
    imtype = kwargs.get('imtype')
    if ( imtype == 'pv'):
        aspect = 'auto'
    else:
        aspect = 'equal'
    fig.show_colorscale(cmap=cmap, vmin=vmin, vmax=vmax, stretch=stretch, aspect=aspect)


###################################################################################################

def _recenter_plot(fitsfile, fig, kwargs):
    recenter = kwargs.get('recenter')
    imtype = kwargs.get('imtype')
    if ( imtype == 'pp' ):
        if ( not recenter is None ) and _test_recenter_format(recenter):
            if (len(recenter) == 2):
                fig.recenter(recenter[0].ra.degree, recenter[0].dec.degree, radius=recenter[1].to(u.degree).value)
            elif (len(recenter) == 3):
                fig.recenter(recenter[0].ra.degree, recenter[0].dec.degree, width=recenter[1].to(u.degree).value, height=recenter[2].to(u.degree).value)
    if ( imtype == 'pv' ):
        if not ( recenter is None ):
            if ( len(recenter) == 4 ):
                header = fits.open(fitsfile)[0].header
                cunit1 = u.Quantity('1'+header['cunit1'])
                cunit2 = u.Quantity('1'+header['cunit2'])
                fig.recenter((recenter[0].to(cunit1.unit)).value, (recenter[1].to(cunit2.unit)).value, width=(recenter[2].to(cunit1.unit)).value, height=(recenter[3].to(cunit2.unit)).value)
            else:
                raise TypeError("Recenter: for a pV diagram specify [offset center, velocity center, width, height] with astropy.units.")


def _test_recenter_format(recenter):
    if not len(recenter) in [2,3]:
        raise TypeError("Recenter: specify SkyCoord(x,y) and either radius or width, height.")
    if not isinstance(recenter[0], SkyCoord):
        raise TypeError("Recenter position is not a SkyCoord object.")
    for x in recenter[1:]:
        if not isinstance(x, u.quantity.Quantity):
            raise TypeError("Recenter size argument(s) is not an astropy.units object.")
    return True


###################################################################################################

def _show_contours(fitsfile, fig, kwargs, panel=None):
    contours = kwargs.get('contours')
    clabel   = kwargs.get('clabel')
    legend   = kwargs.get('legend')
    if contours:
        if panel:
            contours = contours[panel['num']]                                  # get the correct set of contours for this panel
        if contours:
            contournum = 0                          # conting variable for # of contours
            for idx,cont in enumerate(contours):
                if len(cont) == 3:
                    fig.show_contour(data=cont[0], levels=cont[1], colors=cont[2])
                elif len(cont) == 4:
                    # two options when four arguments are given: slice argument (int) as second or kwargs (dict) as last element
                    if isinstance(cont[1],(int, np.int64, np.int32)):
                        fig.show_contour(data=cont[0], slices=[cont[1]], dimensions=[0,1], levels=cont[2], colors=cont[3])
                    elif type(cont[3]) is dict:
                        fig.show_contour(data=cont[0], levels=cont[1], colors=cont[2], **cont[3])
                    else:
                        raise TypeError("Contour: could not interpret contour list.")
                elif len(cont) == 5:
                    fig.show_contour(data=cont[0], slices=cont[1], dimensions=[0,1], levels=cont[2], colors=cont[3], **cont[4])
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


###################################################################################################

def _overplot_regions(fitsfile, fig, kwargs, panel=None):
    regions = kwargs.get('regions')
    if isinstance(regions,(list,tuple)):
        if panel:
            regions = regions[panel['num']]                                  # get the correct set of contours
        for region in regions:
            fig.show_regions(region)


###################################################################################################

def _show_colorbar(fitsfile, fig, kwargs):
    # get the header image unit if not defined by the user
    bunit = 'unknown quantity'
    if not kwargs.get('colorbar'):
        head = fits.open(fitsfile)[0].header
        if ( 'bunit' in head ):
            bunit = head['bunit']
        else:
            print("Header keyword 'BUNIT' not present in "+fitsfile)
    colorbar = kwargs.get('colorbar', ['right',bunit])
    stretch  = kwargs.get('stretch', 'linear')

    if not colorbar is None:
        fig.add_colorbar()
        fig.colorbar.show()
        fig.colorbar.set_location(colorbar[0])
        fig.colorbar.set_axis_label_text(colorbar[1])
        if ( stretch == 'log' ):
            im = fits.open(fitsfile)[0]
            vmin = kwargs.get('vmin', np.nanmin(im.data))
            vmax = kwargs.get('vmax', np.nanmax(im.data))
            log_ticks = [float('{:.2f}'.format(round(x,int(-1*np.log10(vmin))))) for x in np.logspace(np.log10(vmin),np.log10(vmax),num=10, endpoint=True)]
            fig.colorbar.set_ticks(log_ticks)
        fig.colorbar.set_font(size=easy_aplpy.settings.colorbar_ticks_fontsize)
        fig.colorbar.set_axis_label_font(size=easy_aplpy.settings.colorbar_label_fontsize)
        fig.colorbar.set_frame_color(easy_aplpy.settings.frame_color)


###################################################################################################

def _show_grid_colorbar(fitsfile, main_fig, fig, panels, kwargs):
    head = fits.open(fitsfile)[0].header
    if ( 'bunit' in head ):
        bunit = head['bunit']
    else:
        print("Header keyword 'BUNIT' not present in "+fitsfile)
        bunit = 'unknown quantity'
    colorbar = kwargs.get('colorbar', ['right',bunit])        # add the colorbar panel
    cmap     = kwargs.get('cmap', 'viridis')                                   # the recommended cmap
    stretch  = kwargs.get('stretch', 'linear')
    vmin     = kwargs.get('vmin')                                              # no default, aplpy scales automatically
    vmax     = kwargs.get('vmax')
    if not colorbar is None:
        cbpnl = panels[-1]
        ax1 = main_fig.add_axes([cbpnl['x'],cbpnl['y'],cbpnl['width'],cbpnl['height']])
        if ( colorbar[0] == 'last panel' ) or ( colorbar[0] == 'top' ):
            orientation = 'horizontal'
        elif ( colorbar[0] == 'right' ):
            orientation = 'vertical'
        else:
            raise NotImplementedError("Only colorbar in 'last panel' and 'right' of the last panel is supported at the moment.")

        if (stretch == 'linear'):
            mplcolorbar = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax), orientation=orientation)
        elif (stretch == 'log'):
            log_ticks = [float('{:.2f}'.format(x)) for x in np.logspace(np.log10(kwargs['vmin']),np.log10(kwargs['vmax']),num=5, endpoint=True)]
            mplcolorbar = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=mpl.colors.LogNorm(vmin=vmin, vmax=vmax), ticks=log_ticks, orientation=orientation)
            mplcolorbar.set_ticks(log_ticks)
            mplcolorbar.set_ticklabels(['{:.2f}'.format(x) for x in log_ticks])
        else:
            raise NotImplementedError("Scalings other than 'linear' and 'log' are not supported yet for grid plots.")

        mplcolorbar.set_label(colorbar[1])
        mplcolorbar.outline.set_edgecolor(easy_aplpy.settings.frame_color)
        mplcolorbar.ax.tick_params(labelsize=easy_aplpy.settings.colorbar_ticks_fontsize)

        from distutils.version import LooseVersion
        if ( LooseVersion(mpl.__version__) < LooseVersion('1.3') ):
            mplcolorbar.outline.set_color(easy_aplpy.settings.frame_color)
        else:
            mplcolorbar.outline.set_edgecolor(easy_aplpy.settings.frame_color)


###################################################################################################

def _show_scalebar(fitsfile, fig, kwargs, panel=None):
    scalebar = kwargs.get('scalebar')
    if isinstance(scalebar,list) and ( len(scalebar) == 3 ):
        if ( panel == None ) or ( ('left' in panel['position']) and ('bottom' in panel['position']) ):
            #TODO allow user to define in which panel the scalebar should be drawn
            fig.add_scalebar(length=scalebar[0].to(u.degree).value, label=scalebar[1], corner=scalebar[2], frame=easy_aplpy.settings.scalebar_frame)
            fig.scalebar.set_font(size=easy_aplpy.settings.scalebar_fontsize)
            fig.scalebar.set_linestyle(easy_aplpy.settings.scalebar_linestyle)
            fig.scalebar.set_linewidth(easy_aplpy.settings.scalebar_linewidth)
            fig.scalebar.set_color(easy_aplpy.settings.scalebar_color)


###################################################################################################

def _show_beam(fitsfile, fig, kwargs):
    beam = kwargs.get('beam', 'bottom right')
    imtype = kwargs.get('imtype')
    if not ( beam is None ) and not ( imtype == 'pv' ):
        fig.add_beam()
        fig.beam.show()
        fig.beam.set_corner(beam)
        fig.beam.set_frame(easy_aplpy.settings.beam_frame)
        fig.beam.set_color(easy_aplpy.settings.beam_color)


###################################################################################################

def _show_label(fitsfile, fig, kwargs):
    label = kwargs.get('label')
    if label:
        if isinstance(label, str):
            lbl = [[0.5,0.9],label, {}]
        elif isinstance(label,(list,tuple)):
            if ( len (label) == 2 ):
                lbl = [label[0], label[1], {}]
        fig.add_label(label[0][0], label[0][1], label[1].replace('_','$\_$'), color='black', relative=True, size=easy_aplpy.settings.velo_fontsize, **label[2])


###################################################################################################

def _show_overlays(fitsfile, fig, kwargs, panel=None):
    circles = kwargs.get('circles')
    if circles:
        if all(isinstance(x,(list,tuple)) for x in circles):
            if panel:
                circles = circles[panel['num']]                                # get the correct set of circles
            for circle in circles:
                fig.show_circles(xw=circle[0].ra.degree, yw=circle[0].dec.degree, radius=circle[1].to(u.degree).value, **circle[2])
        else:
            raise TypeError("Overlays: Must be list of lists. I.e. a single overlay needs double brackets [[]].")

    markers = kwargs.get('markers')
    if markers:
        if all(isinstance(x,(list,tuple)) for x in markers):
            if panel:
                markers = markers[panel['num']]                                # get the correct set of markers
            for marker in markers:
                fig.show_markers(xw=marker[0].ra.degree, yw=marker[0].dec.degree, **marker[1])
        else:
            raise TypeError("Overlays: Must be list of lists. I.e. a single overlay needs double brackets [[]].")

    polygons = kwargs.get('polygons')
    if polygons:
        if all(isinstance(x,(list,tuple)) for x in polygons):
            if panel:
                polygons = polygons[panel['num']]                              # get the correct set of polygons
            for polygon in polygons:
                fig.show_polygons(polygon_list=[np.array([[p.ra.degree, p.dec.degree] for p in polygon[0]])], **polygon[1])
        else:
            raise TypeError("Overlays: Must be list of lists. I.e. a single overlay needs double brackets [[]].")

    arrows = kwargs.get('arrows')
    if arrows:
        if all(isinstance(x,(list,tuple)) for x in arrows):
            if panel:
                arrows = arrows[panel['num']]                                # get the correct set of circles
            for arrow in arrows:
                fig.show_arrows(x=arrow[0].ra.degree, y=arrow[0].dec.degree, dx=arrow[1][0].to(u.degree).value, dy=arrow[1][1].to(u.degree).value, **arrow[2])
        else:
            raise TypeError("Overlays: Must be list of lists. I.e. a single overlay needs double brackets [[]].")

    ellipses = kwargs.get('ellipses')
    if ellipses:
        raise NotImplementedError("Overplotting ellipses is not supported yet.")

    lines = kwargs.get('lines')
    if lines:
        if all(isinstance(x,(list,tuple)) for x in lines):
            if panel:
                lines = lines[panel['num']]                              # get the correct set of lines
            for line in lines:
                fig.show_lines(line_list=[np.array([[i.value for i in line[0]], [j.value for j in line[1]]])], **line[2])
        else:
            raise TypeError("Overlays: Must be list of lists. I.e. a single overlay needs double brackets [[]].")

    rectangles = kwargs.get('rectangles')
    if rectangles:
        raise NotImplementedError("Overplotting rectangles is not supported yet.")

    texts = kwargs.get('texts')
    if texts:
        if all(isinstance(x,(list,tuple)) for x in texts):
            if panel:
                texts = texts[panel['num']]
            for text in texts:
                if not ( text == [] ):
                    if isinstance(text[0], SkyCoord):
                        x = text[0].ra.degree
                        y = text[0].dec.degree
                        relative = False
                    elif isinstance(text[0],(list,tuple)):
                        x = text[0][0]
                        y = text[0][1]
                        relative = True
                    else:
                        raise TypeError("Texts need to be a list of [position, text, kwargs] and position must be either a SkyCoord or a list of relative positions [x,y].")
                    fig.add_label(x=x, y=y, text=text[1], relative=relative, **text[2])
        else:
            raise TypeError("Overlays: Must be list of lists. I.e. a single overlay needs double brackets [[]].")


###################################################################################################

def _format_grid_ticksNlabels(panel, fig, kwargs):
    imtype = kwargs.get('imtype')
    if ( imtype == 'pp' ):
        fig.ticks.show()
        fig.ticks.set_xspacing((easy_aplpy.settings.ticks_xspacing).to(u.degree).value)
        fig.ticks.set_yspacing((easy_aplpy.settings.ticks_yspacing).to(u.degree).value)
        fig.ticks.set_minor_frequency(easy_aplpy.settings.ticks_minor_frequency)
        fig.ticks.set_color(easy_aplpy.settings.ticks_color)
        fig.frame.set_color(easy_aplpy.settings.frame_color)
        fig.axis_labels.set_font(size=easy_aplpy.settings.tick_label_fontsize)
        fig.axis_labels.hide()
        fig.tick_labels.hide()
    if ( imtype == 'pv' ):
        labels = kwargs.get('labels')
        if labels:
            fig.set_axis_labels(labels[0],labels[1])
        fig.ticks.show()
        fig.ticks.set_minor_frequency(easy_aplpy.settings.ticks_minor_frequency)
        fig.ticks.set_color(easy_aplpy.settings.ticks_color)
        fig.frame.set_color(easy_aplpy.settings.frame_color)
        fig.axis_labels.set_font(size=easy_aplpy.settings.tick_label_fontsize)
        fig.axis_labels.hide()
        fig.tick_labels.hide()

    if ( 'left' in panel['position'] ):
        fig.axis_labels.show_y()
        fig.tick_labels.show_y()
    if ( 'bottom' in panel['position'] ):
        fig.axis_labels.show_x()
        fig.tick_labels.show_x()


###################################################################################################

def _show_ticksNlabels(fitsfile, fig, kwargs):
    imtype = kwargs.get('imtype')
    if ( imtype == 'pp' ):
        fig.tick_labels.show()
        fig.tick_labels.set_xformat(easy_aplpy.settings.tick_label_xformat)
        fig.tick_labels.set_yformat(easy_aplpy.settings.tick_label_yformat)
        fig.tick_labels.set_font(size=easy_aplpy.settings.tick_label_fontsize)
        fig.ticks.show()
        fig.ticks.set_xspacing((easy_aplpy.settings.ticks_xspacing).to(u.degree).value)
        fig.ticks.set_yspacing((easy_aplpy.settings.ticks_yspacing).to(u.degree).value)
        fig.ticks.set_minor_frequency(easy_aplpy.settings.ticks_minor_frequency)
        fig.ticks.set_color(easy_aplpy.settings.ticks_color)
        fig.frame.set_color(easy_aplpy.settings.frame_color)
        fig.axis_labels.set_font(size=easy_aplpy.settings.tick_label_fontsize)
    if ( imtype == 'pv' ):
        labels = kwargs.get('labels')
        if labels:
            fig.set_axis_labels(labels[0],labels[1])
        fig.tick_labels.show()
        fig.tick_labels.set_font(size=easy_aplpy.settings.tick_label_fontsize)
        fig.ticks.show()
        fig.ticks.set_minor_frequency(easy_aplpy.settings.ticks_minor_frequency)
        fig.ticks.set_color(easy_aplpy.settings.ticks_color)
        fig.frame.set_color(easy_aplpy.settings.frame_color)
        fig.axis_labels.set_font(size=easy_aplpy.settings.tick_label_fontsize)


###################################################################################################

def _show_legend(fitsfile, fig, kwargs):
    legend = kwargs.get('legend')
    if legend:
        fig._ax1.legend(loc=0, fontsize=easy_aplpy.settings.colorbar_label_fontsize)
        if isinstance(legend,dict):
            fig._ax1.legend(fontsize=easy_aplpy.settings.colorbar_label_fontsize, **legend)


###################################################################################################

def _show_channel_label(panel, fig, kwargs):
    channel_label = kwargs.get('channel_label','physical')
    if isinstance(channel_label,str):
        if ( channel_label == 'physical' ):
            label = ((easy_aplpy.settings.grid_label_format).format(panel['physical'])).replace('km / s','km\,s$^{-1}$')
        elif ( channel_label == 'number' ):
            label = '{:d}'.format(panel['channel'])
    elif ( channel_label == None ):
        label = False
    elif isinstance(channel_label,(list,tuple)):
        label = str(channel_label[panel['num']])
    else:
        raise TypeError("Unrecognized type of channel_label. Must be an instruction ('physical', 'number', None) or a list of strings.")

    if label:
        fig.add_label(easy_aplpy.settings.grid_label_pos[0],
            easy_aplpy.settings.grid_label_pos[0],
            label,
            color    = easy_aplpy.settings.grid_label_color,
            relative = True,
            size     = easy_aplpy.settings.grid_label_fontsize
            )


###################################################################################################

def _execute_code(fitsfile, fig, kwargs, panel=None):
    execute_code = kwargs.get('execute_code')
    if execute_code:
        if isinstance(execute_code, (list,tuple)):
            if panel:
                execute_code = execute_code[panel['num']]
            for codes in execute_code:
                exec(codes) in locals()
        else:
            raise TypeError("Execute code: Code to execute must be given in a list of strings")


###################################################################################################

def _save_figure(fitsfile, fig, kwargs):
    out = kwargs.get('out',os.path.splitext(fitsfile)[0]+'.png')
    if isinstance(out,str):
        fig.savefig(out, dpi=300, transparent=True, adjust_bbox=True)
        print("\x1b[0;34;40m[easy_aplpy]\x1b[0m saved plot as "+out)
    if isinstance(out,dict):
        fig.savefig(**out)
        print("\x1b[0;34;40m[easy_aplpy]\x1b[0m saved plot as "+out['filename'])


###################################################################################################
