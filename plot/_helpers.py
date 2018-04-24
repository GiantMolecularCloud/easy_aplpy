#####################################################################
#                          APLPY PLOTTING                           #
#####################################################################
# These functions will produce plots of channel maps, moment maps,  #
# pV diagrams, ... in a quality that (hopefully) allows publishing. #
#####################################################################

__all__ = ['_check_image_type','_set_up_figure','_show_map','_recenter_plot','_test_recenter_format','_show_contours','_overplot_regions','_show_colorbar','_show_scalebar','_show_beam','_show_label','_show_overlays','_show_ticksNlabels','_show_legend','_execute_code','_save_figure']


###################################################################################################

# common imports
################

import os
import aplpy
import numpy as np
from astropy.coordinates import SkyCoord as SkyCoord
from astropy import units as u
from astropy.coordinates import Angle as Angle
from astropy.io import fits
from matplotlib import rc as rc
rc('text',usetex=True)

import easy_aplpy


###################################################################################################
# plotting subfunctions
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
    figsize = kwargs.get('figsize', (8.267,11.692))     # A4 in inches
    channel = kwargs.get('channel', None)
    if ( channel == None ):
        fig = aplpy.FITSFigure(fitsfile, figsize=figsize)
    else:
        if isinstance(channel, int):
            pass
        elif isinstance(channel, u.quantity.Quantity):
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
            channel = int(((channel-crval)/cdelt)+crpix)                     # convert to channel number
        else:
            raise TypeError("channel needs to be channel number (int) or velocity/frequency (astropy.units object.)")
        fig = aplpy.FITSFigure(fitsfile, slices=[channel], figsize=figsize)
    return fig


###################################################################################################

def _show_map(fitsfile, fig, kwargs):
    cmap     = kwargs.get('cmap', 'viridis')             # the recommended cmap
    stretch  = kwargs.get('stretch', 'linear')
    vmin     = kwargs.get('vmin')                        # no default, aplpy scales automatically
    vmax     = kwargs.get('vmax')                        # no default, aplpy scales automatically
    fig.show_colorscale(cmap=cmap, vmin=vmin, vmax=vmax, stretch=stretch)


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

def _show_contours(fitsfile, fig, kwargs):
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


###################################################################################################

def _overplot_regions(fitsfile, fig, kwargs):
    regions = kwargs.get('regions')
    if isinstance(regions,(list,tuple)):
        for region in regions:
            fig.show_regions(region)

###################################################################################################

def _show_colorbar(fitsfile, fig, kwargs):
    colorbar = kwargs.get('colorbar', ['right',fits.open(fitsfile)[0].header['bunit']])
    stretch  = kwargs.get('stretch', 'linear')
    if not colorbar is None:
        fig.add_colorbar()
        fig.colorbar.show()
        fig.colorbar.set_location(colorbar[0])
        fig.colorbar.set_axis_label_text(colorbar[1])
        if ( stretch == 'log' ):
            log_ticks = [float('{:.2f}'.format(round(x,int(-1*np.log10(kwargs['vmin']))))) for x in np.logspace(np.log10(kwargs['vmin']),np.log10(kwargs['vmax']),num=10, endpoint=True)]
            fig.colorbar.set_ticks(log_ticks)
        fig.colorbar.set_font(size=easy_aplpy.settings.colorbar_fontsize)
        fig.colorbar.set_axis_label_font(size=easy_aplpy.settings.colorbar_fontsize)
        fig.colorbar.set_frame_color(easy_aplpy.settings.frame_color)


###################################################################################################

def _show_scalebar(fitsfile, fig, kwargs):
    scalebar = kwargs.get('scalebar')
    if isinstance(scalebar,list) and ( len(scalebar) == 3 ):
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

def _show_overlays(fitsfile, fig, kwargs):
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
        fig._ax1.legend(loc=0, fontsize=easy_aplpy.settings.colorbar_fontsize)
        if isinstance(legend,dict):
            fig._ax1.legend(fontsize=easy_aplpy.settings.colorbar_fontsize, **legend)


###################################################################################################

def _execute_code(fitsfile, fig, kwargs):
    execute_code = kwargs.get('execute_code')
    if execute_code:
        if isinstance(execute_code, (list,tuple)):
            for codes in execute_code:
                exec(codes) in locals()
        else:
            raise TypeError("Execute code: Code to execute must be given in a list of strings")


###################################################################################################

def _save_figure(fitsfile, fig, kwargs):
    out = kwargs.get('out',os.path.splitext(fitsfile)[0]+'.png')
    if isinstance(out,str):
        fig.save(out, dpi=300, transparent=True, adjust_bbox=True)
        print("\x1b[0;34;40m[easy_aplpy]\x1b[0m saved plot as "+out)
    if isinstance(out,dict):
        fig.save(**out)
        print("\x1b[0;34;40m[easy_aplpy]\x1b[0m saved plot as "+out['filename'])


###################################################################################################
