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

    ...
    """

    kwargs = _check_image_type(fitsfile, kwargs)
    main_fig = _set_up_grid(fitsfile, shape, kwargs)
    panels = _grid_panels(fitsfile, shape, channels, kwargs)

    for panel in panels:

        fig = _set_up_panel_figure(main_fig, panel, kwargs)
        _show_map(panel['file'], fig, kwargs)
        _recenter_plot(panel['file'], fig, kwargs)
        _show_beam(panel['file'], fig, kwargs)
        _format_grid_ticksNlabels(panel, fig, kwargs)
        #_show_channel_label(panel, fig, kwargs)
        _show_panel_contours(panel, fig, kwargs)
        #_overplot_panel_regions(panel, fig, kwargs)
        #_show_panel_overlays(panel, fig, kwargs)
        #_show_grid_scalebar(panel, fig, kwargs)
        #_show_legend(panel['file'], fig, kwargs)
        #_execute_code(panel['file'], fig, kwargs)

    #_show_grid_colorbar(panel, fig, kwargs)
    _save_figure(fitsfile, main_fig, kwargs)

    return main_fig


###################################################################################################
