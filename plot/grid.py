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
