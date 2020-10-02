#####################################################################
#                            EASY APLPY                             #
#####################################################################
# These functions will produce plots of channel maps, moment maps,  #
# pV diagrams, ... in a quality that (hopefully) allows publishing. #
#####################################################################

###################################################################################################

from astropy.coordinates import Angle as Angle
from astropy import units as u

__all__ = ['tick_label_xformat',
           'tick_label_yformat',
           'tick_labelpad',
           'ticks_xspacing',
           'ticks_yspacing',
           'ticks_minor_frequency',
           'colorbar_label_fontsize',
           'colorbar_width',
           'colorbar_labelpad',
           'scalebar_frame',
           'scalebar_linestyle',
           'scalebar_linewidth',
           'scalebar_color',
           'scalebar_fontsize',
           'beam_kwargs',
           'ticks_color',
           'frame_color',
           'tick_label_fontsize',
           'axis_label_fontsize',
           'tick_length',
           'grid_label_pos',
           'grid_label_color',
           'grid_label_fontsize',
           'grid_label_format',
           'grid_label_all',
           'grid_label_bbox',
           'margins',
           'props'
          ]


###################################################################################################

# settings that are dataset dependend, i.e. need to be changed for a new dataset
# these are just some defaults
tick_label_xformat    = 'hh:mm:ss.s'
tick_label_yformat    = 'dd:mm:ss'
tick_labelpad         = 0
ticks_xspacing        = Angle('0 1 0', unit='hourangle')
ticks_yspacing        = 10.0*u.arcsec
ticks_minor_frequency = 5

# fixed settings, i.e. settings that should be the same for each and every plot
# these can still be changed if necessary
colorbar_label_fontsize = 10      		# unit: point
colorbar_width        = 0.10	        # relative to panel height
colorbar_labelpad     = 0
scalebar_frame        = False
scalebar_linestyle    = 'solid'	        # or any other plt.plot linestyle
scalebar_linewidth    = 2			    # unit: point
scalebar_color        = 'red'		    # any named color or mpl.color instance
scalebar_fontsize     = 10.0
beam_kwargs           = {'loc': 'lower left', 'facecolor': 'k', 'edgecolor': 'k', 'frame': True, 'pad': 0.5, 'borderpad': 1.0, 'filled': True, 'linewidth': 1.0, 'zorder': 99}
ticks_color           = 'black'	        # this setting overrules the matplotlibrc defaults
frame_color           = 'black'
tick_label_fontsize   = 12              # unit: point
axis_label_fontsize   = 12              # unit: point
tick_length           = 8               # unit: point
grid_label_pos        = [0.95,0.95]     # relative position in panel
grid_label_color      = 'black'
grid_label_fontsize   = 10              # unit: point
grid_label_format     = '{:3.1f}'       # .format string
grid_label_all        = False           # label all panels or just the bottom left?
grid_label_bbox       = False           # show a box (see props below) around the labels in grid
margins               = [0.10,0.10,0.05,0.10]    # margins around figure (left, right, top, bottom)
props                 = {'boxstyle': "round", 'facecolor': "white", 'edgecolor': "black", 'linewidth': 0.5, 'alpha': 0.8}
