from astropy.coordinates import Angle as Angle
from astropy import units as u

__all__ = ['tick_label_xformat','tick_label_yformat','ticks_xspacing','ticks_yspacing','ticks_minor_frequency','colorbar_fontsize','colorbar_width','scalebar_frame','scalebar_linestyle','scalebar_linewidth','scalebar_color','scalebar_fontsize','beam_frame','beam_color','ticks_color','frame_color','tick_label_fontsize','axis_label_fontsize','grid_label_pos','grid_label_color','grid_label_fontsize','grid_label_format','margins','props']

###################################################################################################

# settings that are dataset dependend, i.e. need to be changed for a new dataset
# these are just some defaults
tick_label_xformat    = 'hh:mm:ss.s'
tick_label_yformat    = 'dd:mm:ss'
ticks_xspacing        = Angle('0 1 0', unit='hourangle')
ticks_yspacing        = 10.0*u.arcsec
ticks_minor_frequency = 5

# fixed settings, i.e. settings that should be the same for each and every plot
# these can still be changed if necessary
colorbar_fontsize     = 10      		# unit: point
colorbar_width        = 0.15	        # relative to panel height
scalebar_frame        = False
scalebar_linestyle    = 'solid'	        # or any other plt.plot linestyle
scalebar_linewidth    = 2			    # unit: point
scalebar_color        = 'red'		    # any named color or mpl.color instance
scalebar_fontsize     = 10.0
beam_frame            = False           # plot a frame around the beam ellipse
beam_color            = 'black'
ticks_color           = 'black'	        # this setting overrules the matplotlibrc defaults
frame_color           = 'black'
tick_label_fontsize   = 12              # unit: point
axis_label_fontsize   = 12              # unit: point
grid_label_pos        = [0.8,0.8]       # relative position in panel
grid_label_color      = 'black'
grid_label_fontsize   = 10              # unit: point
grid_label_format     = '{:3.1f}'       # .format string
margins               = [0.1,0.05,0.05,0.05]    # margins around figure (left, right, top, bottom)
props                 = {'boxstyle': "round", 'facecolor': "w", 'edgecolor': "black", 'linewidth': 0.5, 'alpha': 0.8}
