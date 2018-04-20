from astropy.coordinates import Angle as Angle
from astropy import units as u

__all__ = ['tick_label_xformat','tick_label_yformat','ticks_xspacing','ticks_yspacing','ticks_minor_frequency','velo_fontsize','colorbar_fontsize','colorbar_width','scalebar_frame','scalebar_linestyle','scalebar_linewidth','scalebar_color','scalebar_fontsize','beam_frame','beam_color','ticks_color','frame_color','tick_label_fontsize','axis_label_fontsize','props']

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
velo_fontsize         = 10.0		# unit: point
colorbar_fontsize     = 10.0		# unit: point
colorbar_width        = 0.15	    # relative to panel size
scalebar_frame        = False
scalebar_linestyle    = 'solid'	    # or any other plt.plot linestyle
scalebar_linewidth    = 2			# unit: points
scalebar_color        = 'red'		# any named color or mpl.color instance
scalebar_fontsize     = 10.0    	# only used in channel map to prevent bar sliding over map
beam_frame            = False
beam_color            = 'black'
ticks_color           = 'black'	    # this setting overrules the matplotlibrc defaults
frame_color           = 'black'
tick_label_fontsize   = 12          # unit: point
axis_label_fontsize   = 12          # unit: point
props                 = {'boxstyle': "round", 'facecolor': "w", 'edgecolor': "black", 'linewidth': 0.5, 'alpha': 0.8}
