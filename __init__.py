#####################################################################
#                          APLPY PLOTTING                           #
#####################################################################

"""
These functions will produce plots of channel maps, moment maps,
pV diagrams, ... in a quality that (hopefully) allows publishing.
"""


###################################################################################################

# import modules
################

from . import helpers
from . import plot
from . import settings
from . import custom_colormaps

###################################################################################################

# set default settings
######################

helpers.hide_deprecationWarnings()
helpers.hide_nonfunctionalWarnings()
helpers.hide_FITSwarnings()
helpers.hide_ComparisonWarnings()


###################################################################################################
