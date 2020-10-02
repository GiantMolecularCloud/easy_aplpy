###################################################################################################
#                                           EASY APLPY                                            #
###################################################################################################

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

# warn that APLpy is a mess
###########################

import aplpy
if aplpy.version.major!=1:
    print("\x1b[0;31;40m"+"easy_aplpy only partially works with APLpy 2 and the old APLpy 0 version. Some functionality will raise unsolvable errors."+"\x1b[0m")
# if aplpy.version.version=='1.1.1':
#     print("\x1b[0;31;40m"+"APLpy 1.1.1 has problems plotting the beam but newer and older version have even more serious bugs."+"\x1b[0m")


###################################################################################################

# easy_aplpy version
####################

version = '0.2'


###################################################################################################
