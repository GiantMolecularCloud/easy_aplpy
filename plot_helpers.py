#####################################################################
#                          APLPY PLOTTING                           #
#####################################################################
# These functions will produce plots of channel maps, moment maps,  #
# pV diagrams, ... in a quality that (hopefully) allows publishing. #
#####################################################################

__all__ = ['test_recenter_format']


###################################################################################################

# common import
###############

import os
import aplpy
import numpy
from astropy.coordinates import SkyCoord as SkyCoord
from astropy import units as u
from astropy.coordinates import Angle as Angle
from matplotlib import rc as rc
rc('text',usetex=True)


###################################################################################################

# plot helper functions
#######################

def test_recenter_format(recenter):
    """
    Test if a given recenter list of correct shape and type.

    Parameters
    ----------
    recenter : recenter list

    Returns
    -------
    Bool
    """
    if not len(recenter) in [2,3]:
        raise TypeError("Recenter: specify SkyCoord(x,y) and either radius or width, height.")
    if not isinstance(recenter[0], SkyCoord):
        raise TypeError("Recenter position is not a SkyCoord object.")
    for x in recenter[1:]:
        if not isinstance(x, u.quantity.Quantity):
            raise TypeError("Recenter size argument(s) is not an astropy.units object.")

    return True


###################################################################################################
