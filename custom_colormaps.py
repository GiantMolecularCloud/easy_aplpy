import numpy as np
import matplotlib.colors as colors
import matplotlib.cm as cm

__all__ = ['viridis_cropped']

###################################################################################################

# define custom colormaps
#########################

# viridis with less indistinguishable dark colors
viridis_cropped = colors.ListedColormap(cm.viridis(np.linspace(0.1,1.0,100)))


###################################################################################################
