import numpy as np
import matplotlib as mpl
import matplotlib.colors as colors
import matplotlib.cm as cm

__all__ = ['viridis_cropped']

###################################################################################################

# define custom colormaps
#########################

# viridis with less indistinguishable dark colors
if mpl.__version__ > '1.5.0':
    viridis_cropped = colors.ListedColormap(cm.viridis(np.linspace(0.1,1.0,100)))
else:
    print("Your version of Matplotlib does not contain the viridis colormap. Please update if possible.")
    viridis_cropped = 'Your version of Matplotlib does not contain the viridis colormap. Please update if possible.'


###################################################################################################
