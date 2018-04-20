# easy_aplpy
Baking a publication quality plot can involve a fair bit of code. Easy APLpy can make this process much simpler by offering wrapper functions for the most often used types of plots.

For details on APLpy see https://github.com/aplpy/aplpy



~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# old description below, needs to be rewritten
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## available wrappers:

### simple map plotting `ap.aplpy_plot`
This is a one-liner to quickly plot a fits file that consists of only on channel/plane. It just needs the name of a fits file to render a figure without having to setup the aplpy figure using several lines of code. Further options are available to customize the plot.

### simple channel plotting `ap.aply_plot_slice`
Basically the same as `ap.aplpy_plot` but for image cube with multiple channels. Beside the name of the file to plot, it also needs the number of the channel that should be plotted.

### channel map `ap.aplpy_channel_map`
Wrapper to plot "channel maps": a grid of plots showing a selection of channels of an image cube. File name, number of rows and columns, first channel to plot and step width need to be given. Further options allow customisation.
The color bar will always be plotted in the last (bottom right) panel, axis labels always in the bottom left panel. This may cause the labels to be missing when to few panels are plotted (e.g. 6 channels in a 3x3 grid).

### map grid `ap.aplpy_map_grid`
Plot a list of single plane fits file in a grid. Rerquires a list of file names and number of rows/columns.

### position-velocity diagram `ap.aply_plot_pv`
Position velocity diagrams slice a 3D image cube along the spectral axis and are sort of a spectrum along a line instead of just a single pixel. This functions sets up the figure correctly, so aplpy can understand that one axis of the image is not spatial but a spectral quantity. Needs just a file name to render a plot.

### position-velocity grid `ap.aply_pv_grid`
A combination of "map grid" and "pv diagram": Plot a grid of multiple pV diagrams side-by-side to allow direct comparison. Requires a list of file names and number of columns/rows to plot.

## Plot customisartion
These wrappers were originally written to simplify my plotting scripts to contain only functions with a handful of arguments instead of hundreds of lines of aplpy code. They may not do exactly what you want to do, so please point out bugs or where functionality is missing.

Each function offers the most often used customisation options directly through arguments. For details see the inline help that can be called using a `?`, e.g. `?ap.aplpy_plot`.

Beside the options that are directly available through arguments, most functions also contain a `execute_code` argument. Strings passed to `execute_code` will be executed within the wrapper function which allows you to directly reach the figure's internals like the figure object or even the underlying matplotlib objects. Through this argument you can modify the figure in any way but still keep your plotting code contained within one relatively short function.
Be careful with this argument. You may crash the wrapper function by redefining objects.


## Installation:
- clone this repository or download the scripts
- add the path where the file are located to your python path
    ```
    import sys
    sys.path.append('/your/path/to/aplpy_plotting')
    ```
- import the wrapper functions
    ```
    import aplpy_plotting as ap
    ```
- call the plotting wrappers with `ap.xxxxx`

## known issues:
- Not all possible combinations of keyword arguments work (e.g. giving just vmin and vmax but not cmap).
- Channel maps: first channel to plot and step need to be given as channel number, giving velocity or frequency is not supported yet.
- Channel map and map grids: color map has to be given twice. Once for using the correct one while plotting the panels and once again to draw the color bar in the last panel. Will be fixed soon!
- The functions sometimes offer more parameters than listed in the inline help (`?ap.aplpy_channel_map`). These are options that I introduced for specific plots but didn't have the time yet to update the docstring.

## some examples:

### A simple figure
```
aplpy.plot('map.fits')
```
![moment map](http://www2.mpia-hd.mpg.de/homes/krieger/images/SWAG_moment_map.png)

### A channel map
It's also very easy to get a channel map with the function aplpy_channel_map. No need to type ~100 lines of code anymore.
```
ap.aplpy_channel_map('datacube.fits',
    2,      # how many columns?
    5,      # how many rows?
    65,     # first channel to plot
    15,     # plot every n-th channel
    cmap            = 'viridis',
    contour         = ['pixelmask', 'mask.fits', [0.5], 'grey'],
    scalebar_length = 0.3452,
    scalebar_label  = '50\,pc',
    scalebar_corner = 'bottom'
    )
```
![channel map](http://www2.mpia-hd.mpg.de/homes/krieger/images/SWAG_channelmap.png)

### A position velocity diagram
It gets a bit more tricky when not both axis are spatial coordinates but aplpy_plot_pv does the job for you.
```
aplpy_plot_pv('pv_file.fits',
    xlabel = 'Galactic Longitude',
    ylabel = 'optical velocity [km\,s$^{-1}$]'
    )
```
![pV diagram](http://www2.mpia-hd.mpg.de/homes/krieger/images/SWAG_pV-l.png)
