# easy_aplpy
Baking a publication quality plot can involve a fair bit of code. Easy APLpy can make this process much simpler by offering wrapper functions for the most often used types of plots.

For details on APLpy see https://github.com/aplpy/aplpy

## available wrappers:

### easy_aplpy.plot.map
A one-liner to quickly plot fits image that also look good. This method also provides extensive options to costumize the plot but still keep the code simple.

Running this code:
```
import easy_aplpy
easy_aplpy.plot.map('test_data/map.fits')
```
Produces this plot:
![simple map](https://github.com/GiantMolecularCloud/easy_aplpy/blob/master/test_data/map.png)


### easy_aplpy.plot.grid
A function to easily plot a grid of images, e.g. the channel maps that are often often in radio astronomy.

Running this code:
```
import easy_aplpy
easy_aplpy.plot.grid('cube.fits', [2,3], [150,200,250,300,350,400])
```
Produces this plot:
![simple channel map](https://github.com/GiantMolecularCloud/easy_aplpy/blob/master/test_data/cube.png)



## Installation:
- clone this repository or download the scripts
- add the path where the file are located to your python path
    ```
    import sys
    sys.path.append('/your/path/to/easy_aplpy')
    ```
- import the wrapper functions
    ```
    import easy_aplpy
    ```
- call the plotting wrappers with `easy_aplpy.plot`


# more details coming soon ...
