from __future__ import print_function
#####################################################################
#                          APLPY PLOTTING                           #
#####################################################################
# Helper functions to accomodate plotting in a nice way.            #
#####################################################################

__all__ = ['hide_deprecationWarnings','hide_nonfunctionalWarnings','hide_FITSwarnings','hide_ComparisonWarnings','m_to_km','trim_whitespace']


###################################################################################################

# common imports
################

import os
import subprocess
from astropy.io import fits


###################################################################################################

# switch of deprecation warnings
################################

def hide_deprecationWarnings():
    try:
        from matplotlib.cbook import MatplotlibDeprecationWarning
        import warnings
        warnings.simplefilter('ignore', MatplotlibDeprecationWarning)
    except:
        raise Warning("easy_aplpy: Could not switch off MatplotlibDeprecationWarnings.These are raised because of aplpy but cause no harm.")


def hide_nonfunctionalWarnings():
    try:
        import warnings
        warnings.filterwarnings("ignore", message="This method is not functional at this time")
        warnings.filterwarnings("ignore", message="WARNING: Requested tick spacing format cannot be shown by current label format. The tick spacing will not be changed.")
        warnings.filterwarnings("ignore", message="No labelled objects found.")
    except:
        raise Warning("easy_aplpy: Could not switch off aplpy warnings about non-functional methods. They do NOT affect easy_aplpy.")


def hide_FITSwarnings():
    try:
        import warnings
        warnings.filterwarnings("ignore", message=".*indices in parameterized keywords must not have leading zeroes.*")
    except:
        raise Warning("easy_aplpy: Could not switch off aplpy warnings about leading zeros in FITS header. They do NOT affect plotting at all.")


def hide_ComparisonWarnings():
    try:
        import warnings
        warnings.filterwarnings("ignore", message="invalid value encountered in less")
        warnings.filterwarnings("ignore", message="Unicode equal comparison failed to convert both arguments to Unicode - interpreting them as being unequal")
        warnings.filterwarnings("ignore", message="comparison to `None` will result in an elementwise object comparison in the future.")
    except:
        raise Warning("easy_aplpy: Could not switch off aplpy comparison warnings. These are raised within aplpy and do NOT affect plotting at all.")


###################################################################################################

# correct velocity header
#########################

def m_to_km(pv, overwrite=False, out=None):

    """
    convert_velo_info: scale the velocity axis in pV diagrams
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    CASA position-velocity diagrams are exported in m/s which is inconvenient
    because the tick labels get very long. This function converts the velocity
    axis to km/s by changing the header.
    !!! Note that aplpy in order to display the velocity correctly requires
    !!! overcorrecting. This seems to be a bug somewhere in aplpy or the FITS
    !!! images. The resulting images are off by a factor of 1000, 100km/s are
    !!! given as 0.1km/s in the FITS header but aplpy shows it as 100km/s.
    It is recommended to use the "corrected" images for plotting only.

    Mandatory arguments:
        pv          Path and file name of a single pV diagram as FITS file or a list
                    of multiple FITS images.

    Optional arguments:
        overwrite   True or False. Default: True
        out         Image name or list of names for output. Default: None
                    If neither overwrite or out are given ".corrected" is appended
                    to the file name.

    example:

    correct_velo_info(['file1.fits', 'file2.fits', file3.fits'],
        overwrite = False
        )
    """

    out = kwargs.get('out', pv)
    if isinstance(pv, str) and isinstance(out, str):
        pv_list = [pv]
        out_list = [out]
    elif isinstance(pv, (list,tuple)) and isinstance(out, (list,tuple)):
        pv_list = pv
        out_list = out
    else:
        raise TypeError('Input needs to be a single FITS file or list of FITS files. (And out parameter accordingly.)')

    for this_pv, this_out in zip(pv_list,out_list):

        # get current velocity unit
        velounit = fits.open(this_pv)[0].header['cunit2']

        # convert to km/s if necessary
        if (velounit == 'm/s'):
            if (overwrite == False):
                if not ( this_out == None ):
                    os.system('cp -r '+this_pv+' '+this_out)
                    this_pv = this_out
                else:
                    os.system('cp -r '+this_pv+' '+this_pv+'.corrected')
                    this_pv = this_pv+'.corrected'

            print('Unit is m/s. Corrrecting to km/s: '+this_pv)
            im = fits.open(this_pv)[0]
            im.header['cdelt2'] /= 1e6
            im.header['crval2'] /= 1e6
            im.header['cunit2'] = 'km/s'
            fits.writeto(this_pv, data=im.data, header=im.header, overwrite=True)
            print('\tNote that the header is now wrong by a factor of 1e3!')
        elif (velounit == 'km/s\n'):
            print('Unit is km/s already: '+this_pv)
        else:
            raise TypeError('Unrecognized velocity unit "'+str(velounit)+'": '+this_pv)


###################################################################################################

# trim whitespace in images
###########################

def trim_whitespace(images, dpi=300):
    """Short summary.

    Parameters
    ----------
    images : str or list
        Single file name or list of file names of the images to be trimmed.
    dpi : int
        The image resolution (caled 'density' in imagemagick). Defaults to 300.
    """

    import random

    # try to import tqdm to get a nice progress bar
    try:
        import tqdm
        check_tqdm = True
    except:
        check_tqdm = False

    # check input(s)
    if isinstance(images, str):
        imlist = [images]
    elif isinstance(images, (list,tuple)):
        imlist = images
    else:
        raise TypeError("Input needs to be a file name or list thereof.")


    # the function that does the actual work
    def imagemagick_convert_trim(img):
        # generate temporary file name
        tmp = int(1e6*random.random())
        name, fext = os.path.splitext(img)
        tmp_fname = 'tmp_'+str(tmp)+fext

        try:
            os.system('convert -density '+str(dpi)+' -fuzz 1% -trim +repage '+img+' '+tmp_fname)
            os.system('rm -f '+img)
            os.system('cp '+tmp_fname+' '+img)
            os.system('rm -rf '+tmp_fname)
        except:
            raise Exception("`convert` as part of the Imagemagick package is not installed or not recognized.\nMake sure the correct program reacts to `convert` in the default shell.")


    # use the tqdm progress bar if possible, otherwise count up
    if check_tqdm:
        for img in tqdm.tqdm(imlist):
            imagemagick_convert_trim(img)
    else:
        for idx,img in enumerate(imlist):
            print("triming image "+str(idx+1)+" of "+str(len(imlist))) #, end='\r')
            imagemagick_convert_trim(img)
        print("\n")



###################################################################################################
