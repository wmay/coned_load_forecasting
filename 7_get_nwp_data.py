#!/<path-to-venv>/bin/python

'''Download and organize NWP forecasts.

This script relies on the wgrib2 system package, which can be installed with
spack or conda. (I recommend spack.)
'''

# pip install 'herbie-data[extras]'
import os
import pandas as pd
import xarray as xr
from herbie import Herbie, wgrib2
from rechunker import rechunk

# first step: download files with Herbie, subset with wgrib2, organize

# to define a dataset, we need: model, parameters, extent

def download_extent_only(archive, params, extent, name='region'):
    '''Download a subset of a grib2 file, and spatially subset the file.
    '''
    grib_file = archive.download(params, verbose=True)
    subset_file = wgrib2.region(grib_file, extent, name=name)
    os.remove(grib_file) # delete larger file
    return subset_file


# GEFS test run

# https://herbie.readthedocs.io/en/latest/gallery/noaa_models/gefs.html

# another complication for GEFS is that there are both .5 and .25 degree
# products

archive = Herbie("2023-01-01", model='gefs', save_dir='nwp_data', member=0)

archive2 = Herbie("2023-01-01", model='gefs', save_dir='nwp_data', member=1)

# figure out the search strings
inv1 = archive.inventory()

rasp_params = pd.read_csv('nwp_code/config/gefs_rasp.csv')
# concatenate the search strings
search_string = '|'.join(rasp_params['search'][:2])

# define the spatial extent for subsetting
project_extent = (285.5, 286.5, 40, 41.5)

f1 = download_extent_only(archive, search_string, project_extent, name='nyc')

grib_file = archive.download(r"TMP:2 m above ground", verbose=True)



# second step: modify as needed


# third step: combine with rechunker
