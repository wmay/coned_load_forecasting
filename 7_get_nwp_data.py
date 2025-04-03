'''Download and organize NWP forecasts.

This script relies on the wgrib2 system package, which can be installed with
spack or conda. (I recommend spack.)
'''

# pip install 'herbie-data[extras]'
import os, glob
import numpy as np
import pandas as pd
import xarray as xr
#from herbie import Herbie, FastHerbie, wgrib2
from herbieplus.herbie import FastRegionHerbie, EnsHerbie
from herbieplus.xarray import open_herbie_dataset
#from rechunker import rechunk

# first step: download files with Herbie, subset with wgrib2, organize

# GEFS test run

# https://herbie.readthedocs.io/en/latest/gallery/noaa_models/gefs.html

# another complication for GEFS is that there are both .5 and .25 degree
# products

# parameters from Rasp and Lerch 2018
rasp_params = pd.read_csv('nwp_code/config/gefs_rasp.csv')

# define the spatial extent for subsetting
nyc_extent = (285.5, 286.5, 40, 41.5)

#archive = Herbie("2023-01-01", model='gefs', save_dir='data/herbie', member=0)

# archive.PRODUCTS:
# {'atmos.5': 'Half degree atmos PRIMARY fields (pgrb2ap5); ~83 most common variables.',
#  'atmos.5b': 'Half degree atmos SECONDARY fields (pgrb2bp5); ~500 least common variables',
#  'atmos.25': 'Quarter degree atmos PRIMARY fields (pgrb2sp25); ~35 most common variables',
#  'wave': 'Global wave products.',
#  'chem.5': 'Chemistry fields on 0.5 degree grid',
#  'chem.25': 'Chemistry fields on 0.25 degree grid'}

# # figure out the search strings
# inv1 = archive.inventory()
# inv3.sort_values(by=['variable']).iloc[:, 6:]
# with pd.option_context('display.max_rows', None):
#         print(inv3.sort_values(by=['variable']).iloc[:, 6:])
# inv4.loc[inv4['variable'] == 'SHTFL', :].iloc[:, 6:]

# # grab example files to get info about the parameters
# import cfgrib
# archive = Herbie("2023-01-01", model='gefs', product='atmos.25',
#                  save_dir='data/herbie', member=0)
# search_string = '|'.join(rasp_params['search'][rasp_params['product'] == 'atmos.25'])
# archive.download(search_string)
# inv1 = archive.inventory(search_string)
# with pd.option_context('display.max_rows', None):
#         print(inv1.sort_values(by=['variable']).iloc[:, 6:])

# try FastHerbie
runs = pd.date_range(start="2022-03-01 12:00", periods=2, freq='D')
fxx = range(0, 12, 3)
members = range(0, 3)

# sl_archive = FastRegionHerbie(members, DATES=runs, model='gefs',
#                               product='atmos.25', fxx=fxx, save_dir='data/herbie',
#                               member=0, extent=nyc_extent, extent_name='nyc')

ens_archive = EnsHerbie(members, DATES=runs, model='gefs', product='atmos.25',
                        fxx=fxx, save_dir='data/herbie', extent=nyc_extent,
                        extent_name='nyc')

search_string = '|'.join(rasp_params['search'][rasp_params['product'] == 'atmos.25'])
ens_archive.download(search_string)

# I think the final outcome will be that we loop through and download each
# model/product combo

# would be nice if we could get a progress bar here



# second step: modify as needed

# have to get the 0.5 and 0.25 resolution variables separately
params_0p5 = rasp_params.loc[np.isin(rasp_params['product'], ['atmos.5', 'atmos.5b']), :]
gefs_0p5 = get_nwp_dataset(params_0p5, 'nwp_data/gefs/20230101', [0])

params_0p25 = rasp_params.loc[rasp_params['product'] == 'atmos.25', :]
gefs_0p25 = get_nwp_dataset(params_0p25, 'nwp_data/gefs/20230101', [0])

# third step: combine with rechunker
