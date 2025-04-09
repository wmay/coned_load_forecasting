'''Download and organize NWP forecasts.

This script relies on the wgrib2 system package, which can be installed with
spack or conda. (I recommend spack.)
'''

# pip install 'herbie-data[extras]'
#import os, glob
import numpy as np
import pandas as pd
import xarray as xr
from herbieplus.dataset import HerbieCollection
from herbieplus.xarray import open_herbie_dataset
#from rechunker import rechunk

# first step: download files with Herbie, subset with wgrib2, organize

# GEFS test run

# https://herbie.readthedocs.io/en/latest/gallery/noaa_models/gefs.html

# another complication for GEFS is that there are both .5 and .25 degree
# products

# archive = Herbie("2023-01-01", model='gefs', product='atmos.25',
#                  save_dir='data/herbie', fxx=3, member='avg')

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
#         print(inv1.sort_values(by=['variable']).iloc[:, 6:])
# inv4.loc[inv4['variable'] == 'SHTFL', :].iloc[:, 6:]

# parameters from Rasp and Lerch 2018
rasp_params = pd.read_csv('config/gefs.csv')
rasp_params = rasp_params.loc[rasp_params['fileType'] == 'forecast', :]
# define the spatial extent for subsetting
nyc_extent = (285.5, 286.5, 40, 41.5)
# fxx = range(0, 24 * 8, 3)
# members = range(0, 31)
fxx = range(3, 24 * 8, 3)
members = ['avg']

# would be nice if we could get a progress bar here
for y in range(2021, 2025):
    print(f'--- Starting year {y}')
    y_runs = pd.date_range(start=f"{y}-04-01 12:00", periods=183, freq='D')
    for product in list(rasp_params['product'].unique()):
        print(f'-- Starting product {product}')
        product_params = rasp_params.loc[rasp_params['product'] == product, :]
        product_search = '|'.join(product_params['search'])
        downloader = HerbieCollection(y_runs, 'gefs', product, product_search, fxx,
                                      members=members, save_dir='data/herbie',
                                      extent=nyc_extent, extent_name='nyc')
        print(f'-- search hash {downloader.search_hash}')
        downloader.download(threads=5)


# second step: modify as needed

# have to get the 0.5 and 0.25 resolution variables separately
params_0p5 = rasp_params.loc[np.isin(rasp_params['product'], ['atmos.5', 'atmos.5b']), :]
gefs_0p5 = get_nwp_dataset(params_0p5, 'nwp_data/gefs/20230101', [0])

params_0p25 = rasp_params.loc[rasp_params['product'] == 'atmos.25', :]
gefs_0p25 = get_nwp_dataset(params_0p25, 'nwp_data/gefs/20230101', [0])

# third step: combine with rechunker
