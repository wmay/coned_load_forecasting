'''Download and organize NWP forecasts.

This script relies on the wgrib2 system package, which can be installed with
spack or conda. (I recommend spack.)
'''

# pip install 'herbie-data[extras]'
# pip install 'dask[distributed]'
# pip install graphviz
import glob
import numpy as np
import pandas as pd
import xarray as xr
# from herbieplus.xarray import open_herbie_dataset
from herbieplus.collection import NwpCollection
# from rechunker import rechunk
from dask.distributed import Client

# following https://examples.dask.org/applications/embarrassingly-parallel.html#Start-Dask-Client-for-Dashboard
client = Client()
# client = Client(threads_per_worker=4, n_workers=1)
client
client.dashboard_link

# first step: download files with herbie, herbieplus

# GEFS

# parameters from Rasp and Lerch 2018
rasp_params = pd.read_csv('config/gefs.csv')
rasp_params = rasp_params.loc[~pd.isnull(rasp_params['product']), :]
rasp_params = rasp_params.loc[rasp_params['fileType'] == 'forecast', :]
# skip atmos.5b for now because there's no average file
rasp_params = rasp_params.loc[rasp_params['product'] != 'atmos.5b', :]
# define the spatial extent for subsetting
nyc_extent = (285.5, 286.5, 40, 41.5)
# April through Sept., starting 2021
runs = pd.date_range(start=f"{2021}-04-01 12:00", periods=183, freq='D')
for y in range(2022, 2025):
    y_runs = pd.date_range(start=f"{y}-04-01 12:00", periods=183, freq='D')
    runs = runs.union(y_runs)
# analysis file contains different variables, so it must be downloaded
# separately
# fxx = range(0, 24 * 8, 3)
fxx = range(3, 24 * 8, 3)
# members = range(0, 31)
members = ['avg']



# # would be nice if we could get a progress bar here
# for y in range(2021, 2025):
#     print(f'--- Starting year {y}')
#     y_runs = pd.date_range(start=f"{y}-04-01 12:00", periods=183, freq='D')
#     for product in list(rasp_params['product'].unique()):
#         print(f'-- Starting product {product}')
#         product_params = rasp_params.loc[rasp_params['product'] == product, :]
#         product_search = '|'.join(product_params['search'])
#         downloader = HerbieCollection(y_runs, 'gefs', product, product_search, fxx,
#                                       members=members, save_dir='data/herbie',
#                                       extent=nyc_extent, extent_name='nyc')
#         print(f'-- search hash {downloader.search_hash}')
#         downloader.download(threads=5)


# testing new NwpCollection
y_runs = pd.date_range(start=f"{2021}-04-01 12:00", periods=5, freq='D')
product_params = rasp_params.loc[rasp_params['product'] == 'atmos.5', :]
product_search = '|'.join(product_params['search'])
fxx = range(3, 24, 3)

downloader = NwpCollection(y_runs, 'gefs', 'atmos.5', product_search, fxx,
                           members=members, save_dir='nwp_test',
                           extent=nyc_extent)

downloader.get_status()

downloader.collection_size()

x = downloader.download()




# how big is the whole thing?

runs = pd.date_range(start=f"{2021}-04-01 12:00", periods=183, freq='D')
for y in range(2022, 2025):
    y_runs = pd.date_range(start=f"{y}-04-01 12:00", periods=183, freq='D')
    runs = runs.union(y_runs)
fxx = range(3, 24 * 8, 3)
members = ['avg']

downloader2 = NwpCollection(y_runs, 'gefs', 'atmos.5', product_search, fxx,
                            members=members, save_dir='nwp_test',
                            extent=nyc_extent)

downloader2.collection_size()
downloader2.file_size

product_params = rasp_params.loc[rasp_params['product'] == 'atmos.25', :]
product_search = '|'.join(product_params['search'])

downloader3 = NwpCollection(y_runs, 'gefs', 'atmos.25', product_search, fxx,
                            members=members, save_dir='nwp_test',
                            extent=nyc_extent)

downloader3.collection_size()
downloader2.file_size




# can use this to remove troublesome incomplete files
# import os, glob
# from herbieplus.xarray import IncompleteDataException
# date_dirs = glob.glob('data/herbie/gefs/*')
# date_dirs.sort()
# params_atmosp5 = rasp_params.loc[rasp_params['product'] == 'atmos.5', :]

# runs = []
# for y in range(2021, 2025):
#         y_runs = pd.date_range(start=f"{y}-04-01 12:00", periods=183, freq='D')
#         runs.extend(list(y_runs))

# for r in runs:
#     success = False
#     while not success:
#         try:
#             gefs_0p25 = open_herbie_dataset(params_0p25, 'data/herbie/gefs',
#                                             ['avg'], runs=[r])
#             success = True
#         except IncompleteDataException as e:
#             inc_file = str(e).split(':')[0]
#             print(f'Removing {inc_file}')
#             os.remove(inc_file)




# second step: modify as needed

# try xarray parallel reading

from herbieplus.xarray import get_nwp_paths

nwp_files = get_nwp_paths('nwp_test/gefs', 'pgrb2a', ['avg'], runs=downloader.DATES)

filters = {'typeOfLevel': 'isobaricInhPa', 'level': 500}
backend_kwargs = {'filter_by_keys': filters, 'indexpath': ''}

nwp_m = xr.open_mfdataset(nwp_files, concat_dim=['time', 'number', 'step'],
                          engine='cfgrib', decode_timedelta=True,
                          combine='nested', parallel=True,
                          backend_kwargs=backend_kwargs)



# because there are so many files, it makes more sense to aggregate them into
# daily files first
runs = []
for y in range(2021, 2025):
    y_runs = pd.date_range(start=f"{y}-04-01 12:00", periods=183, freq='D')
    runs.extend(list(y_runs))

params_0p5 = rasp_params.loc[np.isin(rasp_params['product'], ['atmos.5', 'atmos.5b']), :]
for r in runs:
    print(r)
    r_ds = open_herbie_dataset(params_0p5, 'data/herbie/gefs', ['avg'],
                               runs=[r])
    r_nc = f'results/get_nwp_data/daily_netcdf/gefs_0p5/{r:%Y%m%d}.nc'
    r_ds.to_netcdf(r_nc)
    r_ds.close()
gefs_0p5_files = glob.glob('results/get_nwp_data/daily_netcdf/gefs_0p5/*.nc')
gefs_0p5_files.sort()
gefs_0p5 = xr.open_mfdataset(gefs_0p5_files, decode_timedelta=True)
gefs_0p5.to_netcdf('results/get_nwp_data/gefs_0p5.nc')
gefs_0p5.close()

params_0p25 = rasp_params.loc[rasp_params['product'] == 'atmos.25', :]
for r in runs[136:]:
    print(r)
    r_ds = open_herbie_dataset(params_0p25, 'data/herbie/gefs', ['avg'],
                               runs=[r])
    r_nc = f'results/get_nwp_data/daily_netcdf/gefs_0p25/{r:%Y%m%d}.nc'
    r_ds.to_netcdf(r_nc)
    r_ds.close()
gefs_0p25_files = glob.glob('results/get_nwp_data/daily_netcdf/gefs_0p25/*.nc')
gefs_0p25_files.sort()
gefs_0p25 = xr.open_mfdataset(gefs_0p25_files, decode_timedelta=True)
gefs_0p25.to_netcdf('results/get_nwp_data/gefs_0p25.nc')
gefs_0p25.close()

# have to get the 0.5 and 0.25 resolution variables separately
# params_0p5 = rasp_params.loc[np.isin(rasp_params['product'], ['atmos.5', 'atmos.5b']), :]
# gefs_0p5 = open_herbie_dataset(params_0p5, 'data/herbie/gefs', ['avg'])

# params_0p25 = rasp_params.loc[rasp_params['product'] == 'atmos.25', :]
# gefs_0p25 = open_herbie_dataset(params_0p25, 'data/herbie/gefs', ['avg'])

# import cfgrib

# x1 = cfgrib.open_datasets('data/herbie/gefs/20210401/nyc_subset_fc1be7eb__geavg.t12z.pgrb2s.0p25.f033')

# x2 = cfgrib.open_datasets('data/herbie/gefs/20210401/nyc_subset_fc86e7eb__geavg.t12z.pgrb2s.0p25.f021')


# third step: combine with rechunker -- not necessary because the data is so
# small!

# gefs_0p5_rechunker = rechunk(
#     gefs_0p5,
#     {'step': 16, "time": 30, "latitude": 4, "longitude": 3},
#     "1GB",
#     "results/get_nwp_data/gefs_0p5.zarr"
# )
# gefs_0p5_rechunker.execute()
