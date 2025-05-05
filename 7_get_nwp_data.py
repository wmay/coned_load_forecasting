'''Download and organize NWP forecasts.
'''

import glob
import numpy as np
import pandas as pd
import xarray as xr
from nwpdownload import NwpCollection
from nwpdownload.kubernetes import k8s_download_cluster
from dask.distributed import Client
# from rechunker import rechunk

# In general, the approach here is to get the forecast variables going out to 8
# days ahead.

# However, because we're predicting based on 3 days of weather, for the shorter
# forecasts we need a way to fill in the weather that's already happened. It
# gets complicated, because the analysis files (f00) don't contain all of the
# predictor variables. So the best we can do is to get all of the 00, 03, and 06
# files to fill in the missing data with the shortest-term, and presumably most
# accurate, forecast.

# GEFS adds another twist, because I need individual ensemble members to
# calculate TV statistics. So for TV variables get the ensemble members, while
# for the rest get only the means and spreads. And it's also split into
# different "product" files. They should probably make it even more complicated
# IMO, otherwise a few people may actually use the data.

# Because there are so many sets of data to collect, I think it's best to
# collect the details in **kwargs dictionaries, then iterate over them.

# (For now just get the 8 days of forecasts, to get started.)

cluster = k8s_download_cluster('wmay-download',
                               '/rdma/hulk/coe/Will/nwp_download',
                               n_workers=100, threads_per_worker=6,
                               port_forward_cluster_ip=True)
client = Client(cluster)

# common parameters
nyc_extent = (285.5, 286.5, 40, 41.5)
runs_daily_12utc = pd.date_range(start=f"{2021}-04-01 12:00", periods=183, freq='D')
for y in range(2022, 2025):
    y_runs = pd.date_range(start=f"{y}-04-01 12:00", periods=183, freq='D')
    runs_daily_12utc = runs_daily_12utc.union(y_runs)

# For variables, I think the cleanest option would be to list the desired
# variables, then check an inventory to decide which product to fetch them from.
# Rather than manually record the product for each variable


# GEFS

# Files split into atmos.25, atmos.5a, and atmos.5b. atmos.5b only available as
# individual ensemble members.

# TV variables (temperature, humidity) must be downloaded as individual ensemble
# members for calculating TV summary statistics.

# f00 is missing some variables. So to fill in data, get what we can from f00.
# Get everything from f03. And get what's missing from f00 from f06 to fill that
# in as well. Phwew! Since we're getting everything for f03 anyway, that should
# be simple enough. For f06, need to get 12utc and other runs separately,
# because they won't have the same variables.

# I'll just get 12utc f03, f06 stuff normally. Then have a separate 00/06/18utc
# set for the historical backfilling. No, this is getting ridiculous. Get 12utc
# f06+, then f00 and f03 all runs, then only f06 00/08/18. That requires fewest
# collections.

# Instead of creating separate collections, why not use same search for f00, and
# what's missing will be excluded automatically?

# So the set of downloads are--

# regular forecasts-- f00+ (12utc, each product, mean+spread -- 3 collections)
# TV forecasts-- f00+ (12utc, atmos.25, ensemble members -- 1 collection)
# backfill-- f00, f03 (00/08/18utc, regular forecasts -- 3 collections)
# backfill-- f00, f03 (00/08/18utc, TV ensemble members -- 1 collection)
# backfill-- f06 (00/08/18utc, missing f00 mean+spread -- 1? collection)

# need to generate the 3 collections from the list of desired variables

# would be nice to generate the backfill collections as well

gefs_fxx_fct = range(3, 24 * 8, 3) # out to 8 days ahead

gefs_tv_fct = {
    'model': 'gefs',
    'product': 'atmos.25',
    'DATES': runs_daily_12utc,
    'fxx': gefs_fxx_fct,
    'members': range(0, 31),
    'extent': nyc_extent,
    'search': '|'.join([
        ':TMP:2 m above ground:',
        ':DPT:2 m above ground:'
    ])
}

gefs_0p25_fct = {
    'model': 'gefs',
    'product': 'atmos.25',
    'DATES': runs_daily_12utc,
    'fxx': gefs_fxx_fct,
    'members': ['avg', 'spread'],
    'extent': nyc_extent,
    'search': '|'.join([
        ':TMP:2 m above ground:',
        ':DPT:2 m above ground:'
    ])
}

gefs_0p5a_fct = {
    'model': 'gefs',
    'product': 'atmos.25',
    'DATES': runs_daily_12utc,
    'fxx': gefs_fxx_fct,
    'members': ['avg', 'spread'],
    'extent': nyc_extent,
    'search': '|'.join([
        ':TMP:2 m above ground:',
        ':DPT:2 m above ground:'
    ])
}

gefs_0p5b_fct = {
    'model': 'gefs',
    'product': 'atmos.25',
    'DATES': runs_daily_12utc,
    'fxx': gefs_fxx_fct,
    'members': range(0, 31),
    'extent': nyc_extent,
    'search': '|'.join([
        ':TMP:2 m above ground:',
        ':DPT:2 m above ground:'
    ])
}

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

downloader = NwpCollection(y_runs, 'gefs', 'atmos.5', product_search, fxx,
                           members=members, save_dir='nwp_test',
                           extent=nyc_extent)



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
