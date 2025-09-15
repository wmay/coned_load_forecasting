'''Download and organize NWP forecasts.
'''

import os, dask, time, warnings
# os.chdir('..')
# https://docs.xarray.dev/en/stable/user-guide/dask.html#reading-and-writing-data
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
from datetime import datetime, timedelta
import numpy as np
import pandas as pd
import xarray as xr
from herbie import Herbie, HerbieLatest, HerbieWait
from nwpdownload import NwpCollection
from nwpdownload.xarray import merge_nwp_variables
from dask.distributed import Client

# last 4 days
start_date = datetime.today() - timedelta(days=4)
all_runs = pd.date_range(start=start_date.strftime('%Y-%m-%d 00:00'),
                         periods=4 * 4 + 2, freq='6h').to_series().index

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
# for the rest get only the means and spreads (actually never mind the spreads--
# they're useless because we don't know the covariances, so they can't be
# aggregated over time). And it's also split into different "product" files.
# They should probably make it even more complicated IMO, otherwise a few people
# may actually use the data.

# Because there are so many sets of data to collect, I think it's best to
# collect the details in **kwargs dictionaries, then iterate over them. Also
# let's make some functions to automatically find the product for each variable,
# and set up the NwpCollection arguments.

def find_variable_product(DATES, fxx, model, products, search, member=None):
    '''Find the product containing the given variable (represented by the search
    term).
    '''
    product_herbie = {}
    for p in products:
        product_herbie[p] = Herbie(DATES[0], fxx=fxx[0], model=model, product=p,
                                   member=member)
    var_products = {}
    for s in search:
        var_products[s] = None
        for p in products:
            p_inv = product_herbie[p].inventory(search=s)
            if p_inv.shape[0]:
                var_products[s] = p
                break
    return var_products

def make_collection_args(DATES, fxx, model, products, search):
    '''Given variables and products, create the collection arguments for each
    product.
    '''
    if model == 'gefs':
        member = 0
    else:
        member = None
    var_products = find_variable_product(DATES, fxx, model, products, search,
                                         member=member)
    args_list = []
    for p in products:
        # which variables are in this product?
        p_vars = [ v for v in search if var_products[v] == p ]
        if not len(p_vars):
            continue
        members = None
        if model == 'gefs':
            if p == 'atmos.5b':
                # no "avg" dataset is available for this product, so it must be
                # derived from the members
                members = range(31)
            else:
                members = ['avg']
        args = {
            'model': model,
            'product': p,
            'DATES': DATES,
            'fxx': fxx,
            'members': members,
            'extent': nyc_extent,
            'search': '|'.join(p_vars)
        }
        args_list.append(args)
    return args_list

# Some common variables

# just use a LocalCluster
client = Client(processes=False, n_workers=10)
# client.wait_for_workers(1)

# NYC lat/lon boundaries
nyc_extent = (285.5, 286.5, 40, 41.5)


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

# regular forecasts-- f03+ (12utc, each product, mean -- 3 collections)
# TV forecasts-- f03+ (12utc, atmos.25, ensemble members -- 1 collection)
# regular backfill-- f00, f03 (all runs, regular forecasts -- 3 collections)
# TV backfill-- f00, f03 (all runs, TV ensemble members -- 1 collection)
# extra backfill-- f06 (00/08/18utc, missing f00 mean -- 1? collection)

# When can I get these? For 06Z,
# - 0p25 avg, 7:08
# - 0p25 members, 7:07
# - 0p5 avg, 7:08
# - 0p5b members, 7:07

# parameters from Rasp and Lerch 2018
gefs_params = pd.read_csv('config/gefs.csv')
gefs_params = gefs_params.loc[~pd.isnull(gefs_params['search']), :]
gefs_searches = gefs_params['search'].unique()
runs_06utc = all_runs[np.arange(len(all_runs)) % 4 == 1]
runs_not_06utc = all_runs[np.arange(len(all_runs)) % 4 != 1]
runs_fct = runs_06utc[-1:] # only need the most recent forecasts
# runs_fct = all_runs[-1:] # for testing
gefs_fxx_fct = range(3, 24 * 8, 3) # out to 8 days ahead

# we need to time this carefully, downloading the forecast data incrementally so
# we can quickly get the tail end of the lead times after they're posted
def get_fct_collection(fxx):
    '''Get collections for forecasts at the given fxx range.
    '''
    # f03+
    gefs_fct_args_list = make_collection_args(runs_fct, fxx,
                                              'gefs',
                                              ['atmos.25', 'atmos.5', 'atmos.5b'],
                                              search=gefs_searches)
    # also want individual ensembles to calculate TV statistics
    gefs_tv_fct = {
        'model': 'gefs',
        'product': 'atmos.25',
        'DATES': runs_fct,
        'fxx': fxx,
        'members': range(0, 31),
        'extent': nyc_extent,
        'search': '|'.join([
            ':TMP:2 m above ground:',
            ':DPT:2 m above ground:'
        ])
    }
    return gefs_fct_args_list + [gefs_tv_fct]

# now f00
gefs_f00_args_list = make_collection_args(all_runs, [0],
                                          'gefs',
                                          ['atmos.25', 'atmos.5', 'atmos.5b'],
                                          search=gefs_searches)
gefs_tv_f00 = get_fct_collection([0])[-1].copy() # copy the TV settings
gefs_tv_f00['DATES'] = all_runs

# now f03
gefs_f03_args_list = make_collection_args(runs_not_06utc, [3],
                                          'gefs',
                                          ['atmos.25', 'atmos.5', 'atmos.5b'],
                                          search=gefs_searches)
gefs_tv_f03 = gefs_tv_f00.copy()
gefs_tv_f03['fxx'] = [3]
gefs_tv_f03['DATES'] = runs_not_06utc

# now f06 (a few variables missing from the f00 files)
gefs_f00_products = find_variable_product(runs_not_06utc, [0], 'gefs',
                                          ['atmos.25', 'atmos.5', 'atmos.5b'],
                                          gefs_searches, member=0)
# what's missing from f00?
f00_missing = [ v for v in gefs_searches if gefs_f00_products[v] is None ]
gefs_f06_args_list = make_collection_args(runs_not_06utc, [6],
                                          'gefs',
                                          ['atmos.25', 'atmos.5', 'atmos.5b'],
                                          search=f00_missing)

gefs_all_collections_args = gefs_f00_args_list + [gefs_tv_f00] + \
    gefs_f03_args_list + [gefs_tv_f03] + gefs_f06_args_list

def download_collections(all_collections_args):
    '''Download a list of NwpCollections.
    '''
    collections = []
    for collection_args in all_collections_args:
        collection = NwpCollection(**collection_args, save_dir='scripts/data/nwpdownload')
        collections.append(collection)
    # how much data is all of this?
    collection_sizes = []
    for collection in collections:
        collection_sizes.append(collection.collection_size(humanize=False))
        print(collection.collection_size())
    # how big is the whole thing?
    from humanize import naturalsize
    naturalsize(sum(collection_sizes))
    for i, collection in enumerate(collections):
        print(f'Starting collection {i+1} of {len(collections)}')
        collection.download()

print('Starting analysis/backfill collections ...')
download_collections(gefs_all_collections_args)

# now we incrementally download the fxx of the most recent forecast
print('Starting forecast collections ...')

# upload timings:

# pgrb2s.0p25:
# gec00: 2025-09-05 06:19:43
# gep01: 2025-09-05 06:23:39
# gep30: 2025-09-05 06:23:40
# pgrb2a.0p5 spr: 2025-09-05 06:23:41
# pgrb2b.0p5 gep30: 2025-09-05 06:23:41
# geavg: 2025-09-05 06:24:54
# gespr: 2025-09-05 06:24:55

# failed file:
# geavg.t06z.pgrb2s.0p25.f156.idx: 2025-09-07 06:54:16
# gespr.t06z.pgrb2s.0p25.f156.idx: 2025-09-07 06:54:15
# but error occurred at 2025-09-07 10:53:50,090?

# latest_fxx = gefs_fxx_fct[0]
latest_fxx = 0
while latest_fxx < gefs_fxx_fct[-1]:
    next_fxx = latest_fxx + 3
    # check for the spread file, because it's the last to be written
    H = HerbieLatest(model="gefs", priority=['aws'], product='atmos.25',
                     fxx=next_fxx, member='spr')
    if H.date >= runs_fct[0]:
        latest_fxx = next_fxx
    else:
        break

cur_fxx = -1
all_downloaded = False
while not all_downloaded:
    new_fxx = [ fxx for fxx in gefs_fxx_fct
                if fxx > cur_fxx and fxx <= min(latest_fxx, gefs_fxx_fct[-1]) ]
    print(f'Getting fxx up to {str(new_fxx[-1])} at {str(datetime.today())}')
    fct_collections = get_fct_collection(new_fxx)
    # keep trying for a maximum of 10 minutes if needed
    max_tries = 20
    ntry = 1
    while ntry <= max_tries:
        try:
            download_collections(fct_collections)
            break
        except:
            time.sleep(30)
        ntry += 1
    if ntry > max_tries:
        # something went wrong
        warnings.warn(f'Failed to download fxx {str(new_fxx[-1])}')
    cur_fxx = new_fxx[-1]
    if cur_fxx == gefs_fxx_fct[-1]:
        all_downloaded = True
    else:
        # wait for the next fxx
        next_fxx = latest_fxx + 3
        HerbieWait(run=runs_fct[0], model='gefs', priority=['aws'],
                   product='atmos.25', member='spr', fxx=next_fxx,
                   wait_for='1h', check_interval='60s')
        latest_fxx = next_fxx


# second step: combine each collection into a single data file

# make sure we have the full collections
gefs_all_collections_args = gefs_all_collections_args + \
    get_fct_collection(gefs_fxx_fct)

def get_file_name(c):
    '''Generate a filename based on the collection attributes.
    '''
    # is this a forecast file, analysis, or fill-in data?
    if max(c.fxx) > 6:
        fct = 'fct'
    elif max(c.fxx) > 0:
        fill_n = max(c.fxx)
        fct = f'fill{fill_n}'
    else:
        fct = 'anl'
    # which product?
    product = {'atmos.25': 'atmos0p25', 'atmos.5': 'atmos0p5',
               'atmos.5b': 'atmos0p5b'}[c.product]
    # individual members or means?
    members = len(c.members) > 2
    out = f'gefs_{product}_{fct}'
    if members:
        out = out + '_members'
    return out + '.nc'


# https://stackoverflow.com/a/40818232
comp = { 'zlib': True, 'complevel': 4 }
# encoding = { var: comp for var in ds1.data_vars }

gefs_collections = []
for collection_args in gefs_all_collections_args:
    collection = NwpCollection(**collection_args, save_dir='scripts/data/nwpdownload')
    gefs_collections.append(collection)

# write each collection to a netcdf file
for c in gefs_collections:
    out_path = 'scripts/data/nwpdownload/' + get_file_name(c)
    print(f'writing to {out_path}')
    ds_list = c.open_datasets()
    ds_c = merge_nwp_variables(ds_list)
    encoding = { var: comp for var in ds_c.data_vars }
    # When using the dask cluster, there seem to be contradictory expectations
    # about where the netcdf is written. The only way this works is to wrap the
    # whole function in `delayed`, so that every part of it runs on the remote
    # cluster
    delayed_c = dask.delayed(ds_c.to_netcdf)(out_path, encoding=encoding)
    delayed_c.compute()

client.close()
