'''Download and organize NWP forecasts.
'''

import os
os.chdir('..')
# https://docs.xarray.dev/en/stable/user-guide/dask.html#reading-and-writing-data
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
import glob
import numpy as np
import pandas as pd
import xarray as xr
from herbie import Herbie
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

# tools for reading the files once downloaded
import dask, cfgrib, warnings
import dask.array as da
from nwpdownload.nwppath import NwpPath
from nwpdownload.xarray import merge_nwp_variables, append_level_to_repeated_varnames

def get_filter_by_keys(arr):
    '''Given an xarray DataArray, get appropriate filter_by_keys.
    '''
    # from each variable, grab shortName, typeOfLevel, and level (if present)
    filter_by_keys = {}
    attrs = arr.attrs
    filter_by_keys['shortName'] = attrs['GRIB_shortName']
    if 'GRIB_typeOfLevel' in attrs.keys():
        filter_by_keys['typeOfLevel'] = attrs['GRIB_typeOfLevel']
        # that's probably good-- no reason to separate data by level
    return filter_by_keys

class CollectionReader(NwpCollection):
    '''Methods to open an NwpCollection in xarray using dask.
    '''

    def _array_from_coords(self, coords, var_conf):
        '''Get grib array from run/fxx/member coordinates.
        '''
        fcoords = NwpPath(self.DATES[coords[0]], model=self.model, product=self.product,
                          fxx=self.fxx[coords[1]], member=self.members[coords[2]],
                          save_dir=self.save_dir)
        grib_path = fcoords.get_localFilePath()
        backend_kwargs = {'filter_by_keys': var_conf['filter_by_keys'],
                          'indexpath': ''}
        try:
            ds = xr.open_dataset(grib_path, engine='cfgrib',
                                 decode_timedelta=True,
                                 backend_kwargs=backend_kwargs)
            out = ds[var_conf['name']].values
        except Exception as e:
            # convert the error to a warning, and return an empty array
            warnings.warn(str(e))
            out = np.full(var_conf['shape'], np.nan)
        # check that the array has the correct shape
        if out.shape != var_conf['shape']:
            warnings.warn('Array has incorrect shape')
            out = np.full(var_conf['shape'], np.nan)
        return out

    def _members_arr(self, coords, var_conf):
        return np.stack([ self._array_from_coords(coords + (i, ), var_conf)
                          for i in range(len(self.members)) ])

    # def _fxx_arr(self, coords, var_conf):
    #     return np.stack([ self._members_arr(coords + (i, ), var_conf)
    #                       for i in range(len(self.fxx)) ])

    def _fxx_arr_map(self, var_conf, block_id=None, block_info=None):
        # print(block_id)
        out = np.stack([ self._members_arr((block_id[0], i), var_conf)
                         for i in range(len(self.fxx)) ])
        return np.expand_dims(out, 0)

    # Rather than create a dask array for each file, create one for each
    # forecast run. This is much more manageable for the dask scheduler.
    def _delayed_collection_arr(self, var_conf):
        coords = {'time': self.DATES, 'step': self.fxx, 'number': self.members}
        # print(var_conf)
        coords.update(var_conf['dims'])
        fxx_shape = (len(self.fxx), len(self.members)) + var_conf['shape']
        # use `map_blocks` instead of `stack`
        n_runs = len(self.DATES)
        arr = da.map_blocks(self._fxx_arr_map, var_conf,
                            dtype=var_conf['dtype'],
                            chunks=((1, ) * n_runs, *fxx_shape),
                            meta=np.array((), dtype=var_conf['dtype']))
        return xr.DataArray(arr, coords=coords, name=var_conf['name'])

    def open_datasets(self):
        # read a single grib2 file to get coordinates and attributes
        # Getting info about the dataset--
        # - array shape (x/y coordinates)
        # - variable info for filter_by_keys
        # - attributes
        # f0 = NwpPath(self.DATES[0], model=self.model, product=self.product, fxx=self.fxx[0],
        #              member=self.members[0], save_dir='data/nwpdownload')
        f0 = NwpPath(self.DATES[0], model=self.model, product=self.product, fxx=self.fxx[0],
                     member=self.members[0], save_dir=self.save_dir)

        # f0 = dask.delayed(NwpPath)(self.DATES[0], model=self.model,
        #                            product=self.product, fxx=self.fxx[0],
        #                            member=self.members[0],
        #                            save_dir=self.save_dir)
        # f0.compute()
        grib_path = f0.get_localFilePath()
        # to make sure this reads from the correct location, must use
        # `dask.delayed`
        # ds_list = cfgrib.open_datasets(grib_path, decode_timedelta=True)
        ds_list = dask.delayed(cfgrib.open_datasets)(grib_path,
                                                     decode_timedelta=True)
        ds_list = ds_list.compute()
        attr_dict = {}
        outer_list = []
        for ds in ds_list:
            conf_list = []
            for v in ds.data_vars:
                conf = {'name': v}
                arr = ds[v]
                conf['filter_by_keys'] = get_filter_by_keys(arr)
                conf['shape'] = arr.shape
                conf['dtype'] = arr.dtype
                dim_names = list(arr.dims)
                conf['dims'] = { dim: arr[dim] for dim in dim_names }
                conf_list.append(conf)
                attr_dict[v] = arr.attrs
            outer_list.append(conf_list)
        out_list = []
        # as long as we get data separately for each vertical level, we can
        # easily merge the arrays. But I haven't guaranteed that yet!
        for conf_list in outer_list:
            arrs = []
            for conf in conf_list:
                arr = self._delayed_collection_arr(conf)
                arr.attrs = attr_dict[arr.name]
                # print(arr)
                arrs.append(arr)
            ds = xr.merge(arrs, combine_attrs='drop_conflicts')
            out_list.append(ds)
        return xr.merge(out_list, combine_attrs='drop_conflicts')


# Some common variables

# speed this up immensely with the help of a cluster running kubernetes
cluster = k8s_download_cluster('wmay-download',
                               '/rdma/hulk/coe/Will/coned_nwp',
                               n_workers=100, threads_per_worker=6,
                               port_forward_cluster_ip=True)
client = Client(cluster)

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

# parameters from Rasp and Lerch 2018
gefs_params = pd.read_csv('config/gefs.csv')
gefs_params = gefs_params.loc[~pd.isnull(gefs_params['search']), :]
gefs_searches = gefs_params['search'].unique()
# April through Sept., starting 2021
year_runs = [ pd.date_range(start=f"{y}-04-01 00:00", periods=183 * 4,
                            freq='6h').to_series() for y in range(2021, 2025) ]
all_runs = pd.concat(year_runs).index
runs_12utc = all_runs[np.arange(len(all_runs)) % 4 == 2]
runs_not_12utc = all_runs[np.arange(len(all_runs)) % 4 != 2]
gefs_fxx_fct = range(3, 24 * 8, 3) # out to 8 days ahead

# f03+
gefs_fct_args_list = make_collection_args(runs_12utc, gefs_fxx_fct,
                                          'gefs',
                                          ['atmos.25', 'atmos.5', 'atmos.5b'],
                                          search=gefs_searches)
# also want individual ensembles to calculate TV statistics
gefs_tv_fct = {
    'model': 'gefs',
    'product': 'atmos.25',
    'DATES': runs_12utc,
    'fxx': gefs_fxx_fct,
    'members': range(0, 31),
    'extent': nyc_extent,
    'search': '|'.join([
        ':TMP:2 m above ground:',
        ':DPT:2 m above ground:'
    ])
}

# now f00
gefs_f00_args_list = make_collection_args(all_runs, [0],
                                          'gefs',
                                          ['atmos.25', 'atmos.5', 'atmos.5b'],
                                          search=gefs_searches)
gefs_tv_f00 = gefs_tv_fct.copy()
gefs_tv_f00['fxx'] = [0]
gefs_tv_f00['DATES'] = all_runs

# now f03
gefs_f03_args_list = make_collection_args(runs_not_12utc, [3],
                                          'gefs',
                                          ['atmos.25', 'atmos.5', 'atmos.5b'],
                                          search=gefs_searches)
gefs_tv_f03 = gefs_tv_f00.copy()
gefs_tv_f03['fxx'] = [3]
gefs_tv_f03['DATES'] = runs_not_12utc

# now f06 (a few variables missing from the f00 files)
gefs_f00_products = find_variable_product(runs_not_12utc, [0], 'gefs',
                                          ['atmos.25', 'atmos.5', 'atmos.5b'],
                                          gefs_searches, member=0)
# what's missing from f00?
f00_missing = [ v for v in gefs_searches if gefs_f00_products[v] is None ]
gefs_f06_args_list = make_collection_args(runs_not_12utc, [6],
                                          'gefs',
                                          ['atmos.25', 'atmos.5', 'atmos.5b'],
                                          search=f00_missing)

gefs_all_collections_args = gefs_fct_args_list + [gefs_tv_fct] + \
    gefs_f00_args_list + [gefs_tv_f00] + gefs_f03_args_list + [gefs_tv_f03] + \
    gefs_f06_args_list

gefs_collections = []
for collection_args in gefs_all_collections_args:
    collection = CollectionReader(**collection_args, save_dir='/mnt/nwpdownload')
    gefs_collections.append(collection)

# how much data is all of this?
collection_sizes = []
for collection in gefs_collections:
    collection_sizes.append(collection.collection_size(humanize=False))
    print(collection.collection_size())

# how big is the whole thing?
from humanize import naturalsize
naturalsize(sum(collection_sizes))

for i, collection in enumerate(gefs_collections):
    if 1 == 3:
        continue
    print(f'Starting collection {i+1} of {len(gefs_collections)}')
    collection.download()



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



# second step: combine each collection into a single data file

# close the previous cluster
client.close()
cluster.close()


from dask_kubernetes.operator import KubeCluster, make_cluster_spec
# from dask.distributed import LocalCluster
# worker_resources = {'limits': {'cpu': '1', 'memory': worker_memory}}

worker_resources = {'limits': {'cpu': '1'}} # because this is disk-intensive
scheduler_address = f'tcp://wmay-dask-scheduler:8786'
env = {'HDF5_USE_FILE_LOCKING': 'FALSE',
       'DASK_SCHEDULER_ADDRESS': scheduler_address,
       # have to use eccodes <2.39 due to an issue reading the regional subset
       # files
       'EXTRA_PIP_PACKAGES': 'eccodes==2.38 typing_extensions netCDF4 git+https://github.com/ASRCsoft/nwpdownload'}

spec = make_cluster_spec(name='wmay-dask', n_workers=10,
                         resources=worker_resources, env=env)
scheduler = spec['spec']['scheduler']['spec']['containers'][0]
worker = spec['spec']['worker']['spec']['containers'][0]
# add the volume for writing data
volume = {'hostPath': {'path': '/rdma/hulk/coe/Will/coned_nwp'},
          'name': 'nwpout'}
mount = {'name': 'nwpout', 'mountPath': '/mnt/nwpdownload'}
spec['spec']['worker']['spec']['volumes'] = [volume]
worker['volumeMounts'] = [mount]

cluster = KubeCluster(custom_cluster_spec=spec, port_forward_cluster_ip=True)
# cluster = KubeCluster(name='wmay-dask', n_workers=10, env=env,
#                       port_forward_cluster_ip=True)
cluster.dashboard_link
cluster.wait_for_workers(10)

client = Client(cluster)

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


# r1 = CollectionReader(runs_12utc, model='gefs', product='atmos.25',
#                       fxx=gefs_fxx_fct, members=range(0, 31), extent=nyc_extent,
#                       search=':TMP:2 m above ground:', save_dir='/mnt/nwpdownload')
# ds1 = r1.delayed_ds()

# # neat, let's try writing a netcdf file
# ds1.to_netcdf('results/process_nwp_data/gefs_tv.nc')

# https://stackoverflow.com/a/40818232
comp = { 'zlib': True, 'complevel': 4 }
# encoding = { var: comp for var in ds1.data_vars }
# # When using the dask cluster, there seem to be contradictory expectations about
# # where the netcdf is written. The only way this works is to wrap the whole
# # function in `delayed`, so that every part of it runs on the remote cluster
# delayed_obj = dask.delayed(ds1.to_netcdf)('/mnt/nwpdownload/gefs_tv2.nc',
#                                           encoding=encoding)
# delayed_obj.compute()

from nwpdownload.xarray import merge_nwp_variables

for c in gefs_collections:
    # skip larger ones while testing
    if len(c.members) > 2:
        continue
    out_path = '/mnt/nwpdownload/' + get_file_name(c)
    print(f'writing to {out_path}')
    # ds_list = c.open_datasets()
    # print(ds_list)
    # ds_c = merge_nwp_variables(ds_list)
    ds_c = c.open_datasets()
    print(ds_c)
    encoding = { var: comp for var in ds_c.data_vars }
    # When using the dask cluster, there seem to be contradictory expectations
    # about where the netcdf is written. The only way this works is to wrap the
    # whole function in `delayed`, so that every part of it runs on the remote
    # cluster
    delayed_c = dask.delayed(ds_c.to_netcdf)(out_path, encoding=encoding)
    delayed_c.compute()

# # neat, let's try writing a netcdf file
# ds1.to_netcdf('results/process_nwp_data/gefs_tv.nc')
    
ds_test = xr.open_dataset('results/get_nwp_data/gefs_atmos0p25_fct.nc')



# second step: modify as needed

# try xarray parallel reading

# from herbieplus.xarray import get_nwp_paths

# nwp_files = get_nwp_paths('nwp_test/gefs', 'pgrb2a', ['avg'], runs=downloader.DATES)

# filters = {'typeOfLevel': 'isobaricInhPa', 'level': 500}
# backend_kwargs = {'filter_by_keys': filters, 'indexpath': ''}

# nwp_m = xr.open_mfdataset(nwp_files, concat_dim=['time', 'number', 'step'],
#                           engine='cfgrib', decode_timedelta=True,
#                           combine='nested', parallel=True,
#                           backend_kwargs=backend_kwargs)



# # because there are so many files, it makes more sense to aggregate them into
# # daily files first
# runs = []
# for y in range(2021, 2025):
#     y_runs = pd.date_range(start=f"{y}-04-01 12:00", periods=183, freq='D')
#     runs.extend(list(y_runs))

# params_0p5 = rasp_params.loc[np.isin(rasp_params['product'], ['atmos.5', 'atmos.5b']), :]
# for r in runs:
#     print(r)
#     r_ds = open_herbie_dataset(params_0p5, 'data/herbie/gefs', ['avg'],
#                                runs=[r])
#     r_nc = f'results/get_nwp_data/daily_netcdf/gefs_0p5/{r:%Y%m%d}.nc'
#     r_ds.to_netcdf(r_nc)
#     r_ds.close()
# gefs_0p5_files = glob.glob('results/get_nwp_data/daily_netcdf/gefs_0p5/*.nc')
# gefs_0p5_files.sort()
# gefs_0p5 = xr.open_mfdataset(gefs_0p5_files, decode_timedelta=True)
# gefs_0p5.to_netcdf('results/get_nwp_data/gefs_0p5.nc')
# gefs_0p5.close()

# params_0p25 = rasp_params.loc[rasp_params['product'] == 'atmos.25', :]
# for r in runs[136:]:
#     print(r)
#     r_ds = open_herbie_dataset(params_0p25, 'data/herbie/gefs', ['avg'],
#                                runs=[r])
#     r_nc = f'results/get_nwp_data/daily_netcdf/gefs_0p25/{r:%Y%m%d}.nc'
#     r_ds.to_netcdf(r_nc)
#     r_ds.close()
# gefs_0p25_files = glob.glob('results/get_nwp_data/daily_netcdf/gefs_0p25/*.nc')
# gefs_0p25_files.sort()
# gefs_0p25 = xr.open_mfdataset(gefs_0p25_files, decode_timedelta=True)
# gefs_0p25.to_netcdf('results/get_nwp_data/gefs_0p25.nc')
# gefs_0p25.close()

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
