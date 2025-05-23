'''Prepare NWP data for model training
'''

import os
os.chdir('..')
# https://docs.xarray.dev/en/stable/user-guide/dask.html#reading-and-writing-data
os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'
import numpy as np
import pandas as pd
import xarray as xr
from nwpdownload import NwpCollection
from nwpdownload.nwppath import NwpPath
from dask.distributed import Client

# NYC lat/lon boundaries
nyc_extent = (285.5, 286.5, 40, 41.5)
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



# # how long does it take to read a single day's data?

# c1 = NwpCollection(runs_12utc[:5], model='gefs', product='atmos.25',
#                    fxx=gefs_fxx_fct, members=range(0, 31), extent=nyc_extent,
#                    search=':TMP:2 m above ground:', save_dir='data/nwpdownload')

# # get nested list of files for xarray
# nested_files = []
# for i in range(len(c1.DATES)):
#     date_files = []
#     for j in range(len(c1.fxx)):
#         fxx_files = []
#         for k in range(len(c1.members)):
#             fijk = NwpPath(c1.DATES[i], model=c1.model, product=c1.product,
#                            fxx=c1.fxx[j], member=c1.members[k],
#                            save_dir='data/nwpdownload')
#             fxx_files.append(fijk.get_localFilePath())
#         date_files.append(fxx_files)
#     nested_files.append(date_files)

# n1 = xr.open_mfdataset(nested_files, concat_dim=['time', 'step', 'number'],
#                        combine='nested', engine='cfgrib', decode_timedelta=True,
#                        backend_kwargs={'filter_by_keys': {'shortName': '2t'}})
# # 10-ish seconds

# client = Client(n_workers=4, threads_per_worker=1, memory_limit='1GB',
#                 processes=True)
# client.dashboard_link

# n2 = xr.open_mfdataset(nested_files, concat_dim=['time', 'step', 'number'],
#                        combine='nested', engine='cfgrib', decode_timedelta=True,
#                        backend_kwargs={'filter_by_keys': {'shortName': '2t'}},
#                        parallel=True)


# This next section belongs in nwpdownload and script 7. I'll move it when it's
# more organized

import dask, cfgrib, warnings
import dask.array as da

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
        return out

    def _members_arr(self, coords, var_conf):
        return np.stack([ self._array_from_coords(coords + (i, ), var_conf)
                          for i in range(len(self.members)) ])

    # def _fxx_arr(self, coords, var_conf):
    #     return np.stack([ self._members_arr(coords + (i, ), var_conf)
    #                       for i in range(len(self.fxx)) ])

    def _fxx_arr_map(self, var_conf, block_id=None, block_info=None):
        # print(block_id)
        return np.stack([ self._members_arr((block_id[0], i), var_conf)
                          for i in range(len(self.fxx)) ])

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

    def delayed_ds(self):
        # read a single grib2 file to get coordinates and attributes
        # Getting info about the dataset--
        # - array shape (x/y coordinates)
        # - variable info for filter_by_keys
        # - attributes
        f0 = NwpPath(self.DATES[0], model=self.model, product=self.product, fxx=self.fxx[0],
                     member=self.members[0], save_dir='data/nwpdownload')
        grib_path = f0.get_localFilePath()
        ds_list = cfgrib.open_datasets(grib_path, decode_timedelta=True)
        attr_dict = {}
        conf_list = []
        for ds in ds_list:
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
        arrs = []
        for conf in conf_list:
            arr = self._delayed_collection_arr(conf)
            arr.attrs = attr_dict[arr.name]
            # print(arr)
            arrs.append(arr)
        return xr.merge(arrs, combine_attrs='drop_conflicts')


from dask_kubernetes.operator import KubeCluster, make_cluster_spec
# from dask.distributed import LocalCluster
# worker_resources = {'limits': {'cpu': '1', 'memory': worker_memory}}

worker_resources = {'limits': {'cpu': '1'}} # because this is disk-intensive
scheduler_address = f'tcp://wmay-dask-scheduler:8786'
env = {'HDF5_USE_FILE_LOCKING': 'FALSE',
       'DASK_SCHEDULER_ADDRESS': scheduler_address,
       'EXTRA_PIP_PACKAGES': 'netCDF4 git+https://github.com/ASRCsoft/nwpdownload'}

spec = make_cluster_spec(name='wmay-dask', n_workers=2,
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


r1 = CollectionReader(runs_12utc[:30], model='gefs', product='atmos.25',
                      fxx=gefs_fxx_fct, members=range(0, 31), extent=nyc_extent,
                      search=':TMP:2 m above ground:', save_dir='/mnt/nwpdownload')
ds1 = r1.delayed_ds()

# # neat, let's try writing a netcdf file
# ds1.to_netcdf('results/process_nwp_data/gefs_tv.nc')

# https://stackoverflow.com/a/40818232
comp = { 'zlib': True, 'complevel': 4 }
encoding = { var: comp for var in ds1.data_vars }
# When using the dask cluster, there seem to be contradictory expectations about
# where the netcdf is written. The only way this works is to wrap the whole
# function in `delayed`, so that every part of it runs on the remote cluster
delayed_obj = dask.delayed(ds1.to_netcdf)('/mnt/nwpdownload/gefs_tv2.nc',
                                          encoding=encoding)
delayed_obj.compute()


# ds2 = xr.open_dataset('results/process_nwp_data/gefs_tv.nc')
# ds2.close() # looks good


# local version for debugging
# client = Client(n_workers=2, threads_per_worker=1, memory_limit='1GB',
#                 processes=True)
# client.dashboard_link

# r1 = CollectionReader(runs_12utc[:2], model='gefs', product='atmos.25',
#                       fxx=gefs_fxx_fct, members=range(0, 31), extent=nyc_extent,
#                       search=':TMP:2 m above ground:', save_dir='data/nwpdownload')
# ds1 = r1.delayed_ds()

# delayed_obj = dask.delayed(ds1.to_netcdf)('gefs_tv2.nc',
#                                           encoding=encoding)
# delayed_obj.compute()



# now for deriving the predictor variables

def as_fahrenheit(x):
    '''Convert Celsius to Fahrenheit.
    '''
    return x * 9/5 + 32

def coned_wet_bulb(db, dp):
    '''Derive wet bulb temperature from dry bulb temperature and dew point,
    using ConEd's method. Inputs are assumed to be in Kelvin.
    '''
    pr = 1013.25
    dbc = db - 273.15
    dpc = dp - 273.15
    ddepc = (dbc - dpc) * 0.18
    wbc = dbc - (0.035 * ddepc - 0.00072 * ddepc * (ddepc - 1)) *\
        (dbc + dpc + 95.556 - (pr / 30.474))
    return as_fahrenheit(wbc)

def coned_eff_temp(db, dp):
    dbc = db - 273.15
    wbf = coned_wet_bulb(db, dp)
    return (as_fahrenheit(dbc) + wbf) / 2



ds2 = xr.open_dataset('results/process_nwp_data/gefs_tv.nc', mask_and_scale=True)

# Somehow we got temperature values above 29000, which I hope is a misreading of
# missing data.
(ds2['t2m'] > 29000).sum()

ds2['t2m'].values = np.where(ds2['t2m'] > 500, np.nan, ds2['t2m'].values)
ds2['d2m'].values = np.where(ds2['d2m'] > 500, np.nan, ds2['d2m'].values)

# get the hour in EST
est_hour = ((12 - 5) + ds2['step']) % 24
ds2 = ds2.assign_coords({'est_hour': est_hour})
est_day = ((12 - 5) + ds2['step']) // 24
ds2 = ds2.assign_coords({'est_day': est_day})

# use only 9am to 9pm
ds2 = ds2.where((ds2['est_hour'] >= 9) & (ds2['est_hour'] <= 21), drop=True)

ds2['eff_temp'] = coned_eff_temp(ds2['t2m'], ds2['d2m'])
eff_temps = ds2['eff_temp'].groupby('est_day').max()

np.isnan(eff_temps).sum() / eff_temps.size
# somehow 60% of this is missing, omfg


# for other variables, define days as starting 9pm (because that's the latest we
# include for effective temperature) and just get the mean by day
