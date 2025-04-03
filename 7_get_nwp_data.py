#!/<path-to-venv>/bin/python

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

#archive = Herbie("2023-01-01", model='gefs', save_dir='nwp_data', member=0)

# sl_archive = SlimHerbie("2023-01-01", model='gefs', save_dir='nwp_data',
#                         member=0, extent=nyc_extent, extent_name='nyc')

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

# concatenate the search strings
# search_string = '|'.join(rasp_params['search'][rasp_params['product'] == 'atmos.25'])

# # grab example files to get info about the parameters
# import cfgrib
# archive = Herbie("2023-01-01", model='gefs', product='atmos.25',
#                  save_dir='nwp_data', member=0)
# search_string = '|'.join(rasp_params['search'][rasp_params['product'] == 'atmos.25'])
# archive.download(search_string)
# inv1 = archive.inventory(search_string)
# with pd.option_context('display.max_rows', None):
#         print(inv1.sort_values(by=['variable']).iloc[:, 6:])
# nwp1 = cfgrib.open_datasets(str(archive.get_localFilePath(search_string)))

# archive = Herbie("2023-01-01", model='gefs', product='atmos.25',
#                  save_dir='nwp_data', member=0)
# f1 = download_extent_only(archive, search_string, nyc_extent, name='nyc')

# try FastHerbie
runs = pd.date_range(start="2022-03-01 12:00", periods=2, freq='D')
fxx = range(0, 12, 3)
members = range(0, 3)

# archive = FastHerbie(runs, model='gefs', product='atmos.25', fxx=fxx,
#                      save_dir='nwp_data', member=0)

# sl_archive = FastRegionHerbie(members, DATES=runs, model='gefs',
#                               product='atmos.25', fxx=fxx, save_dir='nwp_data',
#                               member=0, extent=nyc_extent, extent_name='nyc')

ens_archive = EnsHerbie(members, DATES=runs, model='gefs', product='atmos.25',
                        fxx=fxx, save_dir='nwp_data', extent=nyc_extent,
                        extent_name='nyc')

search_string = '|'.join(rasp_params['search'][rasp_params['product'] == 'atmos.25'])
# sl_archive.download(search_string)

ens_archive.download(search_string)

# will need to loop through each product
for f in archive.objects:
    print('getting', f.__repr__())
    download_extent_only(f, search_string, nyc_extent, name='nyc')

for f in sl_archive.objects:
    print('getting', f.__repr__())
    print(f.extent)
    f.download(search_string)

# would be nice if we could get a progress bar here






# second step: modify as needed

def append_level_to_varnames(ds):
    '''Add the vertical coordinate to the variable names.
    '''
    z_coord = list(ds.coords.keys())[3]
    z_value = float(ds[z_coord].values)
    for v in list(ds.data_vars):
        ds = ds.rename({v: v + '_' + str(z_value)})
    return ds

def remove_z_coordinate(ds):
    '''Remove the vertical coordinate from a dataset, and add it as an attribute
    instead.
    '''
    # cfgrib stores the vertical coordinate as the fourth value
    z_coord = list(ds.coords.keys())[3]
    z_value = float(ds[z_coord].values)
    vars = list(ds.data_vars)
    for v in vars:
        ds[v].attrs['typeOfLevel'] = z_coord
        ds[v].attrs[z_coord] = z_value
    return ds.drop_vars(z_coord)

def merge_nwp_variables(ds_list):
    '''Standardize variable coordinates and merge them into a single dataset.
    '''
    out_list = []
    for ds in ds_list:
        # print(ds.data_vars.keys())
        ds = append_level_to_varnames(ds)
        out_list.append(remove_z_coordinate(ds))
    return xr.merge(out_list)


# f1 = 'nwp_data/gefs/20220301/nyc_subset_32ef462d__gec00.t06z.pgrb2s.0p25.f000'
# nwp1 = cfgrib.open_datasets(f1)
# # works but names don't match NCEP docs
# nwp1 = [ remove_z_coordinate(ds) for ds in nwp1 ]
# nwp1 = xr.merge(nwp1)
# seems like I'll have to open each type with open_mfdataset, then merge them


# try the same with grib2io -- maybe not needed here

# can't get paths from Herbie because the subsetting changed the names :(
# grib_paths = [ obj.get_localFilePath() for obj in archive.objects ]

def get_nwp_paths(dir, product, members):
    '''Get the paths to the NWP files for a given date and product.
    '''
    out = []
    for m in members:
        m = str(m).zfill(2)
        m_glob = f'{dir}/subset_*__ge*{m}.t00z.{product}.*.f*'
        m_files = glob.glob(m_glob)
        m_files = [ f for f in m_files if f[-4:] != '.idx' ]
        # sort the files by forecast hour
        m_files.sort(key=lambda x: int(x[-3:]))
        out.append(m_files)
    return out

def get_nwp_product_vars(var_df, dir, product, members):
    '''Get a list of variables for a given NWP product.
    '''
    nwp_files = get_nwp_paths(dir, product, members)
    nwp_ds = []
    for row in var_df.drop_duplicates(['typeOfLevel', 'level']).itertuples():
        if pd.isnull(row.typeOfLevel):
            continue
        filters = {'typeOfLevel': row.typeOfLevel}
        if not pd.isnull(row.level):
            filters['level'] = row.level
        backend_kwargs = {'filter_by_keys': filters}
        nwp_m = xr.open_mfdataset(nwp_files, concat_dim=['number', 'step'],
                                  engine='cfgrib', decode_timedelta=True,
                                  combine='nested',
                                  backend_kwargs=backend_kwargs)
        nwp_ds.append(nwp_m)
    return nwp_ds

def get_nwp_dataset(var_df, dir, members):
    '''Get the NWP dataset for a given date and ensemble members.
    '''
    var_df = var_df.loc[~pd.isnull(var_df['product']), :]
    nwp_ds = []
    for product, product_vars in var_df.groupby('product'):
        product_txt = {'atmos.5': 'pgrb2a', 'atmos.25': 'pgrb2s', 'atmos.5b':
                       'pgrb2b'}[product]
        product_vars = get_nwp_product_vars(product_vars, dir, product_txt,
                                            members)
        nwp_ds.extend(product_vars)
    return merge_nwp_variables(nwp_ds)


# have to get the 0.5 and 0.25 resolution variables separately

ds3 = get_nwp_product_vars(rasp_params.loc[rasp_params['product'] == 'atmos.5b', :],
                           'nwp_data/gefs/20230101', 'pgrb2b', [0])


# backend_kwargs={'filter_by_keys': {'typeOfLevel': 'surface'}}
# backend_kwargs={'filter_by_keys': {'typeOfLevel': 'heightAboveGround', 'level': 2}}

# nwp_files = get_nwp_paths('nwp_data/gefs/20220301', range(0, 2))

# ds2 = xr.open_mfdataset(nwp_files, concat_dim=['number', 'step'],
#                         engine='cfgrib', decode_timedelta=True,
#                         combine='nested', backend_kwargs=backend_kwargs)

params_0p5 = rasp_params.loc[np.isin(rasp_params['product'], ['atmos.5', 'atmos.5b']), :]
gefs_0p5 = get_nwp_dataset(params_0p5, 'nwp_data/gefs/20230101', [0])

params_0p25 = rasp_params.loc[rasp_params['product'] == 'atmos.25', :]
gefs_0p25 = get_nwp_dataset(params_0p25, 'nwp_data/gefs/20230101', [0])

# third step: combine with rechunker
