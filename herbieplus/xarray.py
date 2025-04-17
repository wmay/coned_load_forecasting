'''Open a collection of Herbie files with xarray.
'''

import glob
import numpy as np
import pandas as pd
import xarray as xr

# Possible extension: add support for grib2io

class IncompleteDataException(Exception):
    pass

def get_nwp_paths(model_dir, product, members, runs=None):
    '''Get the paths to the NWP files for a given date and product.
    '''
    if runs is None:
        subdirs = glob.glob(f'{model_dir}/*')
    else:
        subdirs = [ f'{model_dir}/{d:%Y%m%d}' for d in runs ]
        subdirs.sort()
    out = []
    for d in subdirs:
        d_out = []
        for m in members:
            if type(m) == int:
                m = str(m).zfill(2)
            # m_glob = f'{d}/*subset_*__ge*{m}.t12z.{product}.*.f*'
            m_glob = f'{d}/*ge*{m}.t12z.{product}.*.f*'
            m_files = glob.glob(m_glob)
            m_files = [ f for f in m_files if f[-4:] != '.idx' ]
            # sort the files by forecast hour
            m_files.sort(key=lambda x: int(x[-3:]))
            d_out.append(m_files)
        out.append(d_out)
    return out

def append_level_to_varname(ds, v):
    '''Add the vertical coordinate to a variable name.
    '''
    z_coord_idx = 3 if 'number' in ds.coords.keys() else 2
    z_coord = list(ds.coords.keys())[z_coord_idx]
    if len(ds[z_coord].values.shape):
        print(ds)
        raise Exception(f'Unexpected coordinates for {z_coord}: {ds[z_coord].values}')
    z_value = float(ds[z_coord].values)
    z_units = ds[z_coord].attrs['units']
    return ds.rename({v: f'{v}_{z_value:g}{z_units}'})

def append_level_to_repeated_varnames(ds_list):
    '''Avoid name conflicts by adding the vertical coordinate to repeated
    variable names.
    '''
    varnames = [ list(ds.data_vars.keys()) for ds in ds_list ]
    # flatten
    varnames = [ v for v_list in varnames for v in v_list ]
    varnames, counts = np.unique(varnames, return_counts=True)
    repeated_varnames = varnames[counts > 1]
    if not repeated_varnames.shape[0]:
        return ds_list
    for v in repeated_varnames:
        for i, ds in enumerate(ds_list):
            if v in list(ds.data_vars):
                ds_list[i] = append_level_to_varname(ds, v)
    return ds_list

def remove_z_coordinate(ds):
    '''Remove the vertical coordinate from a dataset, and add it as an attribute
    instead.
    '''
    # cfgrib stores the vertical coordinate as the fourth value
    z_coord_idx = 3 if 'number' in ds.coords.keys() else 2
    z_coord = list(ds.coords.keys())[z_coord_idx]
    if len(ds[z_coord].values.shape):
        print(ds)
        raise Exception(f'Unexpected coordinates for {z_coord}: {ds[z_coord].values}')
    z_value = float(ds[z_coord].values)
    vars = list(ds.data_vars)
    for v in vars:
        ds[v].attrs['typeOfLevel'] = z_coord
        ds[v].attrs[z_coord] = z_value
    return ds.drop_vars(z_coord)

def merge_nwp_variables(ds_list):
    '''Standardize variable coordinates and merge them into a single dataset.
    '''
    out_list = append_level_to_repeated_varnames(ds_list)
    out_list = [ remove_z_coordinate(ds) for ds in out_list ]
    return xr.merge(out_list)

def get_nwp_product_vars(var_df, dir, product, members, runs=None):
    '''Get a list of variables for a given NWP product.
    '''
    nwp_files = get_nwp_paths(dir, product, members, runs=runs)
    nwp_ds = []
    for row in var_df.drop_duplicates(['typeOfLevel', 'level']).itertuples():
        if pd.isnull(row.typeOfLevel):
            continue
        filters = {'typeOfLevel': row.typeOfLevel}
        if not pd.isnull(row.level):
            filters['level'] = row.level
        backend_kwargs = {'filter_by_keys': filters, 'indexpath': ''}
        try:
            nwp_m = xr.open_mfdataset(nwp_files,
                                      concat_dim=['time', 'number', 'step'],
                                      engine='cfgrib',
                                      decode_timedelta=True, combine='nested',
                                      backend_kwargs=backend_kwargs)
        except:
            # see which dataset is missing the time coordinate
            for f_list in nwp_files:
                for f_list_inner in f_list:
                    for f in f_list_inner:
                        ds_f = xr.open_dataset(f, engine='cfgrib',
                                               decode_timedelta=True,
                                               backend_kwargs=backend_kwargs)
                        if 'time' not in ds_f.coords.keys():
                            print(ds_f)
                            message = f'{f}: {product}, {row.typeOfLevel}, {row.level}'
                            raise IncompleteDataException(message)
        if not len(nwp_m.coords.keys()):
            print(nwp_files)
            raise Exception(f'Missing data for {product}, {row.typeOfLevel}, {row.level}')
        nwp_ds.append(nwp_m)
    return nwp_ds

def open_herbie_dataset(var_df, dir, members, runs=None):
    '''Open a collection of Herbie files with xarray. All variables should have
    the same horizontal coordinates.
    '''
    #var_df = var_df.loc[~pd.isnull(var_df['product']), :]
    nwp_ds = []
    for product, product_vars in var_df.groupby('product'):
        product_txt = {'atmos.5': 'pgrb2a', 'atmos.25': 'pgrb2s', 'atmos.5b':
                       'pgrb2b'}[product]
        product_ds_list = get_nwp_product_vars(product_vars, dir, product_txt,
                                               members, runs=runs)
        nwp_ds.extend(product_ds_list)
    return merge_nwp_variables(nwp_ds)
