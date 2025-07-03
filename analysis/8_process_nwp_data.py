'''Prepare NWP data for model training.

The output is two files of daily predictors, one for both .5 and .25 degree
resolution.
'''

import os
os.chdir('..')
import numpy as np
import xarray as xr


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
    dbf = as_fahrenheit(db - 273.15)
    wbf = coned_wet_bulb(db, dp)
    return (dbf + wbf) / 2


# We start by creating 4 datasets-- forecast and analysis, .5 and .25 degrees.
# Then we will combine the forecast and analysis into a dataset ranging from -2
# to +7 days, to cover all the needs of the load forecasting models.

# Ensembles can be represented by the mean, except for TV/effective temp, which
# must be calculated individually.

def get_all_0p5(ds_string):
    '''Combine 0p5 and 0p5b datasets
    '''
    nc_0p5b = f'results/get_nwp_data/gefs_atmos0p5b_{ds_string}_members.nc'
    nc_0p5 = f'results/get_nwp_data/gefs_atmos0p5_{ds_string}.nc'
    ds_0p5b = xr.open_dataset(nc_0p5b).\
        mean('number')
    ds_0p5 = xr.open_dataset(nc_0p5).\
        squeeze().\
        drop_vars('number')
    return xr.merge([ds_0p5, ds_0p5b])

def get_pre_weather(anl, t):
    '''Given a forecast time, get the analysis dataset for weather occurring in
    the previous 2 days.
    '''
    # two days back, plus the offset from 8am EDT (12pm UTC) to the previous 9pm
    t_start = t - np.timedelta64(48 + 11, 'h')
    t_end = t
    out = anl.sel(time=slice(t_start, t_end))
    # fix the coordinates so the result can be easily merged
    # swap time with step (negative numbers)
    out['step'] = ((out['time'] - t) / np.timedelta64(1, 'h')).astype(int)
    out = out.swap_dims({'time': 'step'}).drop_vars('time').\
        expand_dims({'time': [t]}) # use t as time
    return out

def backfill_fct(fct, anl):
    '''Add analysis data to a forecast dataset, to add negative forecast steps
    '''
    pre_weathers = [ get_pre_weather(anl, t) for t in fct['time'].values ]
    pre_weather = xr.concat(pre_weathers, 'time')
    return xr.concat([pre_weather, fct], 'step')

def get_edt_day(step):
    '''Get day offset, with days going from 9pm-9pm EDT.
    '''
    # step 0 is 12UTC, -4 offset to EDT, +3 offset to shift to 9pm bounds
    edt9pm_day = (12 - 4 + 3 + step) // 24
    return edt9pm_day

def get_edt_hour(step):
    '''Get EDT hour from forecast step.
    '''
    # step 0 is 12UTC, -4 offset to EDT
    edt_hour = ((12 - 4) + step) % 24
    return edt_hour

def get_daily_avg(ds):
    '''Get daily averages, with days going from 9pm-9pm EDT.
    '''
    edt9pm_day = get_edt_day(ds['step'])
    ds = ds.assign_coords({'edt9pm_day': edt9pm_day})
    return ds.groupby('edt9pm_day').mean().sel(edt9pm_day=slice(None, 7))


# Let's start with the non-TV forecast data
gefs_0p5_fct = get_all_0p5('fct')

# organize the analysis
gefs_0p5_anl = get_all_0p5('anl').\
    squeeze().\
    drop_vars('step')
gefs_0p5_fill3 = get_all_0p5('fill3').\
    squeeze().\
    drop_vars('step')
# convert from forecast run time to valid time
gefs_0p5_fill3['time'] = gefs_0p5_fill3['time'] + np.timedelta64(3, 'h')
gefs_0p5_anl = xr.concat([gefs_0p5_anl, gefs_0p5_fill3], dim='time').\
    sortby('time')

# now connect the fct and anl datasets into a continuous -2 to +7 day dataset
gefs_0p5 = backfill_fct(gefs_0p5_fct, gefs_0p5_anl)

# now that we have all the 0p5 data, get daily averages

# ---------------
# TO DO: We are currently missing 12z from the analysis files. This needs to be
# fixed in script 7
# ---------------

def get_3day_wmean(ds):
    '''Get a 3-day weighted mean, using ConEd's TV weights.
    '''
    # get everything to line up
    ds_m1 = ds.copy()
    ds_m2 = ds.copy()
    ds_m1['edt9pm_day'] = ds['edt9pm_day'] + 1
    ds_m2['edt9pm_day'] = ds['edt9pm_day'] + 2
    return ds * .7 + ds_m1 * .2 + ds_m2 * .1

gefs_0p5_daily = get_daily_avg(gefs_0p5)
gefs_0p5_daily.to_netcdf('results/process_nwp_data/gefs_0p5_daily.nc')
gefs_0p5_3day_wmean = get_3day_wmean(gefs_0p5_daily)
gefs_0p5_3day_wmean.to_netcdf('results/process_nwp_data/gefs_0p5_3day_wmean.nc')


# now the same for 0p25
gefs_0p25_fct = xr.open_dataset('results/get_nwp_data/gefs_atmos0p25_fct.nc').\
    squeeze().\
    drop_vars('number')

# organize the analysis
gefs_0p25_anl = xr.open_dataset('results/get_nwp_data/gefs_atmos0p25_anl.nc').\
    squeeze().\
    drop_vars('step')
# add the variables not in the analysis files
gefs_0p25_fill6 = xr.open_dataset('results/get_nwp_data/gefs_atmos0p25_fill6.nc').\
    squeeze().\
    drop_vars('step')
gefs_0p25_fill6['time'] = gefs_0p25_fill6['time'] + np.timedelta64(6, 'h')
gefs_0p25_anl = xr.merge([gefs_0p25_anl, gefs_0p25_fill6])
# add the 3-hour values
gefs_0p25_fill3 = xr.open_dataset('results/get_nwp_data/gefs_atmos0p25_fill3.nc').\
    squeeze().\
    drop_vars('step')
gefs_0p25_fill3['time'] = gefs_0p25_fill3['time'] + np.timedelta64(3, 'h')
gefs_0p25_anl = xr.concat([gefs_0p25_anl, gefs_0p25_fill3], dim='time').\
    sortby('time').\
    drop_vars('number')

# now connect the fct and anl datasets into a continuous -2 to +7 day dataset
gefs_0p25 = backfill_fct(gefs_0p25_fct, gefs_0p25_anl)
# now that we have all the 0p25 data, get daily averages
gefs_0p25_daily = get_daily_avg(gefs_0p25)
# gefs_0p25_daily.to_netcdf('results/process_nwp_data/gefs_0p25_daily.nc')
# still need to add effective temperature



# No TV needed from analysis, because it's filled in with observation data. Only
# need forecast TV values
gefs_atmos0p25_fct_members = xr.open_dataset('results/get_nwp_data/gefs_atmos0p25_fct_members.nc')

# Somehow we got temperature values above 29000, which I hope is a misreading of
# missing data.
# (ds2['t2m'] > 29000).sum()
# ds2['t2m'].values = np.where(ds2['t2m'] > 500, np.nan, ds2['t2m'].values)
# ds2['d2m'].values = np.where(ds2['d2m'] > 500, np.nan, ds2['d2m'].values)
# fixed!

# # how many missing values did we end up with?
# np.isnan(gefs_atmos0p25_fct_members['t2m']).any(['latitude', 'longitude']).sum()
# # 23 files-- almost none! Not even one set of ensemble members?

edt_coords = {'edt9pm_day': get_edt_day(gefs_atmos0p25_fct_members['step']),
              'edt_hour': get_edt_hour(gefs_atmos0p25_fct_members['step'])}
gefs_atmos0p25_fct_members = gefs_atmos0p25_fct_members.assign_coords(edt_coords)

# # what hours did missing data occur?
# missing_files = np.isnan(gefs_atmos0p25_fct_members['t2m']).any(['latitude', 'longitude'])
# missing_files.sum(['time', 'number']).\
#     set_index(step=['edt_hour', 'edt9pm_day']).unstack().\
#     sel(edt_hour=slice(9, 21)).sum('edt9pm_day').\
#     astype(int).to_dataframe()
# # we have 6 missing files between 9am and 9pm

# use only 9am to 9pm
gefs_atmos0p25_fct_members = gefs_atmos0p25_fct_members.\
    where((gefs_atmos0p25_fct_members['edt_hour'] >= 9) &
          (gefs_atmos0p25_fct_members['edt_hour'] <= 21), drop=True)

gefs_fct_eff_temp = coned_eff_temp(gefs_atmos0p25_fct_members['t2m'],
                                   gefs_atmos0p25_fct_members['d2m']).\
    groupby('edt9pm_day').max(skipna=False)

# np.isnan(gefs_fct_eff_temp).sum() / gefs_fct_eff_temp.size

gefs_0p25_daily['eff_temp'] = gefs_fct_eff_temp.mean('number').transpose('edt9pm_day', 'time', 'latitude', 'longitude')
gefs_0p25_daily.to_netcdf('results/process_nwp_data/gefs_0p25_daily.nc')
gefs_0p25_3day_wmean = get_3day_wmean(gefs_0p25_daily)
gefs_0p25_3day_wmean.to_netcdf('results/process_nwp_data/gefs_0p25_3day_wmean.nc')


# Also need to do 3-day summaries here, to get the correct standard deviations

# pad with zeros to simplify the code
padding = gefs_fct_eff_temp.isel(edt9pm_day=slice(None, 2)).copy()
padding['edt9pm_day'] = padding['edt9pm_day'] - 2
padding[:,:,:,:,:] = 0

gefs_eff_temp = xr.concat([padding, gefs_fct_eff_temp], 'edt9pm_day')

gefs_tv_members = gefs_eff_temp.sel(edt9pm_day=slice(0, None)) * .7 +\
    gefs_eff_temp.sel(edt9pm_day=slice(-1, 6)).values * .2 +\
    gefs_eff_temp.sel(edt9pm_day=slice(None, 5)).values * .1
gefs_tv_members = gefs_tv_members.rename('TV')

gefs_tv = xr.merge([gefs_tv_members.mean('number'),
                    gefs_tv_members.std('number').rename('TV_sd')])
gefs_tv.to_netcdf('results/process_nwp_data/gefs_tv.nc')

# # let's see how much forecast variation there is, split by forecast day
# gefs_tv_members.std('number').mean(['time', 'latitude', 'longitude']).to_dataframe()
#          eff_temp
# est_day          
# 2        1.347529
# 3        1.660509
# 4        1.954607
# 5        2.288087
# 6        2.632163
# 7        2.987738

# import matplotlib.pyplot as plt
# fig, axs = plt.subplots(ncols=6)
# for i in range(2, 8):
#     tv.sel(est_day=i).std('number').plot(ax=axs[i - 2])
