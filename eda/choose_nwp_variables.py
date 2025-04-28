'''See what NWP variables are available and if they can be downloaded in a
practical way.
'''

import os
os.chdir('..')
import numpy as np
import pandas as pd
from herbie import Herbie
from nwpdownload.collection import NwpCollection
from dask.distributed import Client

# GEFS

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
# inv = archive.inventory()
# with pd.option_context('display.max_rows', None):
#         print(inv.sort_values(by=['variable']).iloc[:, 6:])
# inv.loc[inv4['variable'] == 'SHTFL', :].iloc[:, 6:]

# check to see the download requirements for various data collections

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


# what if I want all members for TV calculations?
tv_search = '|'.join([':TMP:2 m above ground:', ':DPT:2 m above ground:'])
tv_collection = NwpCollection(runs, 'gefs', 'atmos.25', tv_search, fxx,
                              members=range(0, 31), save_dir='nwp_test',
                              extent=nyc_extent)
tv_collection.collection_size()
tv_collection.get_status()
# 1.7 TB, 1429596 (1.4 million) files, holy jesus

# what if I want :SPFH:850 mb:, only available as individual members?
spfh850_collection = NwpCollection(runs, 'gefs', 'atmos.5b', ':SPFH:850 mb:', fxx,
                                   members=range(0, 31), save_dir='nwp_test',
                                   extent=nyc_extent)
spfh850_collection.collection_size()
spfh850_collection.get_status()
# 306.2 GB, 1429596 (1.4 million) files

# all 0.25 resolution variables from Rasp and Lerch 2018
params_0p25 = rasp_params.loc[rasp_params['product'] == 'atmos.25', :]
search_0p25 = '|'.join(params_0p25['search'])
collection_0p25 = NwpCollection(runs, 'gefs', 'atmos.25', search_0p25, fxx,
                                members=['avg'], save_dir='nwp_test',
                                extent=nyc_extent)
collection_0p25.collection_size()
collection_0p25.get_status()
# 321.4 GB, 46116 files

# all 0.5 resolution variables from Rasp and Lerch 2018, excluding atmos.5b
# variables since there's no average available
params_0p5 = rasp_params.loc[rasp_params['product'] == 'atmos.5', :]
search_0p5 = '|'.join(params_0p5['search'])
collection_0p5 = NwpCollection(runs, 'gefs', 'atmos.5', search_0p5, fxx,
                                members=['avg'], save_dir='nwp_test',
                                extent=nyc_extent)
collection_0p5.collection_size()
collection_0p5.get_status()
# 59.4 GB, 46116 files


# let's do a test to see what kind of download speed I can achieve

# following https://examples.dask.org/applications/embarrassingly-parallel.html#Start-Dask-Client-for-Dashboard
# client = Client()
client = Client(threads_per_worker=4, n_workers=1)
client
client.dashboard_link
collection_0p5.download() # accidentally used the wrong collection here
client.shutdown()
# 3 to 4.5MB/s

client = Client(threads_per_worker=8, n_workers=1)
client
client.dashboard_link
collection_0p25.download()
client.shutdown()
# 10 to 15MB/s

client = Client(threads_per_worker=16, n_workers=1)
client
client.dashboard_link
collection_0p25.download()
client.shutdown()
# 15 to 17MB/s
# this seems like I'm consistently hitting a bandwidth limit

# UserWarning: Sending large graph of size 14.79 MiB.
# This may cause some slowdown.
# Consider loading the data with Dask directly
#  or using futures or delayed objects to embed the data into the graph without repetition.
# See also https://docs.dask.org/en/stable/best-practices.html#load-data-with-dask for more information.
