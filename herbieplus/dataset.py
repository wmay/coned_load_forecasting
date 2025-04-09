'''Tools for organizing an entire dataset with Herbie.
'''

import os, hashlib, functools
from concurrent.futures import ThreadPoolExecutor, as_completed
from herbie import Herbie
from .herbie import RegionHerbie

class HerbieCollection:
    '''Download a (potentially large) collection of Herbie files. This is
    similar to `FastHerbie`, but instead of downloading all index files, it
    downloads only the files that are currently missing. Unlike `Herbie`, it
    checks for files before downloading anything, which can save large amounts
    of time when a large collection is partially downloaded.
    '''
    def __init__(self, DATES, model, product, search_string, fxx, members=None,
                 save_dir=None, extent=None, extent_name=None):
        self.DATES = DATES
        self.model = model
        self.product = product
        self.search_string = search_string
        self.fxx = fxx
        self.members = members
        self.save_dir = save_dir
        self.extent = extent
        self.extent_name = extent_name

    def download(self, threads=2):
        '''Download the missing files using Herbie.'''
        for date in self.DATES:
            n_date_missing = self.count_run_missing(date)
            if not n_date_missing:
                continue
            print(f'- Downloading {n_date_missing} files for {date}')
            if len(self.members) == 1:
                self.download_fxx(date, threads=threads)
            else:
                for fxx in self.fxx:
                    n_fxx_missing = self.count_fxx_missing(date, fxx)
                    if not n_fxx_missing:
                        continue
                    print(f'Downloading {n_fxx_missing} files for step {fxx}')
                    self.download_members(date, fxx, threads=threads)

    def download_fxx(self, date, threads=2):
        '''Download the missing files for a given date.'''
        outFiles = []
        m = self.members[0]
        with ThreadPoolExecutor(threads) as exe:
            futures = [
                exe.submit(self.download_single, date, fh, m)
                for fh in self.fxx
                if not self.files_exist(date, fh, m)
            ]

    def download_members(self, date, fxx, threads=2):
        '''Download the missing members for a given date and forecast step.'''
        outFiles = []
        with ThreadPoolExecutor(threads) as exe:
            futures = [
                exe.submit(self.download_single, date, fxx, m)
                for m in self.members
                if not self.files_exist(date, fxx, m)
            ]
            # # Return list of Herbie objects in order completed
            # for future in as_completed(futures):
            #     if future.exception() is None:
            #         pass
            #     else:
            #         log.error(f"Exception has occured : {future.exception()}")

    def count_run_missing(self, date):
        '''Count the number of missing files for a given date.'''
        n_missing = 0
        for fxx in self.fxx:
            n_missing += self.count_fxx_missing(date, fxx)
        return n_missing

    def count_fxx_missing(self, date, fxx):
        '''Count the number of missing files for a given date and forecast
        step.'''
        n_missing = 0
        for member in self.members:
            if not self.files_exist(date, fxx, member):
                n_missing += 1
        return n_missing

    def download_single(self, date, fxx, member):
        '''Download a single file using Herbie.'''
        h = RegionHerbie(date, model=self.model, product=self.product, fxx=fxx,
                         save_dir=self.save_dir, member=member,
                         extent=self.extent, extent_name=self.extent_name,
                         verbose=False)
        h.download(self.search_string)

    def file_path(self, date, fxx, member):
        # This should ideally use the model template to create paths, but for
        # now just use string interpolation. Example:
        # 'nwp_data/gefs/20230101/nyc_subset_18efda92__gec00.t00z.pgrb2a.0p50.f000'
        # gefs_member_txt = 'gec00' if member == 0 else f'gep{member:02d}'
        gefs_member_txt = 'geavg'
        product_txt = {'atmos.5': 'pgrb2a.0p50',
                       'atmos.5b': 'pgrb2b.0p50',
                       'atmos.25': 'pgrb2s.0p25'}[self.product]
        fdir = f'{self.save_dir}/{self.model}/{date.strftime("%Y%m%d")}/'
        fhash = self.file_hash(date, fxx)
        fname = f'{self.extent_name}_subset_{fhash}__{gefs_member_txt}.t{date.strftime("%H")}z.{product_txt}.f{fxx:03d}'
        return fdir + fname

    def files_exist(self, date, fxx, member):
        '''Check if the file exists to see if it should be downloaded.'''
        file_path = self.file_path(date, fxx, member)
        idx_path = file_path + '.idx'
        files_exist = os.path.exists(file_path) and os.path.exists(idx_path)
        return files_exist

    @functools.cached_property
    def search_hash(self):
        # get the inventory from the Herbie object for the first file
        h1 = Herbie(self.DATES[0], model=self.model, product=self.product,
                    save_dir=self.save_dir, member=self.members[0])
        # Not using the inventory for the hash because it doesn't seem to be
        # consistent between files:
        # idx_df = h1.inventory(self.search_string)
        # all_grib_msg = "-".join([f"{i:g}" for i in idx_df.index])
        hash_label = hashlib.blake2b(
            #all_grib_msg.encode(), digest_size=2
            # match RegionHerbie
            self.search_string.encode(), digest_size=2
        ).hexdigest()
        return hash_label

    def file_hash(self, date, fxx):
        '''Generate a hash for the file name.'''
        hash_date = hashlib.blake2b(
            f"{date:%Y%m%d%H%M}".encode(), digest_size=1
        ).hexdigest()
        hash_fxx = hashlib.blake2b(
            f"{fxx}".encode(), digest_size=1
        ).hexdigest()
        return f"{hash_date}{hash_fxx}{self.search_hash}"
