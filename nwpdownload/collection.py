'''Organize collections of NWP files.
'''

# import os, tempfile, shutil
import tempfile, shutil
from datetime import datetime
from humanize import naturalsize
import numpy as np
import dask
import dask.bag as db
from herbie import Herbie, wgrib2
from .nwppath import NwpPath
from .nwpdownloader import NwpDownloader

class NwpCollection:
    '''Download a (potentially large) collection of NWP files. This is
    similar to `FastHerbie`, but instead of downloading all index files, it
    downloads only the files that are currently missing. Unlike `Herbie`, it
    checks for files before downloading anything, which can save large amounts
    of time when a large collection is partially downloaded.
    '''
    def __init__(self, DATES, model, product, search_string, fxx, members=None,
                 save_dir=None, extent=None):
        self.DATES = DATES
        self.model = model
        self.product = product
        self.search_string = search_string
        self.fxx = fxx
        self.members = members
        self.save_dir = save_dir
        self.extent = extent

    def get_status(self):
        '''Print the status of the collection.
        '''
        n_files = len(self.DATES) * len(self.fxx) * len(self.members)
        full_download_size = self.file_size * n_files
        print(f'Complete download size (approximate): {naturalsize(full_download_size)}')
        if self.extent is not None:
            print('(Size on disk will be smaller due to regional subsetting.)')
        # make a matrix representing the download status of each file. 0 is
        # missing and 1 is downloaded
        download_status = np.zeros((len(self.DATES), len(self.fxx),
                                    len(self.members)), dtype=int)
        for i, date in enumerate(self.DATES):
            for j, fxx in enumerate(self.fxx):
                for k, member in enumerate(self.members):
                    if self._file_exists(date, fxx=fxx, member=member):
                        download_status[i, j, k] = 1
        n_remaining = n_files - download_status.sum()
        remaining_download_size = self.file_size * n_remaining
        print(f'{n_files - n_remaining} of {n_files} files downloaded')
        print(f'Remaining download size (approximate): {naturalsize(remaining_download_size)}')
        out = {'n_files': n_files,
               'n_remaining': n_remaining,
               'remaining_download_size': remaining_download_size,
               'download_array': download_status}
        return out

    def download(self, overwrite=False):
        '''Download all remaining files in the collection.
        '''
        status = self.get_status()
        if not status['n_remaining']:
            print('Nothing to download.')
            return dask.compute()
        use_bag = status['n_remaining'] > 50000
        remaining_coords = np.stack(np.where(~status['download_array'])).T
        # with tempfile.TemporaryDirectory() as tmp_dir:
        #     print(tmp_dir)
        #     # tasks = [ dask.delayed(self._download_from_coords)(coords, tmp_dir)
        #     #           for coords in list(remaining_coords) ]
        #     tasks = [ self._download_from_coords(coords, tmp_dir)
        #               for coords in list(remaining_coords) ]
        #     return tasks
        #     #return dask.compute(*tasks)
        start_time = datetime.now()
        # with tempfile.TemporaryDirectory() as tmp_dir:
        # print(tmp_dir)
        if use_bag:
            print('organizing dask bag')
            # bag = db.from_sequence(tasks, npartitions=1000)
            # bag.compute()
            bag = db.from_sequence(remaining_coords, npartitions=1000)
            # bag.map(lambda x: self._extract_region(self._download_from_coords(x, tmp_dir), x))
            bag.map(self._download_and_extract)
            bag.compute()
        else:
            tasks = []
            for coords in remaining_coords:
                # full_file = dask.delayed(self._download_from_coords)(coords, tmp_dir)
                # region_file = dask.delayed(self._extract_region)(full_file, coords)
                # tasks.append(region_file)
                out_file = dask.delayed(self._download_and_extract)(coords)
                tasks.append(out_file)
            dask.compute(*tasks)
        time_elapsed = datetime.now() - start_time
        print(f'Total run time: {time_elapsed}')

    # def _download_from_coords(self, coords, tmp_dir):
    #     '''Download a single file using the coordinates from the download
    #     status matrix.
    #     '''
    #     date = self.DATES[coords[0]]
    #     fxx = self.fxx[coords[1]]
    #     member = self.members[coords[2]]
    #     self._download(tmp_dir, date, fxx=fxx, member=member)

    def _download_and_extract(self, coords):
        '''Download a single file using the coordinates from the download
        status matrix.
        '''
        date = self.DATES[coords[0]]
        fxx = self.fxx[coords[1]]
        member = self.members[coords[2]]
        nwp_file = NwpPath(date, model=self.model, product=self.product,
                           save_dir=self.save_dir, fxx=fxx, member=member)
        out_path = nwp_file.get_localFilePath()
        # create an individual directory for each download
        with tempfile.TemporaryDirectory() as tmp_dir:
            full_file = self._download(tmp_dir, date, fxx=fxx, member=member)
            region_file = wgrib2.region(full_file, self.extent)
            if not out_path.parent.is_dir():
                out_path.parent.mkdir(parents=True, exist_ok=True)
            # os.rename(region_file, out_path)
            # os.rename(str(region_file) + '.idx', str(out_path) + '.idx')
            shutil.move(region_file, out_path)
            shutil.move(str(region_file) + '.idx', str(out_path) + '.idx')
        return out_path

    # def _download_from_coords(self, coords, tmp_dir):
    #     '''Download a single file using the coordinates from the download
    #     status matrix.
    #     '''
    #     date = self.DATES[coords[0]]
    #     fxx = self.fxx[coords[1]]
    #     member = self.members[coords[2]]
    #     # create an individual directory for each download
    #     tmp_subdir = tempfile.mkdtemp(dir=tmp_dir)
    #     return self._download(tmp_subdir, date, fxx=fxx, member=member)

    # def _extract_region(self, full_file, coords):
    #     '''Call `wgrib2.region` in a manner that works with `dask.delayed`.
    #     '''
    #     date = self.DATES[coords[0]]
    #     fxx = self.fxx[coords[1]]
    #     member = self.members[coords[2]]
    #     nwp_file = NwpPath(date, model=self.model, product=self.product,
    #                        save_dir=self.save_dir, fxx=fxx, member=member)
    #     out_path = nwp_file.get_localFilePath()
    #     if not out_path.parent.is_dir():
    #         out_path.parent.mkdir(parents=True, exist_ok=True)
    #     region_file = wgrib2.region(full_file, self.extent)
    #     os.rename(region_file, out_path)
    #     os.rename(str(region_file) + '.idx', str(out_path) + '.idx')
    #     shutil.rmtree(full_file.parent)

    def _download(self, tmp_dir, *args, **kwargs):
        '''Download an NWP grib2 file using Herbie, and subset it with wgrib2.
        Discards the full file.
        '''
        # nwp_file = NwpPath(*args, **kwargs, model=self.model,
        #                    product=self.product, save_dir=self.save_dir)
        # out_path = nwp_file.get_localFilePath()
        # if not out_path.parent.is_dir():
        #     out_path.parent.mkdir(parents=True, exist_ok=True)
        # with tempfile.TemporaryDirectory() as tmp_dir:
        # tmp_archive = Herbie(*args, **kwargs, save_dir=tmp_dir,
        #                      model=self.model, product=self.product,
        #                      verbose=False)
        tmp_archive = NwpDownloader(*args, **kwargs, save_dir=tmp_dir,
                                    model=self.model, product=self.product,
                                    verbose=False)
        # herbie will drive me insane with messages if I don't create this
        # directory
        tmp_path = tmp_archive.get_localFilePath()
        tmp_path.parent.mkdir(parents=True, exist_ok=True)
        full_file = tmp_archive.download(self.search_string)
        return full_file
        # region_file = wgrib2.region(full_file, self.extent)
        # full_file = dask.delayed(tmp_archive.download)(self.search_string)
        # region_file = self._wgrib2_region_delayed(full_file, out_path)
        # os.rename(region_file, out_path)
        # os.rename(str(region_file) + '.idx', str(out_path) + '.idx')

    def _file_exists(self, *args, **kwargs):
        '''Check if the file exists for a given date, fxx, and member.
        '''
        nwp_file = NwpPath(*args, **kwargs, model=self.model,
                           product=self.product, save_dir=self.save_dir)
        return nwp_file.get_localFilePath().exists()

    @property
    def file_size(self):
        '''Use the inventory to calculate the file size.
        '''
        # this temp directory prevents Herbie from reading from local index
        # files, which may not match the source file
        with tempfile.TemporaryDirectory() as tmp_dir:
            h = Herbie(self.DATES[0], model=self.model, product=self.product,
                       fxx=self.fxx[0], member=self.members[0],
                       save_dir=tmp_dir, verbose=False)
        inv = h.inventory(search=self.search_string)
        return (inv['end_byte'] - inv['start_byte'] + 1).sum()

    def collection_size(self, humanize=True):
        '''Calculate the size of the collection.
        '''
        n_files = len(self.DATES) * len(self.fxx) * len(self.members)
        full_download_size = self.file_size * n_files
        if humanize:
            return naturalsize(full_download_size)
        else:
            return full_download_size
