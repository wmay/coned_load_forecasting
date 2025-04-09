'''Extensions of the herbie package.

This module relies on the wgrib2 system package, which can be installed with
spack or conda. (I recommend spack.)
'''

import os, tempfile, hashlib, logging
from herbie import Herbie, FastHerbie, wgrib2
from herbie import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

log = logging.getLogger(__name__)

class RegionHerbie(Herbie):
    '''This is just like Herbie but with regional subsetting.
    '''
    def __init__(self, *args, extent=None, extent_name='region', **kwargs):
        self.extent = extent
        self.extent_name = extent_name
        self.download_full = False
        super().__init__(*args, **kwargs)

    def download(self, search=None, save_dir=None, overwrite=None, verbose=None,
                 *args, **kwargs):
        # check for existing file
        outFile = self.get_localFilePath(search)
        if save_dir is not None:
            outFile = (
                self.save_dir.expand()
                / self.model
                / f"{self.date:%Y%m%d}"
                / outFile.name
            )
        if overwrite is not None:
            self.overwrite = overwrite
        if outFile.exists() and not self.overwrite:
            if verbose:
                print(f"🌉 Already have local copy --> {outFile}")
            return outFile
        if not outFile.parent.is_dir():
            outFile.parent.mkdir(parents=True, exist_ok=True)
            print(f"👨🏻‍🏭 Created directory: [{outFile.parent}]")
        # download the full file into a temporary directory and subset with
        # wgrib2
        with tempfile.TemporaryDirectory() as tmp_dir:
            orig_save_dir = self.save_dir
            self.download_full = True # needed to get the correct path
            try:
                tmp_subdir = (
                    Path(tmp_dir)
                    / self.model
                    / f"{self.date:%Y%m%d}"
                )
                tmp_subdir.mkdir(parents=True, exist_ok=True)
                full_file = super().download(search=search, save_dir=tmp_dir,
                                             overwrite=overwrite,
                                             verbose=verbose, *args, **kwargs)
            except e:
                raise e
            finally:
                self.download_full = False
                self.save_dir = orig_save_dir # because download() changes it
            subset_file = wgrib2.region(full_file, self.extent,
                                        name=self.extent_name)
            #print(os.listdir(os.path.dirname(subset_file)))
            os.rename(subset_file, outFile)
            os.rename(str(subset_file) + '.idx', str(outFile) + '.idx')
        return outFile

    # this is the original function, but with the search hash changed
    def _get_localFilePath(self, search = None, *, searchString=None) -> Path:
        """Get full path to the local file."""
        # TODO: Remove this check for searString eventually
        if searchString is not None:
            warnings.warn(
                "The argument `searchString` was renamed `search`. Please update your scripts.",
                DeprecationWarning,
                stacklevel=2,
            )
            search = searchString

        # Predict the localFileName from the first model template SOURCE.
        localFilePath = (
            self.save_dir / self.model / f"{self.date:%Y%m%d}" / self.get_localFileName
        )

        # Check if any sources in a model template are "local"
        # (i.e., a custom template file)
        if any([i.startswith("local") for i in self.SOURCES.keys()]):
            localFilePath = next(
                (
                    Path(self.SOURCES[i])
                    for i in self.SOURCES
                    if i.startswith("local") and Path(self.SOURCES[i]).exists()
                ),
                localFilePath,
            )

        if search is not None:
            # Reassign the index DataFrame with the requested search
            # idx_df = self.inventory(search)

            # ======================================
            # Make a unique filename for the subset

            # Get a list of all GRIB message numbers. We will use this
            # in the output file name as a unique identifier.
            # all_grib_msg = "-".join([f"{i:g}" for i in idx_df.index])

            # To prevent "filename too long" error, create a hash to
            # that represents the file name and subseted variables to
            # shorten the name.

            # I want the files to still be sorted by date, fxx, and
            # subset fields, so include three separate hashes to similar
            # files will be sorted together.

            hash_date = hashlib.blake2b(
                f"{self.date:%Y%m%d%H%M}".encode(), digest_size=1
            ).hexdigest()

            hash_fxx = hashlib.blake2b(
                f"{self.fxx}".encode(), digest_size=1
            ).hexdigest()

            hash_label = hashlib.blake2b(
                search.encode(), digest_size=2
            ).hexdigest()

            # Prepend the filename with the hash label to distinguish it
            # from the full file. The hash label is a cryptic
            # representation of the GRIB messages in the subset.
            localFilePath = (
                localFilePath.parent
                / f"subset_{hash_date}{hash_fxx}{hash_label}__{localFilePath.name}"
            )

        return localFilePath

    def get_localFilePath(self, search=None):
        if self.download_full or self.extent is None:
            return super().get_localFilePath(search=search)
        path0 = self._get_localFilePath(search=search)
        new_name = self.extent_name + '_' + path0.name
        return path0.with_name(new_name)

class FastRegionHerbie(FastHerbie):
    '''FastHerbie, but using RegionHerbie instead of standard Herbie. This allows
    regional subsetting. Include the arguments `extent` and `extent_name` to
    specify regional subsetting.
    '''
    def __init__( self, DATES, fxx, *, max_threads=50, **kwargs):
        # self.DATES = _validate_DATES(DATES)
        # self.fxx = _validate_fxx(fxx)
        self.DATES = DATES
        self.fxx = fxx

        kwargs.setdefault("verbose", False)

        ################
        # Multithreading
        self.tasks = len(DATES) * len(fxx)
        threads = min(self.tasks, max_threads)
        log.info(f"🧵 Working on {self.tasks} tasks with {threads} threads.")

        self.objects = []
        with ThreadPoolExecutor(threads) as exe:
            futures = [
                exe.submit(RegionHerbie, date=DATE, fxx=f, **kwargs)
                for DATE in DATES
                for f in fxx
            ]

            # Return list of Herbie objects in order completed
            for future in as_completed(futures):
                if future.exception() is None:
                    self.objects.append(future.result())
                else:
                    log.error(f"Exception has occured : {future.exception()}")

        log.info(f"Number of Herbie objects: {len(self.objects)}")

        # Sort the list of Herbie objects by lead time then by date
        self.objects.sort(key=lambda H: H.fxx)
        self.objects.sort(key=lambda H: H.date)

        self.objects = self.objects

        # Which files exist?
        self.file_exists = [H for H in self.objects if H.grib is not None]
        self.file_not_exists = [H for H in self.objects if H.grib is None]

        if len(self.file_not_exists) > 0:
            log.warning(
                f"Could not find {len(self.file_not_exists)}/{len(self.file_exists)} GRIB files."
            )

class EnsHerbie:
    '''A tool for downloading ensemble forecasts, relying on `FastRegionHerbie`.
    '''
    
    def __init__(self, members, **kwargs):
        self.members = []
        for i, m in enumerate(members):
            print(f'Initializing ensemble member {str(i)} of {str(len(self.members))}')
            self.members.append(FastRegionHerbie(member=m, **kwargs))

    def download(self, *args, **kwargs):
        for i, m in enumerate(self.members):
            print(f'Downloading ensemble member {str(i)} of {str(len(self.members))}')
            m.download(*args, **kwargs)
