'''Optimized Herbie-based downloader.
'''

import os, warnings, requests, json, functools
import urllib.request
import pandas as pd
from io import StringIO
from typing import Literal, Optional, Union
from herbie.core import log, wgrib2_idx
from herbie.misc import ANSI
from herbie import Path, Herbie
from .nwppath import NwpPath

class NwpDownloader(NwpPath):
    '''Herbie with download optimizations. Skips some unnecessary checks (same
    as NwpPath), and adds timeout limits to requests.
    '''

    # Use timeout options for curl to avoid hanging.
    def download(
        self,
        search: Optional[str] = None,
        *,
        searchString=None,
        source: Optional[str] = None,
        save_dir: Optional[Union[str, Path]] = None,
        overwrite: Optional[bool] = None,
        verbose: Optional[bool] = None,
        errors: Literal["warn", "raise"] = "warn",
    ) -> Path:
        """
        Download file from source.

        TODO: When we download a full file, the value of self.grib and
        TODO:   self.grib_source should change to represent the local file.

        Subsetting by variable follows the same principles described here:
        https://www.cpc.ncep.noaa.gov/products/wesley/fast_downloading_grib.html

        Parameters
        ----------
        search : str
            If None, download the full file. Else, use regex to subset
            the file by specific variables and levels.
            Read more in the user guide:
            https://herbie.readthedocs.io/en/latest/user_guide/tutorial/search.html
        source : {'nomads', 'aws', 'google', 'azure', 'pando', 'pando2'}
            If None, download GRIB2 file from self.grib2 which is
            the first location the GRIB2 file was found from the
            priority lists when this class was initialized. Else, you
            may specify the source to force downloading it from a
            different location.
        save_dir : str or pathlib.Path
            Location to save the model output files.
            If None, uses the default or path specified in __init__.
            Else, changes the path files are saved.
        overwrite : bool
            If True, overwrite existing files. Default will skip
            downloading if the full file exists. Not applicable when
            when search is not None because file subsets might
            be unique.
        errors : {'warn', 'raise'}
            When an error occurs, send a warning or raise a value error.

        """

        def _reporthook(a, b, c):
            """
            Print download progress in megabytes.

            Parameters
            ----------
            a : Chunk number
            b : Maximum chunk size
            c : Total size of the download
            """
            chunk_progress = a * b / c * 100
            total_size_MB = c / 1000000.0
            if verbose:
                print(
                    f"\r🚛💨  Download Progress: {chunk_progress:.2f}% of {total_size_MB:.1f} MB\r",
                    end="",
                )

        def subset(search, outFile):
            """Download a subset specified by the regex search."""
            # TODO: Maybe optimize downloading multiple subsets with MultiThreading

            # TODO An alternative to downloadling subset with curl is
            # TODO  to use the request module directly.
            # TODO  >> headers = dict(Range=f"bytes={start_bytes}-{end_bytes}")
            # TODO  >> r = requests.get(grib_url, headers=headers)

            grib_source = self.grib
            if hasattr(grib_source, "as_posix") and grib_source.exists():
                # The GRIB source is local. Curl the local file
                # See https://stackoverflow.com/a/21023161/2383070
                grib_source = f"file://{str(self.grib)}"
            if verbose:
                print(
                    f"📇 Download subset: {self.__repr__()}{' ':60s}\n cURL from {grib_source}"
                )

            # -----------------------------------------------------
            # Download subsets of the file by byte range with cURL.
            #  Instead of using a single curl command for each row,
            #  group adjacent messages in the same curl command.

            # Find index groupings
            idx_df = self.inventory(search).copy()
            if verbose:
                print(
                    f"Found {ANSI.bold}{ANSI.green}{len(idx_df)}{ANSI.reset} grib messages."
                )
            idx_df["download_groups"] = idx_df.grib_message.diff().ne(1).cumsum()

            # cURL subsets of each group
            for i, curl_group in idx_df.groupby("download_groups"):
                if verbose:
                    print(f"Download subset group {i}")

                if verbose:
                    for _, row in curl_group.iterrows():
                        print(
                            f"  {row.grib_message:<3g} {ANSI.orange}{row.search_this}{ANSI.reset}"
                        )

                range = f"{curl_group.start_byte.min():.0f}-{curl_group.end_byte.max():.0f}".replace(
                    "nan", ""
                )

                if curl_group.end_byte.max() - curl_group.start_byte.min() < 0:
                    # The byte range for GRIB submessages (like in the
                    # RAP model's UGRD/VGRD) need to be handled differently.
                    # See https://github.com/blaylockbk/Herbie/issues/259
                    if verbose:
                        print(f"  ERROR: Invalid cURL range {range}; Skip message.")
                    continue

                if i == 1:
                    # If we are working on the first item, overwrite the existing file...
                    curl = f'''curl -s --connect-timeout 3 --max-time 6 --range {range} "{grib_source}" > "{outFile}"'''
                else:
                    # ...all other messages are appended to the subset file.
                    curl = f'''curl -s --connect-timeout 3 --max-time 6 --range {range} "{grib_source}" >> "{outFile}"'''

                if verbose:
                    print(curl)
                os.system(curl)

            if verbose:
                print(f"💾 Saved the subset to {outFile}")

        # TODO: Remove this eventually
        if searchString is not None:
            warnings.warn(
                "The argument `searchString` was renamed `search`. Please update your scripts.",
                DeprecationWarning,
                stacklevel=2,
            )
            search = searchString

        # This overrides the save_dir specified in __init__
        if save_dir is not None:
            self.save_dir = Path(save_dir).expand()

        if not hasattr(Path(self.save_dir).expand(), "exists"):
            self.save_dir = Path(self.save_dir).expand()

        # If the file exists in the localPath and we don't want to
        # overwrite, then we don't need to download it.
        outFile = self.get_localFilePath(search)

        if save_dir is not None:
            # Looks like the save_dir was changed.
            outFile = (
                self.save_dir.expand()
                / self.model
                / f"{self.date:%Y%m%d}"
                / outFile.name
            )

        # This overrides the overwrite specified in __init__
        if overwrite is not None:
            self.overwrite = overwrite

        # This overrides the verbose specified in __init__
        if verbose is not None:
            self.verbose = verbose

        if outFile.exists() and not self.overwrite:
            if verbose:
                print(f"🌉 Already have local copy --> {outFile}")
            return outFile

        if self.overwrite and self.grib_source.startswith("local"):
            # Search for the grib files on the remote archives again
            self.grib, self.grib_source = self.find_grib(overwrite=True)
            self.idx, self.idx_source = self.find_idx()
            print(f"Overwrite local file with file from [{self.grib_source}]")

        # Check that data exists
        if self.grib is None:
            msg = f"🦨 GRIB2 file not found: {self.model=} {self.date=} {self.fxx=}"
            if errors == "warn":
                log.warning(msg)
                return  # Can't download anything without a GRIB file URL.
            elif errors == "raise":
                raise ValueError(msg)
        if self.idx is None and search is not None:
            msg = f"🦨 Index file not found; cannot download subset: {self.model=} {self.date=} {self.fxx=}"
            if errors == "warn":
                log.warning(
                    msg + " I will download the full file because I cannot subset."
                )
            elif errors == "raise":
                raise ValueError(msg)

        if source is not None:
            # Force download from a specified source and not from first in priority
            self.grib = self.SOURCES[source]

        # Create directory if it doesn't exist
        if not outFile.parent.is_dir():
            outFile.parent.mkdir(parents=True, exist_ok=True)
            print(f"👨🏻‍🏭 Created directory: [{outFile.parent}]")

        # ===============
        # Do the Download
        # ===============
        if search in [None, ":"] or self.idx is None:
            # Download the full file from remote source
            urllib.request.urlretrieve(self.grib, outFile, _reporthook)

            original_source = self.grib

            self.grib = outFile
            # self.grib_source = "local"  # ?? Why did I turn this off?

            if verbose:
                print(
                    f"✅ Success! Downloaded {self.model.upper()} from \033[38;5;202m{self.grib_source:20s}\033[0m\n\tsrc: {original_source}\n\tdst: {outFile}"
                )

        else:
            # Download a subset of the file
            subset(search, outFile)

        return outFile

    # Use timeout options for `requests.get` to avoid hanging.
    @functools.cached_property
    def index_as_dataframe(self) -> pd.DataFrame:
        """Read and cache the full index file."""
        if self.grib_source == "local" and wgrib2:
            # Generate IDX inventory with wgrib2
            self.idx = StringIO(wgrib2_idx(self.get_localFilePath()))
            self.idx_source = "generated"
            self.IDX_STYLE = "wgrib2"
        elif self.idx is None:
            if self.grib_source == "local":
                # Use wgrib2 to get the index file if the file is local
                log.info("🧙🏻‍♂️ I'll use wgrib2 to create the missing index file.")
                self.idx = StringIO(wgrib2_idx(self.get_localFilePath()))
                self.IDX_STYLE = "wgrib2"
            else:
                raise ValueError(
                    f"\nNo index file was found for {self.grib}\n"
                    f"Download the full file first (with `H.download()`).\n"
                    f"You will need to remake the Herbie object (H = `Herbie()`)\n"
                    f"or delete this cached property: `del H.index_as_dataframe()`"
                )
        if self.idx is None:
            raise ValueError(f"No index file found for {self.grib}.")

        if self.IDX_STYLE == "wgrib2":
            # Sometimes idx lines end in ':', other times it doesn't (in some Pando files).
            # https://pando-rgw01.chpc.utah.edu/hrrr/sfc/20180101/hrrr.t00z.wrfsfcf00.grib2.idx
            # https://noaa-hrrr-bdp-pds.s3.amazonaws.com/hrrr.20210101/conus/hrrr.t00z.wrfsfcf00.grib2.idx
            # Sometimes idx has more than the standard messages
            # https://noaa-nbm-grib2-pds.s3.amazonaws.com/blend.20210711/13/core/blend.t13z.core.f001.co.grib2.idx
            if self.idx_source in ["local", "generated"]:
                read_this_idx = self.idx
            else:
                read_this_idx = None
                response = requests.get(self.idx, timeout=3)
                if response.status_code != 200:
                    response.raise_for_status()
                    response.close()
                    raise ValueError(
                        f"\nCant open index file {self.idx}\n"
                        f"Download the full file first (with `H.download()`).\n"
                        f"You will need to remake the Herbie object (H = `Herbie()`)\n"
                        f"or delete this cached property: `del H.index_as_dataframe()`"
                    )
                read_this_idx = StringIO(response.text)
                response.close()

            df = pd.read_csv(
                read_this_idx,
                sep=":",
                names=[
                    "grib_message",
                    "start_byte",
                    "reference_time",
                    "variable",
                    "level",
                    "forecast_time",
                    "?",
                    "??",
                    "???",
                ],
            )

            # Format the DataFrame
            df["reference_time"] = pd.to_datetime(
                df.reference_time, format="d=%Y%m%d%H"
            )
            df["valid_time"] = df["reference_time"] + pd.to_timedelta(f"{self.fxx}h")
            df["start_byte"] = df["start_byte"].astype(int)
            df["end_byte"] = df["start_byte"].shift(-1) - 1
            df["range"] = df.apply(
                lambda x: f"{x.start_byte:.0f}-{x.end_byte:.0f}".replace("nan", ""),
                axis=1,
            )
            df = df.reindex(
                columns=[
                    "grib_message",
                    "start_byte",
                    "end_byte",
                    "range",
                    "reference_time",
                    "valid_time",
                    "variable",
                    "level",
                    "forecast_time",
                    "?",
                    "??",
                    "???",
                ]
            )

            df = df.dropna(how="all", axis=1)

            df["search_this"] = (
                df.loc[:, "variable":]
                .astype(str)
                .apply(
                    lambda x: ":" + ":".join(x).rstrip(":").replace(":nan:", ":"),
                    axis=1,
                )
            )

        if self.IDX_STYLE == "eccodes":
            # eccodes keywords explained here:
            # https://confluence.ecmwf.int/display/UDOC/Identification+keywords

            r = requests.get(self.idx, timeout=3)
            idxs = [json.loads(x) for x in r.text.split("\n") if x]
            r.close()
            df = pd.DataFrame(idxs)

            # Format the DataFrame
            df.index = df.index.rename("grib_message")
            df.index += 1
            df = df.reset_index()
            df["start_byte"] = df["_offset"]
            df["end_byte"] = df["_offset"] + df["_length"]
            df["range"] = df.start_byte.astype(str) + "-" + df.end_byte.astype(str)
            df["reference_time"] = pd.to_datetime(
                df.date + df.time, format="%Y%m%d%H%M"
            )
            df["step"] = pd.to_timedelta(df.step.astype(int), unit="h")
            df["valid_time"] = df.reference_time + df.step

            df = df.reindex(
                columns=[
                    "grib_message",
                    "start_byte",
                    "end_byte",
                    "range",
                    "reference_time",
                    "valid_time",
                    "step",
                    # ---- Used for search ------------------------------
                    "param",  # parameter field (variable)
                    "levelist",  # level
                    "levtype",  # sfc=surface, pl=pressure level, pt=potential vorticity
                    "number",  # model number (used in ensemble products)
                    "domain",  # g=global
                    "expver",  # experiment version
                    "class",  # classification (od=routing operations, rd=research, )
                    "type",  # fc=forecast, an=analysis,
                    "stream",  # oper=operationa, wave=wave, ef/enfo=ensemble,
                ]
            )

            df["search_this"] = (
                df.loc[:, "param":]
                .astype(str)
                .apply(
                    lambda x: ":" + ":".join(x).rstrip(":").replace(":nan:", ":"),
                    axis=1,
                )
            )

        # Attach some attributes
        df.attrs = dict(
            url=self.idx,
            source=self.idx_source,
            description="Inventory index file for the GRIB2 file.",
            model=self.model,
            product=self.product,
            lead_time=self.fxx,
            datetime=self.date,
        )

        return df
