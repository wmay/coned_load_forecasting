'''Get path info from Herbie.
'''

from datetime import timedelta
from typing import Optional, Union
import pandas as pd
from herbie import Path, config
from herbie.help import _search_help
from herbie.core import Datetime
import herbie.models as model_templates
from herbie import Herbie

class NwpPath(Herbie):
    '''Get file paths based on Herbie model templates. Takes the same arguments
    as `Herbie`. Unlike `Herbie`, this does not check for source availability.
    Instead the first source is used.
    '''
    # this is the original init method from Herbie, only without the source
    # check
    def __init__(
        self,
        date: Optional[Datetime] = None,
        *,
        valid_date: Optional[Datetime] = None,
        model: str = config["default"].get("model"),
        fxx: int = config["default"].get("fxx"),
        product: str = config["default"].get("product"),
        priority: Union[str, list[str]] = config["default"].get("priority"),
        save_dir: Union[Path, str] = config["default"].get("save_dir"),
        overwrite: bool = config["default"].get("overwrite", False),
        verbose: bool = config["default"].get("verbose", True),
        **kwargs,
    ):
        """Specify model output and find GRIB2 file at one of the sources."""
        self.fxx = fxx

        if isinstance(self.fxx, (str, pd.Timedelta)):
            # Convert pandas-parsable timedelta string to int in hours.
            self.fxx = pd.to_timedelta(fxx).round("1h").total_seconds() / 60 / 60
            self.fxx = int(self.fxx)

        if date:
            # User supplied `date`, which is the model initialization datetime.
            self.date = pd.to_datetime(date)
            self.valid_date = self.date + timedelta(hours=self.fxx)
        elif valid_date:
            # User supplied `valid_date`, which is the model valid datetime.
            self.valid_date = pd.to_datetime(valid_date)
            self.date = self.valid_date - timedelta(hours=self.fxx)
        else:
            raise ValueError("Must specify either `date` or `valid_date`")

        self.model = model.lower()
        self.product = product

        if self.model == "ecmwf":
            self.model = "ifs"
            warnings.warn(
                "`model='ecmwf'`is deprecated. Please use model='ifs' instead. Also, did you know you can also access `model='aifs'` too!",
                DeprecationWarning,
                stacklevel=2,
            )

        self.priority = priority
        self.save_dir = Path(save_dir).expand()
        self.overwrite = overwrite
        self.verbose = verbose

        # Some model templates may require kwargs not listed (e.g., `nest=`, `member=`).
        for key, value in kwargs.items():
            # TODO: Check if the kwarg is a config default.
            # TODO: e.g. if a user primarily works with RRFS, they may
            # TODO: want to configure "member" as a default argument.
            # You may also set IDX_SUFFIX as an argument.
            setattr(self, key, value)

        # Get details from the template of the specified model.
        # This attaches the details from the `models.<model>.template`
        # class to this Herbie object.
        # This line is equivalent to `model_templates.gfs.template(self)`.
        # I do it this way because the model name is a variable.
        # (see https://stackoverflow.com/a/7936588/2383070 for what I'm doing here)
        getattr(model_templates, self.model).template(self)

        if product is None:
            # The user didn't specify a product, so let's use the first
            # product in the model template.
            self.product = list(self.PRODUCTS)[0]
            log.info(f'`product` not specified. Will use "{self.product}".')
            # We need to rerun this so the sources have the new product value.
            getattr(model_templates, self.model).template(self)

        self.product_description = self.PRODUCTS[self.product]

        # Specify the suffix for the inventory index files.
        # Default value is `.grib2.idx`, but some have weird suffix,
        # like archived RAP on NCEI are `.grb2.inv`.
        self.IDX_SUFFIX = getattr(self, "IDX_SUFFIX", [".grib2.idx"])

        # Specify the index file type. By default, Herbie assumes the
        # index file was created with wgrib2.
        # But for ecmwf files with index files created with eccodes
        # the index files are in a different style.
        self.IDX_STYLE = getattr(self, "IDX_STYLE", "wgrib2")

        self.search_help = _search_help(self.IDX_STYLE)

        # Check the user input
        self._validate()

        # # Ok, now we are ready to look for the GRIB2 file at each of the remote sources.
        # # self.grib is the first existing GRIB2 file discovered.
        # # self.idx is the first existing index file discovered.
        # self.grib, self.grib_source = self.find_grib()
        # self.idx, self.idx_source = self.find_idx()
        # Just use the first souce (AWS)
        self.grib_source = list(self.SOURCES.keys())[0]
        self.grib = self.SOURCES[self.grib_source]
        self.idx_source = list(self.SOURCES.keys())[0]
        self.idx = self.grib + '.idx'

        # if verbose:
        #     # ANSI colors added for style points
        #     if any([self.grib is not None, self.idx is not None]):
        #         print(
        #             "✅ Found",
        #             f"┊ model={self.model}",
        #             f"┊ {ANSI.italic}product={self.product}{ANSI.reset}",
        #             f"┊ {ANSI.green}{self.date:%Y-%b-%d %H:%M UTC}{ANSI.bright_green} F{self.fxx:02d}{ANSI.reset}",
        #             f"┊ {ANSI.orange}{ANSI.italic}GRIB2 @ {self.grib_source}{ANSI.reset}",
        #             f"┊ {ANSI.orange}{ANSI.italic}IDX @ {self.idx_source}{ANSI.reset}",
        #         )
        #     else:
        #         print(
        #             "💔 Did not find",
        #             f"┊ model={self.model}",
        #             f"┊ {ANSI.italic}product={self.product}{ANSI.reset}",
        #             f"┊ {ANSI.green}{self.date:%Y-%b-%d %H:%M UTC}{ANSI.bright_green} F{self.fxx:02d}{ANSI.reset}",
        #         )
