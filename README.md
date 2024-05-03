# Con Edison/NYC weather and load forecasts

Code to analyze and forecast NYC weather and electrical loads. The weather
analysis focuses on the NYC sea breeze due to its affect on daily peak
electrical loads. The load forecasts are constrained to depend on ConEd's
"temperature variable" (TV), because ConEd relies on that number for
decision-making.

The overall forecasting strategy is

- select weather stations used for constructing the TV values for ConEd networks
- predict the TV values with a model output statistics approach
- predict loads based on the TV predictions

An R markdown report summarizes the results, and an R shiny app displays the
current forecasts.

## Dependencies

The required data from ConEd and the NYS Mesonet is not included here.

The maps require API keys for [Mapbox](https://www.mapbox.com/) and [Stadia
Maps](https://stadiamaps.com/).

The easiest way to install the R package dependencies is with the renv package:

```R
install.packages('renv')
pkg_info = renv::dependencies() # run this from the project root folder
pkg_dep = unique(pkg_info$Package)
install.packages(pkg_dep)
```

## Run the code

The entire analysis can be run with bash:

```sh
mapbox_key=your_key stadia_key=your_key ./run_analysis.sh
```
