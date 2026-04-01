
# About the app

TV and daily peak load forecasts come from machine learning models that use
weather forecasts from NOAA (GEFS). Forecasts are updated daily around 7:30am.
At the network level, we are forecasting custom versions of TV tailored to the
network, not the official ConEd TV.

The load forecasts are only valid on business days from May through September.
TV forecasts are valid on all days from May through September.

These forecasts are part of a project with UAlbany evaluating the use of the
[NYC Micronet](https://www.nysmesonet.org/networks/nyc). Earlier in the project,
we created a unique version of TV for each ConEd network, optimized to match the
network loads. The new TVs are derived from weather stations in the NYS Mesonet,
NYC Micronet, and ASOS weather networks (as opposed to the official ConEd TV,
which is derived only from Central Park and LaGuardia stations). We also created
an optimized TV for the NYC system as a whole.

These optimized TVs provide better measurements of local weather patterns,
notably the sea breeze, which has the largest impact on the southeastern coast
of the Long Island portion of NYC.

A previous summary of the project is available at
<https://nyswrcc.org/coeweather/coned/documents/ualbany_micronet_forecasting_2024.pdf>.

The forecast data files can be downloaded from
<https://nyswrcc.org/coeweather/coned/forecasts>.

Please contact William May at <wmay@hey.com> with bug reports about the app.
