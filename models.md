
# Forecast model documentation

Our forecasting approach is based on EMOS (ensemble model output
statistics), sometimes called forecast calibration or ensemble
postprocessing. This approach takes forecasts from a physics-based
weather model, and uses them as predictors in a statistical or machine
learning model, with true values (observations) as the outcome variable.
The trained models can then produce forecasts that remove the systematic
biases from the physics-based forecasts, improving the accuracy of the
forecast.

Ensembles add another layer of sophistication. Physics-based ensembles
produce a collection of plausible forecasts, with the goal of
representing the uncertainty of the forecast. The uncertainty can also
be calibrated, using an approach called *distributional regression*.
Unlike traditional regression, distributional regression allows both the
point prediction and the uncertainty of the prediction to vary based on
the predictors. This is useful for weather forecasting, where some types
of weather are much harder to forecast than others, so the forecast
uncertainty varies from day to day.

We’re following recent EMOS research by using machine learning
distributional regression models, using physics-based ensemble forecasts
from NOAA as predictors. We represent the uncertainty in both TV and
load models with a normal distribution.

## Training data (outcome variable)

Both TV and load models have been trained with daily peak loads, at the
network and NYC system level, from April through September, 2021-2025.

(Although we only use predictions starting in May, April is included so
the load model can learn the trends for a new year, before making the
official forecasts in May.)

## Predictors

Both TV and load models share the same weather predictors, which are
derived from the GEFS forecasts (an ensemble forecast) from NOAA.

We use the 06Z release of the GEFS forecasts since they are the most
recent forecasts available before 7:30am.

The file
<https://github.com/wmay/coned_load_forecasting/blob/main/config/gefs.csv>
describes the variables taken from GEFS. The choice of variables follows
Rasp and Lerch (2018).

For each variable except for TV, we calculate a 3-day weighted average
of the ensemble means, where the daily weights follow the definition of
TV (.1, .2, .7). Days are defined as starting and ending at 9pm New York
time, to align with the calculation for effective temperature (9am-9pm).
For each network, values from GEFS are interpolated using bilinear
interpolation to the centroid of the network area. To create the NYC
system predictors, we use a weighted average of the network predictors,
where the weights are the all-time peak loads of each network. For
forecasts less than 2 days ahead, the 3 days are not all included in the
current GEFS forecast, because they’ve already occurred. In that case,
we backfill those numbers with the GEFS analysis (where available) or
the ensemble mean of the shortest-term GEFS forecast.

We also include the ensemble mean and standard deviation of TV. This
requires calculating TV for each ensemble member (using temperature and
humidity data interpolated to the corresponding network centroid), and
then getting the mean and standard deviation. Note that this is
mathematically different than deriving TV from the ensemble mean
temperature and humidity. For forecasts less than 2 days ahead, where TV
has already been partially observed, we replace that portion of TV with
weather station observations.

## Objective function

For both TV and load, we use CRPS (continuous ranked probability score)
as the objective function. While log loss is more common in machine
learning, CRPS has been more common in weather research, in part because
it is more robust to outliers (Gebetsberger et al. 2018).

## TV models

The TV forecasts use random forest models similar to Taillardat and
Mestre (2020) with minor tweaks. Instead of quantile regression forests,
we use distributional random forests, as implemented in the R package
drf (Cevid et al. 2022). We include all lead times in the same model,
following Wessel, Ferro, and Kwasniok (2024). There are a total of 72
models– 70 individual network TV models, plus a model for the official
ConEd system TV, and one for the new UAlbany system TV.

Rather than predicting TV directly, we subtract the TV forecast from the
observed TV, so that we’re predicting the forecast error. This improves
the ability of the forests to extrapolate to extreme TV values.

On top of the weather predictors, we also include day of year (to
represent the season).

We suggest updating these models once a year.

## Load models

The load forecasts use a neural network based on Rasp and Lerch (2018).
All ConEd networks are included in a single model. The NYC system loads
have their own, separate model. We also include all lead times in the
same model, following Wessel, Ferro, and Kwasniok (2024). (So there are
a total of 2 load forecasting models.)

The neural network is composed of

- an embedding layer for ConEd networks
- one hidden layer
- an output layer with two nodes (one for the mean, one for standard
  deviation)

The neural network for the system load is the same, but without the
embedding layer.

On top of the weather predictors, we also include day of year (to
represent the season) and date (to allow for trends over time). Due to
the date predictor and the time trends it introduces, these models must
be updated regularly – we suggest once a month. In terms of date, all
training data is necessarily in the past. So every forecast of the
future is beyond the range of the training data. That means the models
should be updated relatively frequently, to avoid distortions caused by
extrapolating too far beyond the range of the training data.

The neural networks are fit using the torch package in R, using weight
decay for regularization.

## References

<div id="refs" class="references csl-bib-body hanging-indent"
entry-spacing="0">

<div id="ref-cevid_distributional_2022" class="csl-entry">

Cevid, Domagoj, Loris Michel, Jeffrey Näf, Peter Bühlmann, and Nicolai
Meinshausen. 2022. “Distributional Random Forests: Heterogeneity
Adjustment and Multivariate Distributional Regression.” *Journal of
Machine Learning Research* 23 (333): 1–79.
<http://jmlr.org/papers/v23/21-0585.html>.

</div>

<div id="ref-gebetsberger_estimation_2018" class="csl-entry">

Gebetsberger, Manuel, Jakob W. Messner, Georg J. Mayr, and Achim
Zeileis. 2018. “Estimation Methods for Nonhomogeneous Regression Models:
Minimum Continuous Ranked Probability Score Versus Maximum Likelihood.”
*Monthly Weather Review* 146 (12): 4323–38.
<https://doi.org/10.1175/MWR-D-17-0364.1>.

</div>

<div id="ref-rasp_neural_2018" class="csl-entry">

Rasp, Stephan, and Sebastian Lerch. 2018. “Neural Networks for
Postprocessing Ensemble Weather Forecasts.” *Monthly Weather Review* 146
(11): 3885–3900. <https://doi.org/10.1175/MWR-D-18-0187.1>.

</div>

<div id="ref-taillardat_research_2020" class="csl-entry">

Taillardat, Maxime, and Olivier Mestre. 2020. “From Research to
Applications – Examples of Operational Ensemble Post-Processing in
France Using Machine Learning.” *Nonlinear Processes in Geophysics* 27
(2): 329–47. <https://doi.org/10.5194/npg-27-329-2020>.

</div>

<div id="ref-wessel_lead-time-continuous_2024" class="csl-entry">

Wessel, Jakob Benjamin, Christopher A. T. Ferro, and Frank Kwasniok.
2024. “Lead-Time-Continuous Statistical Postprocessing of Ensemble
Weather Forecasts.” *Quarterly Journal of the Royal Meteorological
Society* 150 (761): 2147–67. <https://doi.org/10.1002/qj.4701>.

</div>

</div>
