# can I do distributional regression with RF-GLS?

setwd('..')
library(psychrolib)
library(magrittr)
library(stars) # also requires ncmeta
# install.packages('ncmeta')
source('R/load_data.R')
source('R/coned_tv.R')
source('R/mlr3_additions.R') # block CV, mean pinball loss, etc.
source('R/mlr3_distr.R')
source('R/forecast_dataset.R')

get_combined_eff_tmp = function(stids) {
  tmp_cols = paste0('tmpf.', stids)
  dwp_cols = paste0('dwpf.', stids)
  tmps = rowMeans(station_obs[, tmp_cols, drop = FALSE], na.rm = TRUE)
  dwps = rowMeans(station_obs[, dwp_cols, drop = FALSE], na.rm = TRUE)
  out = coned_effective_temp(station_obs$time, tmps, dwps)
  row.names(out) = as.character(out$day)
  out
}

# Get combined TV by averaging site dry+wet temperatures first and *then*
# getting the daily maximum of the averaged values. This ensures that the values
# are all taken from the same time of day
# IMPORTANT: Can't remove holidays prior to this calculation bc they're needed for lagged values
get_combined_tv = function(stids) {
  tmp_cols = paste0('tmpf.', stids)
  dwp_cols = paste0('dwpf.', stids)
  tmps = rowMeans(station_obs[, tmp_cols, drop = FALSE], na.rm = TRUE)
  dwps = rowMeans(station_obs[, dwp_cols, drop = FALSE], na.rm = TRUE)
  out = coned_tv(station_obs$time, tmps, dwps)
  row.names(out) = as.character(out$day)
  out
}

# calculate effective temp from 3-hourly forecast data (requires fahrenheit)
coned_effective_temp_3hr = function(t, temp, dwp) {
  data.frame(time = t, temp = temp, dwp = dwp) %>%
    subset(!is.na(temp) & !is.na(dwp)) %>%
    transform(wetbulb = coned_wet_bulb(temp, dwp)) %>%
    transform(eff_temp = (temp + wetbulb) / 2) %>%
    # remove night/morning
    transform(hour_of_day = as.integer(format(time, '%H'))) %>%
    subset(hour_of_day >= 9 & hour_of_day <= 21) %>%
    # get max by day
    transform(day = as.Date(time, tz = 'EST5EDT')) %>%
    aggregate(eff_temp ~ day, ., max, na.rm = T)
}

# ConEd's "temperature variable". `temp` and `dwp` must be in Fahrenheit
coned_tv_nwp = function(t, temp, dwp) {
  coned_effective_temp(t, temp, dwp) %>%
    transform(tv = rolling_mean3(day, eff_temp, 'days', c(.7, .2, .1))) %>%
    transform(tv = coned_round(tv, 1)) %>%
    subset(select = c(day, tv))
}

# strangely `read_mdim` works fine while `read_ncdf` messes up the coordinates
# gefs_tv = read_ncdf('results/process_nwp_data/gefs_tv2.nc')
# gefs_tv2 = read_mdim('results/process_nwp_data/gefs_tv2.nc')
# gefs_tv_members = read_mdim('results/process_nwp_data/gefs_tv_members.nc')

# I don't know why this is such a pain to read
# wrong dimensions:
# read_ncdf('results/process_nwp_data/gefs_0p25_3day_wmean.nc')
# fails, but differently:
# read_mdim('results/process_nwp_data/gefs_0p25_3day_wmean.nc', variable = "?")
# works, but many warnings:
# read_stars('results/process_nwp_data/gefs_0p25_3day_wmean.nc')
# this will fail if some variables have dimensions in a different order!!
read_nwp_nc = function(f) {
  out = read_mdim(f)
  # GEFS reports longitude in a non-standard way, must be corrected to work with
  # other tools
  attr(out, 'dimensions')$longitude$offset =
    attr(out, 'dimensions')$longitude$offset - 360
  out
}
gefs_0p25_3day = read_nwp_nc('results/process_nwp_data/gefs_0p25_3day_wmean.nc')
gefs_0p5_3day = read_nwp_nc('results/process_nwp_data/gefs_0p5_3day_wmean.nc')
gefs_tv = read_nwp_nc('results/process_nwp_data/gefs_tv.nc')

get_valid_day = function(nc, days_ahead) {
  # time is generally read by R as a date, due to the 1 day differences
  attr(nc, 'dimensions')$time$values + days_ahead
}

# get values from stars object at the centroid of a network
wgs84 = 4326
networks = readRDS('results/maps/coned_networks_cleaned.rds')
brownsville = networks[3, ] %>%
  st_centroid %>%
  st_transform(wgs84)

# note that the day returned here is the valid day
extract_values = function(nc, days_ahead) {
  points = merge(brownsville,
                 data.frame(time = get_valid_day(nc, 0)))
  out = nc %>%
    dplyr::slice('edt9pm_day', days_ahead + 1) %>%
    st_extract(points, bilinear = TRUE, time_column = 'time') %>%
    transform(day = time + days_ahead) %>%
    as.data.frame %>%
    subset(select = -c(time, Shape))
}

prepare_tv_task = function(days_ahead = 2, predict_error = TRUE) {
  # Note: netcdf files use forecast day, while `system_tv` uses valid day. Must
  # be consistent!
  out_0p25 = extract_values(gefs_0p25_3day, days_ahead) %>%
    subset(select = -eff_temp) # 3-day eff_temp is same as TV
  # sometimes soilw is all NA? for now just skip it if so
  if (all(is.na(out_0p25$soilw))) {
    warning('All soilw values are missing')
    out_0p25$soilw = NULL
  }
  out_0p5 = extract_values(gefs_0p5_3day, days_ahead) 
  out_tv = extract_values(gefs_tv, days_ahead) 
  # For days_ahead < 2, need to add observed (past) portion of TV to the
  # forecast portion. TV_sd is correct though, since observations have 0 sd
  if (days_ahead == 1) {
    out_tv$TV = out_tv$TV +
      .1 * system_tv$tv[match(out_tv$day - 2, system_tv$day)]
  }
  if (days_ahead == 0) {
    out_tv$TV = out_tv$TV +
      .2 * system_tv$tv[match(out_tv$day - 1, system_tv$day)] +
      .1 * system_tv$tv[match(out_tv$day - 2, system_tv$day)]
  }
  out = merge(out_0p25, out_0p5, by = 'day') %>%
    merge(out_tv, by = 'day') %>%
    transform(doy = as.POSIXlt(day)$yday)
  id = paste0('Forecast TV ', days_ahead)
  if (predict_error) {
    # Predict the GEFS forecast error instead of directly predicting TV
    out_dat = out %>%
      transform(fct_err = TV - system_tv$tv[match(day, system_tv$day)]) %>%
      na.omit
    dates = out_dat$day
    out = out_dat %>%
      subset(select = -day) %>%
      as_task_regr('fct_err', id = id)
  } else {
    out %>%
      transform(tv_obs = system_tv$tv[match(day, system_tv$day)]) %>%
      subset(select = -day) %>%
      na.omit %>%
      as_task_regr('tv_obs', id = id)
  }
  out$extra_args = list(dates = dates)
  out
}


SetUnitSystem('SI')
tv_vars = c('tair', 'relh')
nysm_obs = get_combined_nc_data(tv_vars, hour_to_nc_index(7:21)) %>%
  subset(!(is.na(tair) | is.na(relh))) %>%
  transform(dwp = GetTDewPointFromRelHum(tair, relh / 100)) %>%
  transform(tmpf = as_fahrenheit(tair), dwpf = as_fahrenheit(dwp)) %>%
  subset(select = c(time, stid, tmpf, dwpf)) %>%
  reshape(direction = 'wide', idvar = 'time', timevar = 'stid')

asos_obs = get_hourly_asos_data(7:21) %>%
  transform(stid = station, time = valid_hour) %>%
  subset(!(is.na(tmpf) | is.na(dwpf))) %>%
  subset(select = c(time, stid, tmpf, dwpf)) %>%
  reshape(direction = 'wide', idvar = 'time', timevar = 'stid')

station_obs = merge(nysm_obs, asos_obs, all = T)


# selected_sites = readRDS('results/select_stations/selected_sites.rds')
selected_sites = c('BKLN', 'QNSOZO', 'LDJ', 'QNHBEA', 'LGA', 'JFK', 'SIFKIL')
system_tv = get_combined_tv(selected_sites) %>%
  # not training on 2024 data
  subset(day < '2024-01-01')


library(RandomForestsGLS)

# forecast_tv_tasks = lapply(0:7, prepare_tv_task, predict_error = FALSE)
forecast_tv7 = prepare_tv_task(1)
task = forecast_tv7

y = as.matrix(task$data(cols = task$target_names))
X = as.matrix(task$data(cols = task$feature_names))

      # RFGLS_estimate_timeseries(y, X, Xtest = NULL, nrnodes = NULL,
      #                           nthsize = 20, mtry = 1,
      #                           pinv_choice = 1, n_omp = 1,
      #                           ntree = 50, h = 1, lag_params = 0.5,
      #                           variance = 1,
      #                           param_estimate = FALSE,
      #                           verbose = FALSE)

rf = RFGLS_estimate_timeseries(y, X, param_estimate = TRUE, verbose = TRUE)

# can we treat this as space to get more accurate results?
coords = cbind(as.integer(forecast_tv7$dates), 0)
rf2 = RFGLS_estimate_spatial(coords, y, X, param_estimate = TRUE, verbose = TRUE)

names(rf)
# [1] "P_matrix"         "predicted_matrix" "predicted"        "X"               
# [5] "y"                "RFGLS_object"

# P is the random selection (with replacement)

# let's say we want to get the nodes for X[1, ]

# cutoff decision is <= upper

get_leaf_members = function(fit, x, tree) {
  node_ind = 1
  s_ind = fit$P_matrix[, tree] + 1
  s = fit$X[s_ind,, drop = FALSE] # the sample
  with(fit$RFGLS_object, {
    while (ldaughter[node_ind, tree]) {
      node_mbest = mbest[node_ind, tree]
      node_upper = upper[node_ind, tree]
      x_l = x[node_mbest] <= node_upper
      tryCatch(s[, node_mbest], error = function(e) {
        print(head(s))
        print(dim(s))
        print(node_mbest)
      })
      s_l = s[, node_mbest] <= node_upper
      keep = if (x_l) s_l else !s_l
      s_ind = s_ind[keep]
      s = s[keep,, drop = FALSE]
      if (x_l) {
        node_ind = ldaughter[node_ind, tree]
      } else {
        node_ind = rdaughter[node_ind, tree]
      }
    }
    fit$y[s_ind, ]
  })
}

.predict_rfgls_distr = function(fit, newdata) {
  ntrees = ncol(fit$P_matrix)
  leaves = 1:ntrees %>%
    lapply(function(tree) get_leaf_members(fit, newdata, tree)) %>%
    unlist
  c(mean = mean(leaves), sd = sd(leaves))
}

predict_rfgls_distr = function(fit, newdata) {
  newdata %>%
    apply(1, function(x) .predict_rfgls_distr(fit, x)) %>%
    t
}

.predict_rfgls_distr(rf, newdata)






# test dataset
system_tv = get_combined_tv(selected_sites) %>%
  subset(day > '2024-01-01')

forecast_tv7_test = prepare_tv_task(1)
task_test = forecast_tv7_test$task


newdata = task_test$data(cols = task_test$feature_names)

preds = predict_rfgls_distr(rf, newdata)

# did it work?
test_fct_err = task_test$data(cols = task_test$target_names)$fct_err
plot(preds[, 'mean'], test_fct_err)

mean(abs(preds[, 'mean'] - test_fct_err))
# [1] 0.8242659
# ok that's reasonable


# cool, let's see how well it does

LearnerRegrRfgls = R6::R6Class(
  "LearnerRegrRfgls",
  inherit = LearnerRegr,
  public = list(
    initialize = function() {
      param_set = paradox::ps(
          nrnodes             = paradox::p_int(default = NULL, lower = 1L, special_vals = list(NULL), tags = "train"),
          nthsize             = paradox::p_int(default = NULL, lower = 1L, special_vals = list(NULL), tags = "train"),
          mtry                = paradox::p_int(lower = 1L, special_vals = list(NULL), tags = "train"),
          mtry.ratio          = paradox::p_dbl(lower = 0, upper = 1, tags = "train"),
          ntree               = paradox::p_int(1L, default = 500L, tags = c("train", "predict", "hotstart")),
          h                   = paradox::p_int(1L, default = 1L, tags = "train"),
          lags                = paradox::p_int(0L, default = 1L, tags = "train"),
          lag_params          = paradox::p_uty(default = 0.5, tags = "train"),
          param_estimate      = paradox::p_lgl(default = FALSE, tags = "train"),
          verbose             = paradox::p_lgl(default = FALSE, tags = c("train", "predict"))
      )
      # param_set$set_values(xval = 10L)
      super$initialize(
        id = "regr.rfgls",
        feature_types = c("logical", "integer", "numeric", "factor", "ordered"),
        predict_types = "distr",
        packages = "RandomForestsGLS",
        param_set = param_set,
        properties = c("missings"),
        label = "Non-homogeneous Gaussian Regression"
      )
    }
  ),
  private = list(
    .train = function(task) {
      pv = self$param_set$get_values(tags = "train")
      pv = mlr3learners:::convert_ratio(pv, "mtry", "mtry.ratio", length(task$feature_names))
      if (!is.null(pv$lags)) {
        pv$lag_params = rep(.5, pv$lags)
        pv$lags = NULL
      }
      mlr3misc::invoke(
          RandomForestsGLS::RFGLS_estimate_timeseries,
          y = as.matrix(task$data(cols = task$target_names)),
          X = as.matrix(task$data(cols = task$feature_names)),
          .args = pv
      )
      # RFGLS_estimate_spatial(coords, y, X, param_estimate = TRUE, verbose = TRUE)
      # coords = cbind(as.integer(task$extra_args$dates), 0)
      # mlr3misc::invoke(
      #     RandomForestsGLS::RFGLS_estimate_spatial,
      #     coords,
      #     y = as.matrix(task$data(cols = task$target_names)),
      #     X = as.matrix(task$data(cols = task$feature_names)),
      #     .args = pv
      # )
    },
    .predict = function(task) {
      pv = self$param_set$get_values(tags = "predict")
      newdata = task$data(cols = task$feature_names)
      pred = predict_rfgls_distr(self$model, newdata)
      # need to wrap in a `VectorDistribution` (see `LearnerRegr`)
      distrs = as.data.frame(pred) %>%
        distr6::VectorDistribution$new(distribution = "Normal", params = .)
      list(distr = distrs)
    }
  )
)

# try to simplify this using mlr3
# install.packages(c('ranger', 'drf', 'RandomForestsGLS'))
# remotes::install_github("mlr-org/mlr3temporal")
# distr6_repos =  c(CRAN = 'https://cloud.r-project.org',
#                   raphaels1 = 'https://raphaels1.r-universe.dev')
# install.packages('distr6', repos = distr6_repos)
library(mlr3)
library(mlr3learners)
library(mlr3tuning)
library(mlr3tuningspaces)
library(mlr3temporal)

# see notes in `?benchmark` and `?mlr_tuners_random_search`
lgr::get_logger("mlr3")$set_threshold("warn")
lgr::get_logger("bbotk")$set_threshold("warn")

rangerdist = LearnerRegrRangerDist$new()
# add ranger default search space
rangerdist$param_set$values =
  mlr3misc::insert_named(rangerdist$param_set$values, lts('regr.ranger.default')$values)
rangerdist_at = auto_tuner(
    tuner = tnr("random_search"),
    learner = rangerdist, 
    resampling = rsmp('block_cv', folds = 4),
    measure = msr('regr.crps'),
    term_evals = 6,
    store_tuning_instance = FALSE
)

rfgls = LearnerRegrRfgls$new()
# rfgls$param_set$set_values(param_estimate = TRUE, mtry.ratio = 1/3)
# rfgls$param_set$set_values(param_estimate = TRUE, ntree = 100, mtry.ratio = 1/3,
#                            nthsize = 16, h = 1, lags = 1)
rfgls$param_set$set_values(param_estimate = TRUE, mtry.ratio = to_tune(0, 1),
                           ntree = to_tune(1, 500), nthsize = to_tune(8, 20),
                           h = 1, lags = to_tune(1, 2))
rfgls_at = auto_tuner(
    tuner = tnr("random_search"),
    learner = rfgls, 
    resampling = rsmp('block_cv', folds = 4),
    measure = msr('regr.crps'),
    term_evals = 6,
    store_tuning_instance = FALSE
)

# back to the full dataset
system_tv = get_combined_tv(selected_sites) %>%
  subset(day < '2025-01-01')

future::plan(list('sequential', 'multicore'), workers = 4)
# future::plan(list(
#             tweak('multicore', workers = 2),
#             tweak('multicore', workers = I(3))
#         ))

forecast_tv_tasks = lapply(0:7, prepare_tv_task, predict_error = TRUE)
bgrid = benchmark_grid(
    tasks = forecast_tv_tasks[8],
    learners = c(rangerdist_at, rfgls_at),
    resamplings = rsmp('forecast_holdout')
)
system.time(bres <- progressr::with_progress(benchmark(bgrid)))
bres$aggregate(c(msr('regr.crps'), msr('regr.mae'), msr('regr.rmse')))
# no clear improvement sadly, but it's definitely competitive
