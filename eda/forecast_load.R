# Evaluate daily peak load forecasting methods

# Issues to be addressed:

# - Feature engineering:
#   - how to summarize/aggregate hourly forecasts for daily predictions
#   - including lagged values
#   - summarizing/aggregating gridded values
# - Statistical issues:
#   - regularization due to small data size
#   - forecasting for multiple nearby areas (spatially autocorrelated errors)
#   - serially autocorrelated errors for longer forecast times

# A helpful trick: can subtract the lagged portion of TV from the output so it
# doesn't need to be predicted (and add it back afterward)

# NEXT STEPS:
# - get confidence intervals/p-values (randomize data etc.) May want to shuffle
#   minibatches. See https://arxiv.org/pdf/1206.5533, p. 6

setwd('..')
library(psychrolib)
library(timeDate) # holidays
library(magrittr)
library(stars) # also requires ncmeta
# install.packages('ncmeta')
# library(quantregForest)
# library(crch)
library(scoringRules)
library(beepr)
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

make_predictor_dataset = function(nc, lon, lat, ltime) {
  if (!is.null(attr(nc, 'dimensions')$refDate$values)) {
    refDates = attr(nc, 'dimensions')$refDate$values
  } else {
    refDates = with(attr(nc, 'dimensions')$refDate, {
      offset + (seq(from, to) - 1) * delta
    })
  }
  times = refDates + as.difftime(3 * ltime, units = 'hours')
  attr(times, 'tzone') = 'America/New_York'
  out = lapply(nc, function(x) x[lon, lat, ltime, ]) %>%
    as.data.frame
  cbind(time = times, out)
}

make_predictor_dataset2 = function(nc, init, lon, lat, days_ahead = 1) {
  # want 9am-9pm. 9 (3) is 4pm. 2-14, 3-12, 1-4
  # This should really be a function of the forecast init hour
  ltime = days_ahead * 8 + 1:4
  if (!is.null(attr(nc, 'dimensions')$refDate$values)) {
    refDates = attr(nc, 'dimensions')$refDate$values
  } else {
    refDates = with(attr(nc, 'dimensions')$refDate, {
      offset + (seq(from, to) - 1) * delta
    })
  }
  dates = as.Date(refDates + as.difftime(3 * ltime[1], units = 'hours'),
                  tz = 'America/New_York')
  # times = refDates + as.difftime(3 * ltime, units = 'hours')
  # attr(times, 'tzone') = 'America/New_York'
  init_ind = as.POSIXlt(refDates)$hour == init
  out = lapply(nc, function(x) colMeans(x[lon, lat, ltime, init_ind])) %>%
    as.data.frame
  # let's remove the sd columns, as they aren't very meaningful or useful
  sd_col = endsWith(names(out), '_sd')
  out = out[, !sd_col]
  # I also want max temperature
  tmp_max = apply(nc[['TMP']][lon, lat, ltime, init_ind], 2, max)
  # I also want "effective" temperature
  eff_tmp = with(nc, {
    tmp = TMP[lon, lat, ltime, init_ind]
    dpt = DPT[lon, lat, ltime, init_ind]
    time = refDates[init_ind] %>%
      sapply(function(x) x + as.difftime(3 * ltime, units = 'hours')) %>%
      as.POSIXct(tz = 'America/New_York')
    coned_effective_temp_3hr(c(time), as_fahrenheit(c(tmp) - 273.15),
                             as_fahrenheit(c(dpt) - 273.15))
  })
  cbind(day = dates[init_ind], out, TMP_max = tmp_max) %>%
    merge(eff_tmp, all.x = TRUE)
}

prepare_dataset = function(days_ahead = 1) {
  x = make_predictor_dataset2(gefs, 12, 3, 4, days_ahead = days_ahead)
  # add doy to account for seasonality
  x$doy = as.POSIXlt(x$day)$yday
  # y is the forecast error
  y = x %>%
    transform(obs_eff_temp = system_eff_tmp$eff_temp[match(day, system_eff_tmp$day)]) %>%
    transform(eff_tmp_fc_err = eff_temp - obs_eff_temp) %>%
    subset(select = eff_tmp_fc_err)
  ds = na.omit(cbind(x, y))
  x = subset(ds, select = -c(day, eff_tmp_fc_err))
  y = subset(ds, select = eff_tmp_fc_err)
  list(x = x, y = y, day = ds$day)
}

average_columns = function(dat, x) {
  .1 * dat[, paste0(x, '.1')] + .2 * dat[, paste0(x, '.2')] +
    .7 * dat[, paste0(x, '.3')]
}

average_3days = function(dat) {
  columns = names(dat)[endsWith(names(dat), '.1')]
  for (n in columns) {
    new_name = substr(n, 1, nchar(n) - 2)
    old_names = paste(new_name, 1:3, sep = '.')
    dat[, new_name] = average_columns(dat, new_name)
    dat[, old_names] = NULL
  }
  dat
}

# Just like for TV forecasting, but use GAM prediction error as the outcome and
# add GAM outputs as predictors. Also only use business days
prepare_multiday_dataset = function(days_ahead = 2, predict_error = TRUE) {
  # x: make a predictor dataset for each day, then combine
  d_list = lapply((days_ahead - 2):days_ahead, prepare_dataset)
  # "day" starts off as the valid day for each forecast. Change it to TV day 3,
  # that is, the last of the 3 days used to calculate TV
  for (i in 1:2) d_list[[i]]$day = d_list[[i]]$day + (3 - i)

  # # what was wrong with this?
  # d_list = lapply(days_ahead:(days_ahead - 2), prepare_dataset)
  # for (i in 2:3) d_list[[i]]$day = d_list[[i]]$day + (i - 1)
  
  # for (i in 1:3) d_list[[i]]$tv_day = i
  for (i in 1:3) d_list[[i]]$tv_day = (4 - i)
  v.names = with(d_list[[1]], c(names(x), names(y)))
  gam_loads = gam_fct$forecasts[, days_ahead - 1, ] %>%
    as.data.frame %>%
    cbind(day = gam_fct$model_day + days_ahead) %>%
    subset(isBizday(as.timeDate(day), holidays = holidayNYSE(2021:2024)))
  # gam_loads = data.frame(day = gam_fct$model_day + days_ahead)
  tv = make_nwp_tv_dataset(gefs_tv, 3, 4, days_ahead = days_ahead)
  out = d_list %>%
    lapply(function(x) with(x, cbind(x, y, tv_day, day))) %>%
    do.call(rbind, .) %>%
    reshape(direction = "wide", idvar = 'day', timevar = "tv_day",
            v.names = v.names) %>%
    # subset(select = -c(eff_temp.1, eff_temp.2, eff_temp.3)) %>%
    subset(select = -c(doy.1, doy.2, doy.3)) %>%
    transform(TV = tv$tv[match(day, tv$day)],
              TV_sd = tv$tv_sd[match(day, tv$day)],
              load = gam_loads$fit[match(day, gam_loads$day)],
              load_se = gam_loads$se.fit[match(day, gam_loads$day)],
              load_err = gam_loads$err[match(day, gam_loads$day)],
              load_obs = gam_loads$obs[match(day, gam_loads$day)]) %>%
    # transform(tv_fc_err = TV - system_tv$tv[match(day, system_tv$day)]) %>%
    subset(select = -c(eff_tmp_fc_err.1, eff_tmp_fc_err.2, eff_tmp_fc_err.3)) %>%
    na.omit %>%
    average_3days %>%
    subset(select = -eff_temp)

  if (predict_error) {
    x = subset(out, select = -c(day, load_err, load_obs))
    y = subset(out, select = load_err)
  } else {
    # no GAM prediction in this case
    x = subset(out, select = -c(day, load, load_se, load_err, load_obs))
    y = subset(out, select = load_obs)
  }
  x$doy = as.POSIXlt(out$day)$yday
  list(x = x, y = y, day = out$day)
  # # y: get either TV, or the future portion of TV
  # tv = system_tv[match(out$day, system_tv$day), 'tv']
  # list(x = out, y = tv)
}

# strangely `read_mdim` works fine while `read_ncdf` messes up the coordinates
# gefs_tv = read_ncdf('results/process_nwp_data/gefs_tv2.nc')
gefs_tv = read_mdim('results/process_nwp_data/gefs_tv2.nc')

# I don't know why this is such a pain to read
# wrong dimensions:
# read_ncdf('results/process_nwp_data/gefs_0p25_3day_wmean.nc')
# fails, but differently:
# read_mdim('results/process_nwp_data/gefs_0p25_3day_wmean.nc', variable = "?")
# works, but many warnings:
# read_stars('results/process_nwp_data/gefs_0p25_3day_wmean.nc')
# this will fail if some variables have dimensions in a different order!!
gefs_0p25_3day = read_mdim('results/process_nwp_data/gefs_0p25_3day_wmean.nc')


dim(gefs_0p25_3day['u10', 1, 1, , 1][['u10']])

# st_extract for interpolation


# set up the effective temp data
data_wide = readRDS('results/load_vs_weather/tv_and_load.rds')

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
# system_tv_sites = readRDS('results/select_stations/system_results.rds')$stids
# system_tv_sites = c("BKNYRD", "BKMAPL", "QNDKIL", "JFK", "SIFKIL", "LGA",
#                     "QNSOZO", "MHMHIL")
system_tv_sites = c('NYC', 'LGA')
system_tv = get_combined_tv(system_tv_sites)
# only 400 of these, clearly I should collect more weather data for this step.
# Note that this data is only May-September (technically April 29th)

gefs = read_ncdf('data/gefs_nick.nc')

# predict daily value (effective temp) instead of TV
system_eff_tmp = get_combined_eff_tmp(system_tv_sites)
# only 400 of these, clearly I should collect more weather data for this step

# out = list(model_day = forecast_days, gam_forecasts)
gam_fct = readRDS('results/load_curve_forecast/gam_forecasts.rds')

loads = read.csv('data/coned/Borough and System Data 2020-2024.csv') %>%
  transform(DT = as.POSIXct(DT, tz = 'EST5EDT', '%m/%d/%Y %H:%M'),
            # some inconsistency, but various forms of 'False' mean data is fine
            BAD = !startsWith(BAD, 'F')) %>%
  subset(!is.na(DT) & !BAD)
names(loads) = tolower(names(loads))

peaks = loads %>%
  subset(borough == 'CECONY') %>%
  transform(day = as.Date(dt, tz = 'EST5EDT')) %>%
  aggregate(reading ~ day, ., max) %>%
  subset(isBizday(as.timeDate(day), holidays = holidayNYSE(2021:2024)))

# load_curve_dat = merge(peaks, system_tv) %>%
#   subset(isBizday(as.timeDate(day), holidays = holidayNYSE(2021:2024))) %>%
#   transform(nday = as.integer(day), doy = as.POSIXlt(day)$yday)


# # get baseline accuracies-- must make sure to use exact same testing data!!
# nwp_acc = sapply(2:7, function(i) {
#   tv_fc_err = prepare_multiday_dataset(i)$y$tv_fc_err
#   n_train = ceiling(length(tv_fc_err) * 4 / 5)
#   test_data = tail(tv_fc_err, -n_train)
#   c(day = i, MAE = mean(abs(test_data)), RMSE = sqrt(mean(test_data^2)))
# }) %>%
#   t %>%
#   as.data.frame

# old (calculated from mean TMP and DWP)
#   day      MAE     RMSE
# 1   2 1.536112 2.102919
# 2   3 1.851774 2.469837
# 3   4 2.270325 2.830988
# 4   5 2.624355 3.262561
# 5   6 3.079295 3.875903
# 6   7 3.533117 4.485339

# new (mean TV from ens. members) -- it's better!
#   day      MAE     RMSE
# 1   2 1.374344 1.857518
# 2   3 1.635118 2.104380
# 3   4 2.032850 2.503328
# 4   5 2.392062 2.923228
# 5   6 2.829654 3.531663
# 6   7 3.236689 4.133973


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
# lgr::get_logger("mlr3")$set_threshold("trace")
# lgr::get_logger("bbotk")$set_threshold("trace")


make_tasks = function() {
  sapply(2:7, function(ahead) {
    id = paste0('Forecast TV ', ahead)
    prepare_multiday_dataset(ahead) %>%
      with(cbind(x, y[, 1, drop = FALSE])) %>%
      as_task_regr('load_err', id = id)
  })
}

# make_gefs_tasks = function() {
#   sapply(2:7, function(ahead) {
#     id = paste0('Forecast TV ', ahead)
#     gefs = get_gefs_samples(gefs_tv, days_ahead = ahead)
#     # make sure days are exactly the same as the prediction tasks
#     compare_days = prepare_multiday_dataset(ahead)$day
#     x = gefs$tv[match(compare_days, gefs$day), ]
#     y = system_tv[match(compare_days, system_tv$day), 'tv', drop = FALSE]
#     out = cbind(as.data.frame(x), y) %>%
#       as_task_regr('tv', id = id)
#   })
# }

# gefs_tasks = make_gefs_tasks()
forecast_load_tasks = make_tasks()


# I think the plan is that we run this benchmark over all the models

ngr = LearnerRegrNGR$new()

rangerdist = LearnerRegrRangerDist$new()
# add ranger default search space
rangerdist$param_set$values =
  mlr3misc::insert_named(rangerdist$param_set$values, lts('regr.ranger.default')$values)
rangerdist_at = auto_tuner(
    tuner = tnr("random_search"),
    learner = rangerdist, 
    resampling = rsmp('block_cv'),
    measure = msr('regr.crps'),
    term_evals = 10
)

drf = LearnerRegrDrf$new()
# drf hypertuning should be similar to ranger and grf. See
# `lts('regr.ranger.default')` and
# https://grf-labs.github.io/grf/REFERENCE.html#parameter-tuning
# drf$param_set$set_values(num.trees = to_tune(1, 2000),
#                          #sample.fraction = .8,
#                          sample.fraction = to_tune(0.1, 1),
#                          mtry.ratio = to_tune(0, 1),
#                          honesty.fraction = to_tune(.5, .8),
#                          honesty.prune.leaves = to_tune(),
#                          ci.group.size = 1,
#                          num.threads = 1)
# honesty is causing the code to crash somehow. But results may be better
# without it anyway. Probably should retry it if the data changes
drf$param_set$set_values(num.trees = to_tune(1, 2000),
                         sample.fraction = to_tune(0.1, 1),
                         mtry.ratio = to_tune(0, 1),
                         honesty = FALSE,
                         # honesty.fraction = to_tune(.5, .8),
                         # honesty.prune.leaves = to_tune(),
                         ci.group.size = 1,
                         num.threads = 1)
drf_at = auto_tuner(
    tuner = tnr("random_search"),
    learner = drf, 
    resampling = rsmp('block_cv'),
    measure = msr('regr.crps'),
    term_evals = 10,
    store_tuning_instance = TRUE
)
# system.time(drf_at$train(forecast_tv_tasks[[6]]))
# from the source code:
all.tunable.params <- c("sample.fraction", "mtry", "min.node.size", "honesty.fraction",
                          "honesty.prune.leaves", "alpha", "imbalance.penalty")



future::plan('multicore', workers = 4)
# future::plan("multisession", workers = 4)

# run the inner loop in parallel and the outer loop sequentially
# future::plan(list("sequential", "multisession"), workers = 4)
future::plan(list("sequential", "multicore"), workers = 4)
# future::plan("sequential")

bgrid = benchmark_grid(
    tasks = forecast_load_tasks,
    # learners = ngr,
    learners = c(ngr, rangerdist_at, drf_at),
    resamplings = rsmp('forecast_holdout')
)
system.time(bres <- progressr::with_progress(benchmark(bgrid)))
bres$aggregate(c(msr('regr.crps'), msr('regr.mape'), msr('regr.mae')))
saveRDS(bres, 'benchmarks.rds')


# now we're going to try without using the GAM
make_tasks2 = function() {
  sapply(2:7, function(ahead) {
    id = paste0('Forecast TV ', ahead)
    prepare_multiday_dataset(ahead, predict_error = FALSE) %>%
      with(cbind(x, y[, 1, drop = FALSE])) %>%
      as_task_regr('load_obs', id = id)
  })
}

# gefs_tasks = make_gefs_tasks()
forecast_load_tasks2 = make_tasks2()

bgrid = benchmark_grid(
    tasks = forecast_load_tasks2,
    # learners = ngr,
    # learners = drf_at,
    learners = c(ngr, rangerdist_at, drf_at),
    resamplings = rsmp('forecast_holdout')
)
system.time(bres2 <- progressr::with_progress(benchmark(bgrid, store_models = TRUE)))
bres2$aggregate(c(msr('regr.crps'), msr('regr.mape'), msr('regr.mae')))

# for (i in 6:1) {
#   print(i)
#   drf_at$train(forecast_load_tasks2[[i]])
# }

# drf_at$train(forecast_load_tasks2[[6]])






# get hyperparameter tuning results
extract_inner_tuning_results(bres2)

# slightly different for GEFS benchmarking
gefs_tester = LearnerRegrIdent$new()
gefs_tester$id = 'GEFS'
bgrid = benchmark_grid(
    tasks = gefs_tasks,
    learners = gefs_tester,
    resamplings = rsmp('forecast_holdout')
)
bres_gefs = benchmark(bgrid)
bres_gefs$aggregate(c(msr('regr.crps'), msr('regr.mae'), msr('regr.rmse')))


all_res = rbind(
    bres$aggregate(c(msr('regr.crps'), msr('regr.mae'), msr('regr.rmse'))),
    bres_gefs$aggregate(c(msr('regr.crps'), msr('regr.mae'), msr('regr.rmse')))
)
saveRDS(all_res, 'benchmarks.rds')

all_res %>%
  subset(select = c(task_id, learner_id, regr.crps)) %>%
  reshape(direction = 'wide', idvar = 'task_id', timevar = 'learner_id')





library(ggplot2)
# library(tidyr)

gefs_baseline = all_res %>%
  subset(learner_id == 'GEFS') %>%
  transform(days_ahead = as.integer(substr(task_id, 13, 13))) %>%
  with(setNames(regr.crps, days_ahead))

all_res %>%
  transform(days_ahead = as.integer(substr(task_id, 13, 13))) %>%
  ggplot(aes(x = days_ahead, y = regr.crps, color = learner_id)) +
  geom_line()

all_res %>%
  transform(days_ahead = as.integer(substr(task_id, 13, 13))) %>%
  transform(crps = regr.crps / (gefs_baseline[as.character(days_ahead)])) %>%
  ggplot(aes(x = days_ahead, y = crps, color = learner_id)) +
  geom_line()

comparison_plot = function(metrics, nwp_acc) {
  nwp_acc_long = nwp_acc %>%
    transform(days_ahead = day, model = 'GEFS') %>%
    subset(select = -day) %>%
    gather('metric', 'value', -c(model, days_ahead))
  names(metrics) = sub('regr\\.', '', names(metrics))
  names(metrics)[2:3] = toupper(names(metrics)[2:3])
  metrics_long = metrics %>%
    transform(model = 'Random forest') %>%
    gather('metric', 'value', -c(model, days_ahead))

  rbind(nwp_acc_long, metrics_long) %>%
    ggplot(aes(x = days_ahead, y = value, color = model)) + geom_line() + facet_wrap( ~ metric)
}

comparison_plot(rf_all_results$metrics, nwp_acc)
