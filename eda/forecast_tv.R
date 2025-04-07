# Plan: This script will prepare the machine learning datasets, and then the
# models will be fit and evaluated in python




# run some simple models for baseline comparisons

# use ensemble MOS with neural networks to predict ConEd's TV

# Mostly following Rasp 2018 (10.1175/MWR-D-18-0187.1) and Lerch 2022
# (https://arxiv.org/abs/2204.05102). May also try some sort of RNN for the load
# forecasting part

# Results of this Masters thesis suggest two hidden layers may be better than 1:
# http://hdl.handle.net/2429/78331

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

# To reduce hourly values to daily predictors, can probably get max temperature
# (and other values at the same time?). Eh, I'll just pick a relevant time for
# now

# Fix tensorflow 2.16.1. This must be set in the shell/emacs before running R.
# See
# https://github.com/tensorflow/tensorflow/issues/63362#issuecomment-1988630226.

# LD_LIBRARY_PATH=/home/wmay/.local/lib/python3.10/site-packages/nvidia/cublas/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/cuda_cupti/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/cuda_nvcc/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/cuda_nvrtc/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/cuda_runtime/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/cudnn/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/cufft/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/curand/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/cusolver/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/cusparse/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/nccl/lib:/home/wmay/.local/lib/python3.10/site-packages/nvidia/nvjitlink/lib

# This is fixed in 2.17 but I can't install it because they dropped support for
# my GPU :'(

# set this if you don't want an annoying prompt from reticulate
# Sys.setenv(RETICULATE_PYTHON = '/usr/bin/python3')

# May want to shuffle minibatches. See https://arxiv.org/pdf/1206.5533, p. 6

# NEXT STEPS:
# - get confidence intervals/p-values (randomize data etc.)
# - Calculate effective temperature for each ensemble member for correct averaging and sd

# DONE:
# - fix predictors (correctly calculate daily mean, etc.)
# - test 7-day ahead forecasts


library(psychrolib)
library(magrittr)
library(stars) # also requires ncmeta
# install.packages('ncmeta')
library(quantregForest)
# library(crch)
library(scoringRules)
library(beepr)
source('../R/load_data.R')
source('../R/coned_tv.R')
# block CV, mean pinball loss, etc.
source('../R/mlr3_additions.R')

# get_effective_temp = function(x) {
#   drywet = with(x, dry_wet_ave(tair, relh / 100, pres * 100))
#   eff_temp = rolling_mean3(x$time, drywet, 'hours')
#   hour_of_day = as.integer(format(x$time, '%H'))
#   data.frame(stid = x$stid, time = x$time, eff_temp = eff_temp) %>%
#     subset(hour_of_day >= 9 & hour_of_day <= 21)
# }

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
# get_combined_tv = function(stids) {
#   cols = paste0('eff_temp.', stids)
#   combined_temps = rowMeans(et_wide[, cols, drop = FALSE])
#   out = data.frame(time = et_wide$time, eff_temp = combined_temps) %>%
#     # transform(day = as.Date(time, tz = 'EST5EDT')) %>%
#     # # should `na.rm` be true here? Not sure
#     # aggregate(eff_temp ~ day, ., max, na.rm = TRUE) %>%
#     transform(tv = rolling_mean3(day, eff_temp, 'days', c(.7, .2, .1))) %>%
#     subset(select = c(day, tv))
#   row.names(out) = as.character(out$day)
#   out
# }

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

# based on Rasp+Lerch 2018,
# https://github.com/slerch/ppnn/blob/7af9af3cfb754b8b6da54d2fe5d917cd54e32b9d/nn_postprocessing/nn_src/losses.py#L14
crps_loss = function(y_true, y_pred) {
  # Split y_pred into mean and std dev
  mu = y_pred[, 1]
  sigma = op_abs(y_pred[, 2])
  # Calculate the standardized difference
  z = (y_true[, 1] - mu) / sigma
  # Calculate the CDF and PDF of the standard normal distribution
  pdf = op_exp(-op_square(z) / 2) / sqrt(2 * pi)
  cdf = 0.5 * (1 + op_erf(z / sqrt(2)))
  # Calculate the CRPS components
  term1 = z * (2 * cdf - 1)
  term2 = 2 * pdf - 1 / sqrt(pi)
  crps = sigma * (term1 + term2)
  return(op_mean(crps))
}

# Based on Rasp+Lerch 2018:
# https://github.com/slerch/ppnn/blob/7af9af3cfb754b8b6da54d2fe5d917cd54e32b9d/nn_postprocessing/nn_src/losses.py#L14
# however here I'm using a log link for sigma
crps_loss2 = function(y_true, y_pred) {
  # Split y_pred into mean and std dev
  mu = y_pred[, 1]
  sigma = op_exp(y_pred[, 2])
  # Calculate the standardized difference
  z = (y_true[, 1] - mu) / sigma
  # Calculate the CDF and PDF of the standard normal distribution
  pdf = op_exp(-op_square(z) / 2) / sqrt(2 * pi)
  cdf = 0.5 * (1 + op_erf(z / sqrt(2)))
  # Calculate the CRPS components
  term1 = z * (2 * cdf - 1)
  term2 = 2 * pdf - 1 / sqrt(pi)
  crps = sigma * (term1 + term2)
  return(op_mean(crps))
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

prepare_multiday_dataset = function(days_ahead = 2) {
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
  out = d_list %>%
    lapply(function(x) with(x, cbind(x, y, tv_day, day))) %>%
    do.call(rbind, .) %>%
    reshape(direction = "wide", idvar = 'day', timevar = "tv_day",
            v.names = v.names) %>%
    # subset(select = -c(eff_temp.1, eff_temp.2, eff_temp.3)) %>%
    subset(select = -c(doy.1, doy.2, doy.3)) %>%
    na.omit %>%
    transform(TV = .7 * eff_temp.1 + .2 * eff_temp.2 + .1 * eff_temp.1) %>%
    # find TV forecast error
    transform(tv_fc_err = .7 * eff_tmp_fc_err.1 + .2 * eff_tmp_fc_err.2 +
              .1 * eff_tmp_fc_err.1) %>%
    subset(select = -c(eff_tmp_fc_err.1, eff_tmp_fc_err.2, eff_tmp_fc_err.3))
  x = subset(out, select = -c(day, tv_fc_err))
  x$doy = as.POSIXlt(out$day)$yday
  y = subset(out, select = tv_fc_err)
  list(x = x, y = y, day = out$day)
  # # y: get either TV, or the future portion of TV
  # tv = system_tv[match(out$day, system_tv$day), 'tv']
  # list(x = out, y = tv)
}


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


# get baseline accuracies-- must make sure to use exact same testing data!!
nwp_acc = sapply(2:7, function(i) {
  tv_fc_err = prepare_multiday_dataset(i)$y$tv_fc_err
  n_train = ceiling(length(tv_fc_err) * 4 / 5)
  test_data = tail(tv_fc_err, -n_train)
  c(day = i, MAE = mean(abs(test_data)), RMSE = sqrt(mean(test_data^2)))
}) %>%
  t %>%
  as.data.frame


# try to simplify this using mlr3
library(mlr3)
library(mlr3learners)
# install.packages('ranger')
library(mlr3tuning)
library(mlr3tuningspaces)
# remotes::install_github("mlr-org/mlr3temporal")
library(mlr3temporal)

# see notes in `?benchmark` and `?mlr_tuners_random_search`
lgr::get_logger("mlr3")$set_threshold("warn")
lgr::get_logger("bbotk")$set_threshold("warn")

# quantile regression forest
qrf = lrn('regr.ranger', predict_type = 'quantiles', quantiles = c(.1, .5, .9),
          quantile_response = .5)
# see docs for `LearnerRegr` for some of the quantile settings

# add an autotuner for qrf
qrf_at = auto_tuner(
    tuner = tnr("random_search"),
    learner = lts(qrf), # default search space
    resampling = rsmp('block_cv'),
    measure = msr('regr.mean_pinball'),
    term_evals = 10
)

forecast_tv_ahead2 = prepare_multiday_dataset(2) %>%
  with(cbind(x, y[, 1, drop = FALSE])) %>%
  as_task_regr(target = 'tv_fc_err', id = 'Forecast TV 2')

make_tasks = function() {
  sapply(2:7, function(ahead) {
    id = paste0('Forecast TV ', ahead)
    prepare_multiday_dataset(ahead) %>%
      with(cbind(x, y[, 1, drop = FALSE])) %>%
      as_task_regr('tv_fc_err', id = id)
  })
}

forecast_tv_tasks = make_tasks()

# quick benchmark
bgrid = benchmark_grid(
    tasks = forecast_tv_tasks[1:3],
    learners = qrf_at,
    resamplings = rsmp('forecast_holdout')
)

future::plan('multicore', workers = 3)

bres = benchmark(bgrid) # it works!

bres$aggregate(c(msr('regr.mean_pinball'), msr('regr.rmse')))




# # train the model
# qrf$train(forecast_tv, split$train_set)
# prediction = qrf$predict(forecast_tv, split$test_set)



# # hypertune
# qrf_instance = tune(
#     tuner = tnr("random_search"),
#     task = forecast_tv,
#     learner = lts(qrf), # default search space
#     resampling = rsmp("forecast_holdout"),
#     measures = msr("regr.mean_pinball"),
#     term_evals = 10
# )
# # can also use `ti` and `tnr`

# library(mlr3viz)
# autoplot(qrf_instance)

eval_model3 = function(eval_data, learner = qrf) {
  forecast_tv_all = with(eval_data, cbind(x, y[, 1, drop = FALSE]))
  # tune model on training data
  n_train = ceiling(nrow(forecast_tv_all) * 4 / 5)
  forecast_tv_train = forecast_tv_all %>%
    head(n_train) %>%
    as_task_regr(tv_fc_err ~ ., .)
  objective = ifelse(is.null(learner$quantiles), 'regr.mse',
                     'regr.mean_pinball')
  tune_res = tune(
      tuner = tnr('random_search'),
      task = forecast_tv_train,
      learner = lts(learner), # default search space
      resampling = rsmp('forecast_holdout'),
      measures = msr(objective),
      term_evals = 10
  )
  # now that we've tuned it, evaluate on unseen data
  learner$param_set$values = tune_res$result_learner_param_vals
  forecast_tv_eval = forecast_tv_all %>%
    tail(-n_train) %>%
    as_task_regr(tv_fc_err ~ ., .)
  learner$train(forecast_tv_train)
  prediction = learner$predict(forecast_tv_eval)
  # get pinball and MAE
  measures = c('regr.mae', 'regr.rmse')
  if (!is.null(learner$quantiles)) measures[3] = 'regr.mean_pinball'
  metrics = sapply(measures, function(m) unname(prediction$score(msr(m))))
  # pinball = prediction$score(msr('regr.mean_pinball'))
  # mae = prediction$score(msr('regr.mae'))
  # metrics = c(Mean_Pinball = unname(pinball), MAE = unname(mae))
  list(metrics = metrics)
}


eval_all_days3 = function(learner = qrf) {
  res_list = lapply(2:7, function(ahead) {
    message('starting ', ahead)
    eval_data = prepare_multiday_dataset(ahead)
    res = eval_model3(eval_data, learner)
    res$day = eval_data$day
    message('finished ', ahead)
    print(res$metrics)
    res
  })
  daily_metrics = sapply(res_list, getElement, name='metrics') %>%
    t %>%
    as.data.frame %>%
    transform(days_ahead = 2:7)
  list(metrics = daily_metrics[, 3:1], results = res_list)
}


all_results = eval_all_days3()

rand_for = lrn('regr.ranger')

rf_all_results = eval_all_days3(rand_for)

rf_all_results$metrics

nwp_acc$model = 'GEFS'

cbind(rf_all_results$metrics, nwp_acc[, -1])[, c(1, 3:4, 2, 5)]

library(ggplot2)
library(tidyr)

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


# let's try out torch and sequence-to-vector RNN
library(mlr3torch)
library(torch)

# create an RNN model using the regr.torch_model class

seq2vec_rnn = nn_module(
    "seq2vec_rnn",
    initialize = function(input_size, hidden_size, output_size, num_layers = 1) {
      self$rnn = nn_rnn(input_size, hidden_size, num_layers)
      self$fc = nn_linear(hidden_size, output_size)
    },
    forward = function(x) {
      out = self$rnn(x)
      self$fc(out[ , -1, ]) # only predict for the last RNN values
    }
)

# will need a TorchIngressToken
ingress_token = TorchIngressToken(
    features = features,
    batchgetter = function(x) {
      # convert a data frame to a matrix for the RNN
      
    },
    shape = c(NA, 2)
)


learner = lrn(
    "regr.torch_model",
    network = seq2vec_rnn,
    ingress_tokens = ingress_tokens,
    batch_size = 16,
    epochs = 1,
    device = "cpu"
)
