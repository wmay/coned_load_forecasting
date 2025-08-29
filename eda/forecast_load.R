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

# pkg_deps = c('ncmeta', 'crch', 'ranger', 'drf', 'RandomForestsGLS',
#              'scoringRules')
# install.packages(pkg_deps)
# remotes::install_github("mlr-org/mlr3temporal")
# distr6_repos =  c(CRAN = 'https://cloud.r-project.org',
#                   raphaels1 = 'https://raphaels1.r-universe.dev')
# install.packages('distr6', repos = distr6_repos)
setwd('..')
library(timeDate) # holidays
library(magrittr)
library(stars) # also requires ncmeta
library(future)
# library(future.batchtools)
library(mlr3)
library(torch)
library(mlr3torch)
library(mlr3learners)
library(mlr3tuning)
library(mlr3tuningspaces)
library(mlr3pipelines)
source('R/mlr3_additions.R') # block CV, mean pinball loss, etc.
source('R/mlr3_distr.R')

# see notes in `?benchmark` and `?mlr_tuners_random_search`
lgr::get_logger("mlr3")$set_threshold("warn")
lgr::get_logger("bbotk")$set_threshold("warn")

all_eff_tmps = read.csv('results/process_station_data/eff_tmp.csv')
all_tvs = read.csv('results/process_station_data/tv.csv')

# if coordinates are not in the same order for every variable, `read_mdim` works
# fine while `read_ncdf` messes up the coordinates
gefs_tv_members = read_ncdf('results/process_nwp_data/gefs_tv_members.nc')
# ^ this is only needed for evaluating the GEFS
# gefs_daily = read_ncdf('results/process_nwp_data/gefs_daily.nc')
gefs_3day = read_ncdf('results/process_nwp_data/gefs_3day_wmean.nc')

# get values from stars object at the centroid of a network
networks = readRDS('results/maps/coned_networks_cleaned.rds')
stations = readRDS('results/station_data/stations.rds')
system_results = readRDS('results/select_stations/system_results.rds')
gam_fct = readRDS('results/load_curve_forecast/gam_forecasts.rds')

get_valid_day = function(nc, days_ahead) {
  # time is generally read by R as a date, due to the 1 day differences
  attr(nc, 'dimensions')$time$values + days_ahead
}

# note that the day returned here is the valid day
extract_values = function(network, days_ahead) {
  network_idx = which(attr(gefs_3day, 'dimensions')$network$values == network)
  day_idx = days_ahead + 1
  gefs_3day %>%
    dplyr::slice('edt9pm_day', day_idx) %>%
    dplyr::slice('network', network_idx) %>%
    as.data.frame %>%
    # getElement('time') %>%
    # class
    transform(valid_day = as.Date(time) + days_ahead) %>%
    subset(select = -c(time, network))
}

# add observed effective temps to short-term TV forecasts
fill_incomplete_tv = function(tv, day, network, days_ahead) {
  # return(tv)
  system_eff_tmp = all_eff_tmps[, c('day', paste0('network.', network))]
  names(system_eff_tmp)[2] = 'eff_temp'
  if (days_ahead == 1) {
    tv + .1 * system_eff_tmp$eff_temp[match(day - 2, system_eff_tmp$day)]
  } else if (days_ahead == 0) {
    tv + .2 * system_eff_tmp$eff_temp[match(day - 1, system_eff_tmp$day)] +
      .1 * system_eff_tmp$eff_temp[match(day - 2, system_eff_tmp$day)]
  }
}

# prepare_dataset = function(days_ahead = 1) {
#   x = make_predictor_dataset2(gefs, 12, 3, 4, days_ahead = days_ahead)
#   # add doy to account for seasonality
#   x$doy = as.POSIXlt(x$day)$yday
#   # y is the forecast error
#   y = x %>%
#     transform(obs_eff_temp = system_eff_tmp$eff_temp[match(day, system_eff_tmp$day)]) %>%
#     transform(eff_tmp_fc_err = eff_temp - obs_eff_temp) %>%
#     subset(select = eff_tmp_fc_err)
#   ds = na.omit(cbind(x, y))
#   x = subset(ds, select = -c(day, eff_tmp_fc_err))
#   y = subset(ds, select = eff_tmp_fc_err)
#   list(x = x, y = y, day = ds$day)
# }

# average_columns = function(dat, x) {
#   .1 * dat[, paste0(x, '.1')] + .2 * dat[, paste0(x, '.2')] +
#     .7 * dat[, paste0(x, '.3')]
# }

# average_3days = function(dat) {
#   columns = names(dat)[endsWith(names(dat), '.1')]
#   for (n in columns) {
#     new_name = substr(n, 1, nchar(n) - 2)
#     old_names = paste(new_name, 1:3, sep = '.')
#     dat[, new_name] = average_columns(dat, new_name)
#     dat[, old_names] = NULL
#   }
#   dat
# }

# # Just like for TV forecasting, but use GAM prediction error as the outcome and
# # add GAM outputs as predictors. Also only use business days
# prepare_multiday_dataset = function(days_ahead = 2, predict_error = TRUE) {
#   # x: make a predictor dataset for each day, then combine
#   d_list = lapply((days_ahead - 2):days_ahead, prepare_dataset)
#   # "day" starts off as the valid day for each forecast. Change it to TV day 3,
#   # that is, the last of the 3 days used to calculate TV
#   for (i in 1:2) d_list[[i]]$day = d_list[[i]]$day + (3 - i)

#   # # what was wrong with this?
#   # d_list = lapply(days_ahead:(days_ahead - 2), prepare_dataset)
#   # for (i in 2:3) d_list[[i]]$day = d_list[[i]]$day + (i - 1)
  
#   # for (i in 1:3) d_list[[i]]$tv_day = i
#   for (i in 1:3) d_list[[i]]$tv_day = (4 - i)
#   v.names = with(d_list[[1]], c(names(x), names(y)))
#   gam_loads = gam_fct$forecasts[, days_ahead - 1, ] %>%
#     as.data.frame %>%
#     cbind(day = gam_fct$model_day + days_ahead) %>%
#     subset(isBizday(as.timeDate(day), holidays = holidayNYSE(2021:2024)))
#   # gam_loads = data.frame(day = gam_fct$model_day + days_ahead)
#   tv = make_nwp_tv_dataset(gefs_tv2, 3, 4, days_ahead = days_ahead)
#   out = d_list %>%
#     lapply(function(x) with(x, cbind(x, y, tv_day, day))) %>%
#     do.call(rbind, .) %>%
#     reshape(direction = "wide", idvar = 'day', timevar = "tv_day",
#             v.names = v.names) %>%
#     # subset(select = -c(eff_temp.1, eff_temp.2, eff_temp.3)) %>%
#     subset(select = -c(doy.1, doy.2, doy.3)) %>%
#     transform(TV = tv$tv[match(day, tv$day)],
#               TV_sd = tv$tv_sd[match(day, tv$day)],
#               load = gam_loads$fit[match(day, gam_loads$day)],
#               load_se = gam_loads$se.fit[match(day, gam_loads$day)],
#               load_err = gam_loads$err[match(day, gam_loads$day)],
#               load_obs = gam_loads$obs[match(day, gam_loads$day)]) %>%
#     # transform(tv_fc_err = TV - system_tv$tv[match(day, system_tv$day)]) %>%
#     subset(select = -c(eff_tmp_fc_err.1, eff_tmp_fc_err.2, eff_tmp_fc_err.3)) %>%
#     na.omit %>%
#     average_3days %>%
#     subset(select = -eff_temp)

#   if (predict_error) {
#     x = subset(out, select = -c(day, load_err, load_obs))
#     y = subset(out, select = load_err)
#   } else {
#     # no GAM prediction in this case
#     x = subset(out, select = -c(day, load, load_se, load_err, load_obs))
#     y = subset(out, select = load_obs)
#   }
#   x$doy = as.POSIXlt(out$day)$yday
#   list(x = x, y = y, day = out$day)
#   # # y: get either TV, or the future portion of TV
#   # tv = system_tv[match(out$day, system_tv$day), 'tv']
#   # list(x = out, y = tv)
# }

prepare_load_task = function(network, days_ahead, use_gam = FALSE) {
  # Note: netcdf files use forecast day, while `system_tv` uses valid day. Must
  # be consistent!
  is_system = startsWith(network, 'system')
  extract_name = if (is_system) 'system' else network
  out = extract_values(extract_name, days_ahead) %>%
    subset(select = -eff_temp) # 3-day eff_temp is same as TV
  # sometimes soilw is all NA? for now just skip it if so
  if (all(is.na(out$soilw))) {
    # warning('All soilw values are missing')
    out$soilw = NULL
  }
  # For days_ahead < 2, need to add observed (past) portion of TV to the
  # forecast portion. TV_sd is correct though, since observations have 0 sd
  if (days_ahead < 2) {
    out$TV = with(out, fill_incomplete_tv(TV, valid_day, network, days_ahead))
  }
  out$doy = as.POSIXlt(out$valid_day)$yday
  id = paste('TV', network, days_ahead, sep = '_')
  tv_name = if (is_system) network else paste0('network.', network)
  system_tv = all_tvs[, c('day', tv_name)]
  names(system_tv)[2] = 'tv'
  # get load data
  gam_loads = gam_fct$forecasts[, days_ahead + 1, ] %>%
    as.data.frame %>%
    cbind(day = gam_fct$model_day + days_ahead) %>%
    subset(isBizday(as.timeDate(day), holidays = holidayNYSE(2021:2024)))
  names(gam_loads)[1:4] = paste0('load_', c('pred', 'se', 'err', 'obs'))
  if (use_gam) {
    # Predict the error of the GAM instead of directly predicting load. But we
    # need to leave the observation here so that mlr3 can calculate the MAPE.
    # The error term is calculated later in an mlr3 pipeop.
    # out$load_obs = NULL
    # as_task_regr(out, 'load_err', id = id)
    out %>%
      transform(fct_err = TV - system_tv$tv[match(day, system_tv$day)]) %>%
      subset(select = -day) %>%
      na.omit %>%
      as_task_regr('fct_err', id = id)
  } else {
    out %>%
      transform(load = gam_loads$load_obs[match(valid_day, gam_loads$day)]) %>%
      subset(select = -valid_day) %>%
      na.omit %>%
      as_task_regr('load', id = id)
  }
}


# prepare_load_task = function(days_ahead = 2, use_gam = TRUE) {
#   gam_loads = gam_fct$forecasts[, days_ahead + 1, ] %>%
#     as.data.frame %>%
#     cbind(day = gam_fct$model_day + days_ahead) %>%
#     subset(isBizday(as.timeDate(day), holidays = holidayNYSE(2021:2024)))
#   names(gam_loads)[1:4] = paste0('load_', c('pred', 'se', 'err', 'obs'))
#   out_0p25 = gefs_0p25_3day %>%
#     lapply(function(x) x[3, 4,, days_ahead + 1]) %>%
#     as.data.frame %>%
#     subset(select = -eff_temp) %>% # 3-day eff_temp is same as TV
#     transform(day = get_valid_day(gefs_0p25_3day, days_ahead))
#   out_0p5 = gefs_0p5_3day %>%
#     lapply(function(x) x[1, 1,, days_ahead + 1]) %>%
#     as.data.frame %>%
#     transform(day = get_valid_day(gefs_0p5_3day, days_ahead))
#   out_tv = gefs_tv %>%
#     lapply(function(x) x[3, 4, days_ahead + 1, ]) %>%
#     as.data.frame %>%
#     transform(day = get_valid_day(gefs_tv, days_ahead))
#   # For days_ahead < 2, need to add observed (past) portion of TV to the
#   # forecast portion. TV_sd is correct though, since observations have 0 sd
#   if (days_ahead == 1) {
#     out_tv$TV = out_tv$TV +
#       .1 * system_tv$tv[match(out_tv$day, system_tv$day + 1)]
#   }
#   if (days_ahead == 0) {
#     out_tv$TV = out_tv$TV +
#       .2 * system_tv$tv[match(out_tv$day, system_tv$day + 1)] +
#       .1 * system_tv$tv[match(out_tv$day, system_tv$day + 2)]
#   }
#   out = merge(out_0p25, out_0p5, by = 'day') %>%
#     merge(out_tv, by = 'day') %>%
#     merge(gam_loads, by = 'day') %>%
#     transform(doy = as.POSIXlt(day)$yday) %>%
#     subset(select = -day) %>%
#     na.omit
#   id = paste0('Forecast TV ', days_ahead)
#   if (use_gam) {
#     # Predict the error of the GAM instead of directly predicting load. But we
#     # need to leave the observation here so that mlr3 can calculate the MAPE.
#     # The error term is calculated later in an mlr3 pipeop.
#     # out$load_obs = NULL
#     # as_task_regr(out, 'load_err', id = id)
#     out$load_err = NULL
#   } else {
#     out[, c('load_pred', 'load_se', 'load_err')] = NULL
#   }
#   as_task_regr(out, 'load_obs', id = id)
# }




# out = list(model_day = forecast_days, gam_forecasts)


# loads = read.csv('data/coned/Borough and System Data 2020-2024.csv') %>%
#   transform(DT = as.POSIXct(DT, tz = 'EST5EDT', '%m/%d/%Y %H:%M'),
#             # some inconsistency, but various forms of 'False' mean data is fine
#             BAD = !startsWith(BAD, 'F')) %>%
#   subset(!is.na(DT) & !BAD)
# names(loads) = tolower(names(loads))

# peaks = loads %>%
#   subset(borough == 'CECONY') %>%
#   transform(day = as.Date(dt, tz = 'EST5EDT')) %>%
#   aggregate(reading ~ day, ., max) %>%
#   subset(isBizday(as.timeDate(day), holidays = holidayNYSE(2021:2024)))

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



# Ok forget all the GAM stuff for now. Just running a neural network. Including
# time predictor to allow for nonstationarity. Will have to run incremental
# updates over time. For this version, neither NGR nor random forest make much
# sense (nonlinear main effect, and requires extrapolation)

# just doing an MLP now, but will want to change to DRN. Also why does MLP only
# support dropout? Probably need to extend base torch class to support
# distributional regression
lrn("regr.mlp", ...)

# needed for evaluating nonstationary model
rsmp("forecast_cv")





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
# all.tunable.params <- c("sample.fraction", "mtry", "min.node.size", "honesty.fraction",
#                         "honesty.prune.leaves", "alpha", "imbalance.penalty")



# future::plan('multicore', workers = 4)
# future::plan("multisession", workers = 4)

# run the inner loop in parallel and the outer loop sequentially
# future::plan(list("sequential", "multisession"), workers = 4)
future::plan(list("sequential", "multicore"), workers = 4)
# future::plan("sequential")

# We are going to benchmark the models, first without using the GAM, and then
# with

forecast_load_tasks = lapply(0:7, prepare_load_task, use_gam = FALSE)
bgrid = benchmark_grid(
    tasks = forecast_load_tasks,
    #learners = ngr,
    learners = c(ngr, rangerdist_at, drf_at),
    resamplings = rsmp('forecast_holdout')
)
system.time(bres <- progressr::with_progress(benchmark(bgrid)))
bres$aggregate(c(msr('regr.crps'), msr('regr.mape'), msr('regr.mae')))
# saveRDS(bres, 'benchmarks.rds')

# try out the GAM+pipeop
forecast_load_tasks2 = lapply(0:7, prepare_load_task, use_gam = TRUE)
subtract_gam = PipeOpSubtractColumn$new(param_vals = list(column = 'load_pred'))

model_pipelines = lapply(c(ngr, rangerdist_at, drf_at), function(x) {
  pipeline_targettrafo(PipeOpLearner$new(ngr),
                       trafo_pipeop = subtract_gam)
})

bgrid = benchmark_grid(
    tasks = forecast_load_tasks2,
    #learners = model_pipelines[[1]],
    learners = model_pipelines,
    resamplings = rsmp('forecast_holdout')
)
system.time(bres2 <- progressr::with_progress(benchmark(bgrid)))
bres2$aggregate(c(msr('regr.crps'), msr('regr.mape'), msr('regr.mae')))
# So this is much better for NGR (predictably), but worse for everything else?
# How strange.


# # get hyperparameter tuning results
# extract_inner_tuning_results(bres2)

# benchmark the GAM forecasts for comparison
make_gam_tasks = function() {
  sapply(0:7, function(days_ahead) {
    id = paste0('Forecast TV ', days_ahead)
    gam_loads = gam_fct$forecasts[, days_ahead + 1, ] %>%
      as.data.frame %>%
      cbind(day = gam_fct$model_day + days_ahead) %>%
      subset(isBizday(as.timeDate(day), holidays = holidayNYSE(2021:2024)))
    names(gam_loads)[1:4] = paste0('load_', c('pred', 'se', 'err', 'obs'))
    gam_loads %>%
      subset(select = c(load_pred, load_obs)) %>%
      na.omit %>%
      as_task_regr('load_obs', id = id)
  })
}

gam_tasks = make_gam_tasks()

gam_tester = LearnerRegrIdent$new()
gam_tester$predict_type = 'response'
gam_tester$id = 'GEFS+GAM'
bgrid = benchmark_grid(
    tasks = gam_tasks,
    learners = gam_tester,
    resamplings = rsmp('forecast_holdout')
)
bres_gam = benchmark(bgrid)
bres_gam$aggregate(c(msr('regr.mape'), msr('regr.mae')))

load_measures = c(msr('regr.crps'), msr('regr.mae'), msr('regr.mape'),
                  msr('regr.rmse'))
all_res = rbind(
    bres$aggregate(load_measures),
    bres2$aggregate(load_measures),
    bres_gam$aggregate(load_measures)
)
saveRDS(all_res, 'results/forecast_load/benchmarks.rds')

all_res %>%
  subset(select = c(task_id, learner_id, regr.crps)) %>%
  reshape(direction = 'wide', idvar = 'task_id', timevar = 'learner_id')





library(ggplot2)
# library(tidyr)

gefs_baseline = all_res %>%
  subset(learner_id == 'GEFS+GAM') %>%
  transform(days_ahead = as.integer(substr(task_id, 13, 13))) %>%
  with(setNames(regr.mae, days_ahead))

all_res %>%
  transform(days_ahead = as.integer(substr(task_id, 13, 13))) %>%
  ggplot(aes(x = days_ahead, y = regr.crps, color = learner_id)) +
  geom_line()

all_res %>%
  # remove the GAM-based models
  subset(!startsWith(learner_id, 'subtract')) %>%
  subset(!endsWith(learner_id, 'ngr')) %>%
  transform(days_ahead = as.integer(substr(task_id, 13, 13))) %>%
  transform(mae = regr.mae / (gefs_baseline[as.character(days_ahead)])) %>%
  ggplot(aes(x = days_ahead, y = mae, color = learner_id)) +
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
