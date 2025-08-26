# run some ML distributional regression models for forecasting TV

# use ensemble MOS with neural networks to predict ConEd's TV

# Mostly following Rasp 2018 (10.1175/MWR-D-18-0187.1) and Lerch 2022
# (https://arxiv.org/abs/2204.05102). May also try some sort of RNN for the load
# forecasting part

# A helpful trick: can subtract the lagged portion of TV from the output so it
# doesn't need to be predicted (and add it back afterward)

# pkg_deps = c('ncmeta', 'crch', 'ranger', 'drf', 'RandomForestsGLS',
#              'scoringRules')
# install.packages(pkg_deps)
# remotes::install_github("mlr-org/mlr3temporal")
# distr6_repos =  c(CRAN = 'https://cloud.r-project.org',
#                   raphaels1 = 'https://raphaels1.r-universe.dev')
# install.packages('distr6', repos = distr6_repos)
setwd('..')
library(magrittr)
library(stars) # also requires ncmeta
# library(beepr)
library(future)
library(future.batchtools)
library(mlr3)
library(mlr3learners)
library(mlr3tuning)
library(mlr3tuningspaces)
# library(mlr3temporal)
source('R/mlr3_distr.R')
source('R/mlr3_additions.R') # block CV, mean pinball loss, etc.

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
    transform(day = as.Date(time) + days_ahead) %>%
    subset(select = -c(time, network))
}

# # This is kind of silly, but to get the system (all NYC) predictors, I'm going
# # to find the ConEd networks that contain the weather stations, and just average
# # the networks. Even though it's not the best representation of all NYC, we
# # really are trying to predict the station values.
# extract_system_values = function(version = c('orig', 'new'), days_ahead) {
#   if (version == 'new') {
#     sys_stations = system_results$stids
#   } else {
#     sys_stations = c('NYC', 'LGA')
#   }
#   # which networks are the stations in?
#   network_ids = stations %>%
#     subset(stid %in% sys_stations) %>%
#     st_join(networks) %>%
#     with(setNames(id, stid))
#   # Central Park is not in a ConEd network! Nor is Teterboro
#   missing_dict = c(NYC = '23M', TEB = '1X')
#   missing_idx = which(names(network_ids) %in% names(missing_dict))
#   network_ids[missing_idx] = missing_dict[names(network_ids)[missing_idx]]
#   # get the values at each station, and average them
#   network_predictors =
#     lapply(network_ids, extract_values, days_ahead = days_ahead)
#   parameters = head(names(network_predictors[[1]]), -1)
#   sys_predictors = network_predictors %>%
#     # put into long format
#     lapply(reshape, direction = 'long', varying = parameters, v.names = 'value',
#            idvar = 'day', timevar = 'parameter', times = parameters) %>%
#     do.call(rbind, .) %>%
#     # get averages
#     aggregate(value ~ parameter + day, FUN = mean, na.rm = TRUE) %>%
#     # widen again
#     reshape(direction = 'wide', timevar = 'parameter', idvar = 'day')
#   names(sys_predictors) = sub('value.', '', names(sys_predictors), fixed = TRUE)
#   sys_predictors
# }

# add observed effective temps to short-term TV forecasts
fill_incomplete_tv = function(tv, day, network, days_ahead) {
  return(tv)
  system_eff_tmp = all_eff_tmps[, c('day', paste0('network.', network))]
  names(system_eff_tmp)[2] = 'eff_temp'
  if (days_ahead == 1) {
    tv + .1 * system_eff_tmp$eff_temp[match(day - 2, system_eff_tmp$day)]
  } else if (days_ahead == 0) {
    tv + .2 * system_eff_tmp$eff_temp[match(day - 1, system_eff_tmp$day)] +
      .1 * system_eff_tmp$eff_temp[match(day - 2, system_eff_tmp$day)]
  }
}

prepare_tv_task = function(network, days_ahead, predict_error = TRUE) {
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
    out$TV = with(out, fill_incomplete_tv(TV, day, network, days_ahead))
  }
  out$doy = as.POSIXlt(out$day)$yday
  id = paste('TV', network, days_ahead, sep = '_')
  tv_name = if (is_system) network else paste0('network.', network)
  system_tv = all_tvs[, c('day', tv_name)]
  names(system_tv)[2] = 'tv'
  if (predict_error) {
    # Predict the GEFS forecast error instead of directly predicting TV
    out %>%
      transform(fct_err = TV - system_tv$tv[match(day, system_tv$day)]) %>%
      subset(select = -day) %>%
      na.omit %>%
      as_task_regr('fct_err', id = id)
  } else {
    out %>%
      transform(tv_obs = system_tv$tv[match(day, system_tv$day)]) %>%
      subset(select = -day) %>%
      na.omit %>%
      as_task_regr('tv_obs', id = id)
  }
}

# create a task that combines lead times, as in https://doi.org/10.1002/qj.4701
prepare_tv_task_combined = function(network, predict_error = TRUE) {
  task = prepare_tv_task(network, 0, predict_error = predict_error)
  task$cbind(data = data.frame(days_ahead = rep(0, task$nrow)))
  # append the other lead times
  for (i in 1:7) {
    task_i = prepare_tv_task(network, i, predict_error = TRUE)
    task_i$cbind(data = data.frame(days_ahead = rep(i, task_i$nrow)))
    task$rbind(data = task_i$data())
  }
  task$id = paste('TV', network, sep = '_')
  # make sure resampling is grouped by day_ahead
  task$set_col_roles('days_ahead', add_to = 'stratum')
  task
}


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


# Models

# ngr = LearnerRegrNGR$new()
# rangerdist = LearnerRegrRangerDist$new()
# drf = LearnerRegrDrf$new()

# # How long does it take to run these?
# task0 = prepare_tv_task_combined('3B', predict_error = TRUE)

# system.time(res <- ngr$train(task0))
#   #  user  system elapsed 
#   # 8.383   0.013   8.453 
# system.time(res <- rangerdist$train(task0))
#  #   user  system elapsed 
#  # 11.199   0.112  10.878 
# system.time(res <- drf$train(task0)) # using 4 cores
#  #   user  system elapsed 
#  # 48.991   0.228  14.624 
# drf$param_set$set_values(num.threads = 1)
# system.time(res <- drf$train(task0)) # using 1 core
#  #   user  system elapsed 
#  # 45.852   0.234  46.063

# rfgls = LearnerRegrRfgls$new()
# rfgls$param_set$set_values(param_estimate = TRUE, ntree = 314,
#                            mtry.ratio = 0.3185942, nthsize = 8, h = 1, lags = 1)
# system.time(rfgls$train(forecast_tv_tasks[[6]]))

# blas library makes a huge difference here--
# libblas: 2500 seconds
# openblas: 840 seconds
# mkl: 8500 seconds omfg


# so now we want to autotune these
ngr = LearnerRegrNGR$new()

# ngr_at = auto_tuner(
#     tuner = tnr("random_search"),
#     learner = ngr,
#     resampling = rsmp('block_cv'),
#     measure = msr('regr.crps'),
#     term_evals = 10
# )
# system.time(res <- ngr_at$train(task0)) # using 1 core
# seems to run folds x term_evals = 10 x 10 = 100 individual fits

rangerdist = LearnerRegrRangerDist$new()
# add ranger default search space
rangerdist$param_set$values =
  mlr3misc::insert_named(rangerdist$param_set$values, lts('regr.ranger.default')$values)
rangerdist_at = auto_tuner(
    tuner = tnr("random_search"),
    learner = rangerdist, 
    resampling = rsmp('block_cv', folds = 5),
    measure = msr('regr.crps'),
    term_evals = 10
)

drf = LearnerRegrDrf$new()
# drf hypertuning should be similar to ranger and grf. See
# `lts('regr.ranger.default')` and
# https://grf-labs.github.io/grf/REFERENCE.html#parameter-tuning
drf$param_set$set_values(num.trees = to_tune(1, 2000),
                         sample.fraction = to_tune(0.1, 1),
                         mtry.ratio = to_tune(0, 1),
                         honesty = to_tune(),
                         honesty.fraction = to_tune(.5, .8),
                         honesty.prune.leaves = to_tune(),
                         ci.group.size = 1,
                         num.threads = 1)
drf_at = auto_tuner(
    tuner = tnr("random_search"),
    learner = drf, 
    resampling = rsmp('block_cv', folds = 5),
    measure = msr('regr.crps'),
    term_evals = 10
)
# system.time(drf_at$train(forecast_tv_tasks[[6]]))

# rfgls = LearnerRegrRfgls$new()
# # rfgls$param_set$set_values(param_estimate = TRUE, ntree = 100, mtry = 3, nthsize = 16, h = 1, lags = 1)
# # I have to restrict the ranges of some of these compared to the other models,
# # because this takes so long to fit
# rfgls$param_set$set_values(param_estimate = TRUE, mtry.ratio = to_tune(0, 1),
#                            ntree = to_tune(1, 500), nthsize = to_tune(10, 20),
#                            h = 1, lags = to_tune(1, 2))
# rfgls_at = auto_tuner(
#     tuner = tnr("random_search"),
#     learner = rfgls,
#     resampling = rsmp('block_cv', folds = 4),
#     measure = msr('regr.crps'),
#     term_evals = 5
# )


# # did we split it correctly?
# r1 = rsmp('block_cv')
# r1$instantiate(task0)
# r1$instance %>%
#   subset(task0$data(cols='days_ahead')$days_ahead == 6) %>%
#   tail
# r1$instance %>%
#   subset(task0$data(cols='days_ahead')$days_ahead == 7) %>%
#   head
# # it works!


# It'd take too long to run this on all networks and lead times. So I'm randomly
# selecting a few, running the benchmarks, then choosing whichever model
# performs best on average

set.seed(123)
rand_networks = sample(networks$id, 5, replace = TRUE)
# [1] "4M"  "5Q"  "3X"  "39M" "12M"
rand_days = sample(0:7, 10, replace = TRUE)

# RFGLS specifically takes more time, so worth running in parallel

# run the inner loop in parallel and the outer loop sequentially
# future::plan(list('sequential', 'multicore'), workers = 4)
# future::plan("sequential")

# run outer loop via slurm, inner loop multicore on a slurm node
on_slurm = system('whoami', intern = T) == 'wm177874'
if (on_slurm) {
  # for running on slurm, requesting 10 CPUs per job
  slurm_resources = list(walltime = 2400, memory = 10e3, ncpus = 10)
  plan(list(tweak(batchtools_slurm, resources = slurm_resources, workers = 50),
            multicore))
} else {
  # for testing locally
  plan(list(batchtools_local, multicore), workers = 4)
}

# benchmark_tasks = rand_networks %>%
#   lapply(function(net) {
#     lapply(0:7, function(i) {
#       print(net)
#       print(i)
#       prepare_tv_task(net, i, predict_error = TRUE)
#     })
#   }) %>%
#   do.call(c, .) # put them all in the same list

# tl;dr predicting the forecast error works better than predicting TV directly

# # benchmark methods, first predicting the result directly, then predicting the
# # GEFS forecast error  
# forecast_tv_tasks = lapply(0:7, prepare_tv_task, predict_error = FALSE)
# bgrid = benchmark_grid(
#     tasks = forecast_tv_tasks,
#     learners = c(ngr, rangerdist_at, drf_at),
#     resamplings = rsmp('forecast_holdout')
# )
# system.time(bres <- progressr::with_progress(benchmark(bgrid)))
# bres$aggregate(c(msr('regr.crps'), msr('regr.mae'), msr('regr.rmse')))
# saveRDS(bres, 'benchmarks.rds')

# benchmark_tasks = rand_networks %>%
#   lapply(function(net) {
#     lapply(0:7, function(i) {
#       print(net)
#       print(i)
#       prepare_tv_task(net, i, predict_error = TRUE)
#     })
#   }) %>%
#   do.call(c, .) # put them all in the same list
# # forecast_tv_tasks2 = lapply(0:7, prepare_tv_task, predict_error = TRUE)
# bgrid = benchmark_grid(
#     tasks = benchmark_tasks,
#     learners = c(ngr, rangerdist_at, drf_at, rfgls_at),
#     resamplings = rsmp('forecast_holdout')
# )
# system.time(bres2 <- progressr::with_progress(benchmark(bgrid)))
# bres2$aggregate(c(msr('regr.crps'), msr('regr.mae'), msr('regr.rmse')))
# # this one is consistently a bit better

# bres$aggregate(c(msr('regr.crps'), msr('regr.mae')))
# bres2$aggregate(c(msr('regr.crps'), msr('regr.mae')))


# First step-- model selection

# benchmark_tasks = mapply(prepare_tv_task, rand_networks, rand_days,
#                          MoreArgs = list(predict_error = TRUE))
benchmark_tasks = paste0('network.', rand_networks) %>%
  lapply(prepare_tv_task_combined, predict_error = TRUE)
bgrid = benchmark_grid(
    tasks = benchmark_tasks,
    # learners = c(ngr, rangerdist_at, drf_at, rfgls_at),
    # learners = c(rangerdist_at),
    learners = c(rangerdist_at, drf_at),
    # resamplings = rsmp('forecast_holdout')
    resamplings = rsmp('block_cv', folds = 5)
)
# this requires n_tasks x n_model x n_fold_outer x n_fold_inner x term_evals model fits!
# = 5 * 2 * 5 * 5 * 10 = 2500, yikes!
# n_submissions = n_tasks x n_model x n_fold_outer = 5 * 2 * 5 = 50
# ok that's reasonable
system.time(bres3 <- benchmark(bgrid, store_models = TRUE))
# system.time(bres3 <- progressr::with_progress(benchmark(bgrid, store_models = TRUE)))
saveRDS(bres3, 'results/forecast_tv/benchmarks.rds')
# bres3 = readRDS('results/forecast_tv/benchmarks.rds')
# this is an outrageously large file-- at least 25GB in RAM

bres3$aggregate(c(msr('regr.crps'), msr('regr.mae'), msr('regr.rmse')))

# get the loss by lead time
bres3$obs_loss() %>%
  transform(days_ahead = benchmark_tasks[[1]]$data(rows = row_ids)$days_ahead,
            abs_err = abs(response - truth)) %>%
  aggregate(abs_err ~ days_ahead, FUN = mean, data = .)

# does one stand out as best?
bres3$aggregate(msr('regr.crps')) %>%
  aggregate(regr.crps ~ learner_id, FUN = mean, data = .)
# very very close

# bres3$aggregate(msr('regr.crps')) %>%
#   subset(learner_id != 'regr.ngr') %>%
#   subset(select = c(task_id, learner_id, regr.crps)) %>%
#   reshape(direction = 'wide', idvar = 'task_id', timevar = 'learner_id',
#           v.names = 'regr.crps')

# get the hyperparameters
drf_res_params = extract_inner_tuning_results(bres3) %>%
  subset(learner_id == 'regr.drf.tuned')

# two logical values, honesty and honesty.prune.leaves
table(drf_res_params$honesty)
table(drf_res_params[, c('honesty', 'honesty.prune.leaves')])
# alright, we're going with honesty and honesty.prune.leaves

drf_params = drf_res_params %>%
  subset(honesty & honesty.prune.leaves) %>%
  subset(select = c('mtry.ratio', 'num.trees', 'sample.fraction',
                    'honesty.fraction')) %>%
  colMeans %>%
  as.list
drf_params$num.trees = round(drf_params$num.trees)
drf_params[c('honesty', 'honesty.prune.leaves')] = TRUE
saveRDS(drf_params, 'results/forecast_tv/tv_hyperparameters.rds')
drf_params = readRDS('results/forecast_tv/tv_hyperparameters.rds')

# now we create the model with the hyperparameters, and fit it to all 70
# networks + 2 system TVs
drf_params = c(drf_params, list(ci.group.size = 1, num.threads = 4))
make_drf = function(network) {
  # message('starting ', network)
  drf = LearnerRegrDrf$new()
  do.call(drf$param_set$set_values, drf_params)
  task = prepare_tv_task_combined(network, predict_error = TRUE)
  drf$train(task)
  drf
}

library(pbapply)
tv_models = names(all_tvs)[-1] %>%
  sub('network.', '', ., fixed = TRUE) %>%
  pblapply(make_drf)
names(tv_models) = names(all_tvs)[-1] %>%
  sub('network.', '', ., fixed = TRUE)
saveRDS(tv_models, 'results/forecast_tv/tv_models.rds')


# slightly different for GEFS benchmarking

# get a matrix of GEFS TV forecasts
get_gefs_samples = function(nc, lon = 3, lat = 4, days_ahead = 2) {
  # When read by `read_mdim` time values are read as dates despite having a time
  # as well. However the numeric date will have a fractional component!
  dates = attr(nc, 'dimensions')$time$values + days_ahead
  # days ahead starts at 0 in the array
  # if (days_ahead < 2) stop('TV only available starting day 2')
  tv = t(nc[['TV']][lon, lat,, days_ahead + 1, ])
  # For days_ahead < 2, need to add observed (past) portion of TV to the
  # forecast portion. TV_sd is correct though, since observations have 0 sd
  if (days_ahead == 1) {
    tv = tv +
      .1 * system_tv$tv[match(dates - 2, system_tv$day)]
  }
  if (days_ahead == 0) {
    tv = tv +
      .2 * system_tv$tv[match(dates - 1, system_tv$day)] +
      .1 * system_tv$tv[match(dates - 2, system_tv$day)]
  }
  list(day = dates, tv = tv)
}

make_gefs_tasks = function() {
  sapply(0:7, function(ahead) {
    id = paste0('Forecast TV ', ahead)
    gefs = get_gefs_samples(gefs_tv3, days_ahead = ahead)
    # make sure days are exactly the same as the prediction tasks
    # compare_days = prepare_multiday_dataset(ahead)$day
    compare_days = get_valid_day(gefs_tv, ahead)
    x = gefs$tv
    missing_row = apply(x, 1, function(x) all(is.na(x)))
    y = system_tv[match(gefs$day, system_tv$day), 'tv', drop = FALSE]
    cbind(as.data.frame(x), y) %>%
      subset(!missing_row) %>%
      subset(!is.na(tv)) %>%
      as_task_regr('tv', id = id)
  })
}
gefs_tasks = make_gefs_tasks()

gefs_tester = LearnerRegrIdent$new()
gefs_tester$predict_type = 'distr'
gefs_tester$id = 'GEFS'
bgrid = benchmark_grid(
    tasks = gefs_tasks,
    learners = gefs_tester,
    resamplings = rsmp('forecast_holdout')
)
bres_gefs = benchmark(bgrid)
bres_gefs$aggregate(c(msr('regr.crps'), msr('regr.mae'), msr('regr.rmse')))


all_res = rbind(
    bres2$aggregate(c(msr('regr.crps'), msr('regr.mae'), msr('regr.rmse'))),
    bres_gefs$aggregate(c(msr('regr.crps'), msr('regr.mae'), msr('regr.rmse')))
)
saveRDS(all_res, 'results/forecast_tv/benchmarks.rds')

all_res %>%
  subset(select = c(task_id, learner_id, regr.crps)) %>%
  reshape(direction = 'wide', idvar = 'task_id', timevar = 'learner_id')


# rfgls = LearnerRegrRangerDist$new()

# bgrid = benchmark_grid(
#     tasks = forecast_tv_tasks[5:6],
#     learners = rfgls,
#     resamplings = rsmp('forecast_holdout')
# )
# bres2 = benchmark(bgrid)

# rfgls$train(forecast_tv_tasks[[1]])


# # train the model
# qrf$train(forecast_tv, split$train_set)
# prediction = qrf$predict(forecast_tv, split$test_set)

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




# library(ggplot2)
# library(tidyr)

# comparison_plot = function(metrics, nwp_acc) {
#   nwp_acc_long = nwp_acc %>%
#     transform(days_ahead = day, model = 'GEFS') %>%
#     subset(select = -day) %>%
#     gather('metric', 'value', -c(model, days_ahead))
#   names(metrics) = sub('regr\\.', '', names(metrics))
#   names(metrics)[2:3] = toupper(names(metrics)[2:3])
#   metrics_long = metrics %>%
#     transform(model = 'Random forest') %>%
#     gather('metric', 'value', -c(model, days_ahead))

#   rbind(nwp_acc_long, metrics_long) %>%
#     ggplot(aes(x = days_ahead, y = value, color = model)) + geom_line() + facet_wrap( ~ metric)
# }

# comparison_plot(rf_all_results$metrics, nwp_acc)


# # let's try out torch and sequence-to-vector RNN
# library(mlr3torch)
# library(torch)

# # create an RNN model using the regr.torch_model class

# seq2vec_rnn = nn_module(
#     "seq2vec_rnn",
#     initialize = function(input_size, hidden_size, output_size, num_layers = 1) {
#       self$rnn = nn_rnn(input_size, hidden_size, num_layers)
#       self$fc = nn_linear(hidden_size, output_size)
#     },
#     forward = function(x) {
#       out = self$rnn(x)
#       self$fc(out[ , -1, ]) # only predict for the last RNN values
#     }
# )

# # will need a TorchIngressToken
# ingress_token = TorchIngressToken(
#     features = features,
#     batchgetter = function(x) {
#       # convert a data frame to a matrix for the RNN
      
#     },
#     shape = c(NA, 2)
# )

# learner = lrn(
#     "regr.torch_model",
#     network = seq2vec_rnn,
#     ingress_tokens = ingress_tokens,
#     batch_size = 16,
#     epochs = 1,
#     device = "cpu"
# )
