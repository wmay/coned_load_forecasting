# run some ML distributional regression models for forecasting TV

# use ensemble MOS with neural networks to predict ConEd's TV

# Mostly following Rasp 2018 (10.1175/MWR-D-18-0187.1) and Lerch 2022
# (https://arxiv.org/abs/2204.05102). May also try some sort of RNN for the load
# forecasting part

# A helpful trick: can subtract the lagged portion of TV from the output so it
# doesn't need to be predicted (and add it back afterward)

# install.packages(c('ranger', 'drf', 'RandomForestsGLS'))
# remotes::install_github("mlr-org/mlr3temporal")
# distr6_repos =  c(CRAN = 'https://cloud.r-project.org',
#                   raphaels1 = 'https://raphaels1.r-universe.dev')
# install.packages('distr6', repos = distr6_repos)
setwd('..')
library(magrittr)
library(stars) # also requires ncmeta
# install.packages('ncmeta')
# library(crch)
# library(scoringRules)
# library(beepr)
library(mlr3)
library(mlr3learners)
library(mlr3tuning)
library(mlr3tuningspaces)
library(mlr3temporal)
source('R/mlr3_distr.R')
source('R/mlr3_additions.R') # block CV, mean pinball loss, etc.
source('R/forecast_dataset.R')

all_eff_tmps = read.csv('results/process_station_data/eff_tmp.csv')
all_tvs = read.csv('results/process_station_data/tv.csv')

# see notes in `?benchmark` and `?mlr_tuners_random_search`
lgr::get_logger("mlr3")$set_threshold("warn")
lgr::get_logger("bbotk")$set_threshold("warn")

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

get_valid_day = function(nc, days_ahead) {
  # time is generally read by R as a date, due to the 1 day differences
  attr(nc, 'dimensions')$time$values + days_ahead
}

# note that the day returned here is the valid day
extract_values = function(nc, network, days_ahead) {
  network_centroid = networks %>%
    subset(id == network) %>%
    st_centroid %>%
    st_transform(wgs84) %>%
    merge(data.frame(time = attr(nc, 'dimensions')$time$values))
  nc %>%
    dplyr::slice('edt9pm_day', days_ahead + 1) %>%
    st_extract(network_centroid, bilinear = TRUE, time_column = 'time') %>%
    transform(day = time + days_ahead) %>%
    as.data.frame %>%
    subset(select = -c(time, Shape))
}

# add observed effective temps to short-term TV forecasts
fill_incomplete_tv = function(tv, day, network, days_ahead) {
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
  out_0p25 = extract_values(gefs_0p25_3day, network, days_ahead) %>%
    subset(select = -eff_temp) # 3-day eff_temp is same as TV
  # sometimes soilw is all NA? for now just skip it if so
  if (all(is.na(out_0p25$soilw))) {
    warning('All soilw values are missing')
    out_0p25$soilw = NULL
  }
  out_0p5 = extract_values(gefs_0p5_3day, network, days_ahead)
  out_tv = extract_values(gefs_tv, network, days_ahead)
  # For days_ahead < 2, need to add observed (past) portion of TV to the
  # forecast portion. TV_sd is correct though, since observations have 0 sd
  if (days_ahead < 2) {
    out_tv$TV = with(out_tv, fill_incomplete_tv(TV, day, network, days_ahead))
  }
  out = merge(out_0p25, out_0p5, by = 'day') %>%
    merge(out_tv, by = 'day') %>%
    transform(doy = as.POSIXlt(day)$yday)
  id = paste('TV', network, days_ahead, sep = '_')
  system_tv = all_tvs[, c('day', paste0('network.', network))]
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
  task
}

# strangely `read_mdim` works fine while `read_ncdf` messes up the coordinates
# gefs_tv = read_ncdf('results/process_nwp_data/gefs_tv2.nc')
gefs_tv_members = read_mdim('results/process_nwp_data/gefs_tv_members.nc')
# (this is only needed for evaluating the GEFS)

gefs_0p25_3day = read_nwp_nc('results/process_nwp_data/gefs_0p25_3day_wmean.nc')
gefs_0p5_3day = read_nwp_nc('results/process_nwp_data/gefs_0p5_3day_wmean.nc')
gefs_tv = read_nwp_nc('results/process_nwp_data/gefs_tv.nc')

# get values from stars object at the centroid of a network
wgs84 = 4326
networks = readRDS('results/maps/coned_networks_cleaned.rds')


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
    resampling = rsmp('block_cv'),
    measure = msr('regr.crps'),
    term_evals = 10
)
# system.time(drf_at$train(forecast_tv_tasks[[6]]))

rfgls = LearnerRegrRfgls$new()
# rfgls$param_set$set_values(param_estimate = TRUE, ntree = 100, mtry = 3, nthsize = 16, h = 1, lags = 1)
# I have to restrict the ranges of some of these compared to the other models,
# because this takes so long to fit
rfgls$param_set$set_values(param_estimate = TRUE, mtry.ratio = to_tune(0, 1),
                           ntree = to_tune(1, 500), nthsize = to_tune(10, 20),
                           h = 1, lags = to_tune(1, 2))
rfgls_at = auto_tuner(
    tuner = tnr("random_search"),
    learner = rfgls,
    resampling = rsmp('block_cv', folds = 4),
    measure = msr('regr.crps'),
    term_evals = 5
)


# It'd take too long to run this on all networks and lead times. So I'm randomly
# selecting a few, running the benchmarks, then choosing whichever model
# performs best on average

set.seed(123)
rand_networks = sample(networks$id, 5, replace = TRUE)
rand_days = sample(0:7, 10, replace = TRUE)

# # how long might this take?
# t_all=prepare_tv_task_combined(rand_networks[1])
# system.time(rangerdist$train(t_all))
#  #   user  system elapsed 
#  # 10.119   0.187   9.961 
# system.time(drf$train(t_all))
#  #   user  system elapsed 
#  # 43.316   0.288  13.013 


# RFGLS specifically takes more time, so worth running in parallel

# run the inner loop in parallel and the outer loop sequentially
future::plan(list('sequential', 'multicore'), workers = 4)
# future::plan("sequential")

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
benchmark_tasks =
  lapply(rand_networks, prepare_tv_task_combined, predict_error = TRUE)
bgrid = benchmark_grid(
    tasks = benchmark_tasks,
    # learners = c(ngr, rangerdist_at, drf_at, rfgls_at),
    learners = c(ngr, rangerdist_at, drf_at),
    resamplings = rsmp('forecast_holdout')
)
system.time(bres3 <- progressr::with_progress(benchmark(bgrid)))
bres3$aggregate(c(msr('regr.crps'), msr('regr.mae'), msr('regr.rmse')))
extract_inner_tuning_results(bres3)

# does one stand out as best?
bres3$aggregate(msr('regr.crps')) %>%
  aggregate(regr.crps ~ learner_id, FUN = mean, data = .)
# ranger really beat DRF???

bres3$aggregate(msr('regr.crps')) %>%
  subset(learner_id != 'regr.ngr') %>%
  reshape(direction = 'wide', idvar = 'task_id', timevar = 'learner_id',
          v.names = 'regr.crps')
# total blowout, ranger is dominating

bres3$obs_loss(measures = msr('regr.crps')) %>%
  subset(select = -distr) %>%
  head

bres3$obs_loss() %>%
  subset(select = -distr) %>%
  getElement('se') %>%
  hist()

bres3$obs_loss() %>%
  subset(select = -distr) %>%
  getElement('se') %>%
  range()

msr('regr.crps')

bres3$score(measures = msr('regr.crps')) %>%
  nrow


# rfgls = LearnerRegrRfgls$new()
# rfgls$param_set$set_values(param_estimate = TRUE, ntree = 314,
#                            mtry.ratio = 0.3185942, nthsize = 8, h = 1, lags = 1)

# system.time(rfgls$train(forecast_tv_tasks[[6]]))

# blas library makes a huge difference here--
# libblas: 2500 seconds
# openblas: 840 seconds
# mkl: 8500 seconds omfg

# prediction = rfgls$predict(forecast_tv_tasks[[6]])

# # rfgls$train(forecast_tv_tasks[[1]])

# rfgls$param_set$set_values(param_estimate = TRUE, mtry.ratio = to_tune(0, 1),
#                            ntree = to_tune(1, 500), nthsize = to_tune(8, 20),
#                            h = 1, lags = to_tune(1, 2))


# future::plan('sequential') # for debugging
# bgrid = benchmark_grid(
#     tasks = forecast_tv_tasks,
#     learners = rfgls,
#     resamplings = rsmp('forecast_holdout')
# )
# system.time(bres0 <- benchmark(bgrid, store_models = TRUE))
# # bres0 = benchmark(bgrid, store_models = TRUE)
# bres0$aggregate(c(msr('regr.crps'), msr('regr.mae'), msr('regr.rmse')))

# get hyperparameter tuning results
extract_inner_tuning_results(bres)

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
