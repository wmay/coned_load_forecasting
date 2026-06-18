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
# setwd('..')
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
all_tvs = read.csv('results/process_station_data/tv.csv') %>%
  transform(day = as.Date(day))

# if coordinates are not in the same order for every variable, `read_mdim` works
# fine while `read_ncdf` messes up the coordinates
gefs_tv_members = read_ncdf('results/process_nwp_data/gefs_tv_members.nc')
# ^ this is only needed for evaluating the GEFS
# gefs_daily = read_ncdf('results/process_nwp_data/gefs_daily.nc')
gefs_3day = read_ncdf('results/process_nwp_data/gefs_3day_wmean.nc')

# get values from stars object at the centroid of a network
networks = readRDS('results/maps/coned_networks_cleaned.rds')
stations = readRDS('results/station_data/stations.rds')
# system_results = readRDS('results/select_stations/system_results.rds')


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
  is_system = startsWith(network, 'system')
  eff_tmps_name = if (is_system) network else paste0('network.', network)
  system_eff_tmp = all_eff_tmps[, c('day', eff_tmps_name)]
  names(system_eff_tmp)[2] = 'eff_temp'
  if (days_ahead == 1) {
    tv + .1 * system_eff_tmp$eff_temp[match(day - 2, system_eff_tmp$day)]
  } else if (days_ahead == 0) {
    tv + .2 * system_eff_tmp$eff_temp[match(day - 1, system_eff_tmp$day)] +
      .1 * system_eff_tmp$eff_temp[match(day - 2, system_eff_tmp$day)]
  } else {
    tv
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


# Models

# ngr = LearnerRegrNGR$new()
# rangerdist = LearnerRegrRangerDist$new()
# drf = LearnerRegrDrf$new()

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
# DRF sometimes fails during training, so it needs a fallback to make
# autotuning/benchmarking run reliably
drf$encapsulate('try', fallback = lrn("regr.featureless_distr"))
drf_at = auto_tuner(
    tuner = tnr("random_search"),
    learner = drf, 
    resampling = rsmp('block_cv', folds = 5),
    measure = msr('regr.crps'),
    term_evals = 5
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


# It'd take too long to run this on all networks and lead times. So I'm randomly
# selecting a few, running the benchmarks, then choosing whichever model
# performs best on average

set.seed(123)
rand_networks = sample(networks$id, 5, replace = TRUE)
# [1] "4M"  "5Q"  "3X"  "39M" "12M"
rand_days = sample(0:7, 10, replace = TRUE)

# run outer loop via slurm, inner loop multicore on a slurm node
on_slurm = system('whoami', intern = T) == 'wm177874'
if (on_slurm) {
  # for running on slurm, requesting 10 CPUs per job
  slurm_resources = list(walltime = 2400, memory = 10e3, ncpus = 17)
  plan(list(
      tweak(batchtools_slurm, resources = slurm_resources, workers = 50),
      multicore
  ))
  # plan(list(
  #     tweak(multicore, workers = availableCores() %/% 10),
  #     tweak(multicore, workers = I(10))
  # ))
} else {
  # for testing locally
  n_workers = parallel::detectCores() - 1
  # plan(list(batchtools_local, multicore), workers = 4)
  plan(list('sequential', 'multicore'), workers = n_workers)
}

# First step-- model selection

benchmark_tasks = lapply(rand_networks, prepare_tv_task_combined,
                         predict_error = TRUE)
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

# does one stand out as best?
bres3$aggregate(msr('regr.crps')) %>%
  aggregate(regr.crps ~ learner_id, FUN = mean, data = .)
# very very close

bres3_summ = bres3$aggregate() %>%
  transform(network = sub('.*_', '', task_id))
drf_samples = bres3_summ %>%
  subset(learner_id == 'regr.drf.tuned') %>%
  getElement('nr')
# summarize drf benchmarking results for plotting
tv_losses = bres3$obs_loss(msr('regr.mae')) %>%
  subset(resample_result %in% drf_samples, select = -distr) %>%
  # have to calculate crps manually for now
  transform(regr.crps = scoringRules::crps_norm(truth, response, se)) %>%
  # get the loss by lead time
  transform(days_ahead = benchmark_tasks[[1]]$data(rows = row_ids)$days_ahead) %>%
  by(.[, c('days_ahead', 'resample_result')], function(x) {
    network = bres3_summ$network[bres3_summ$nr == x$resample_result[1]]
    data.frame(days_ahead = x$days_ahead[1],
               network = network,
               regr.mae = mean(x$regr.mae),
               regr.crps = mean(x$regr.crps))
  }) %>%
  as.list %>%
  do.call(rbind, .)
write.csv(tv_losses, file = 'results/forecast_tv/tv_losses.csv',
          row.names = F)

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
# alright, we're going with no honesty

drf_params = drf_res_params %>%
  subset(select = c('mtry.ratio', 'num.trees', 'sample.fraction')) %>%
  colMeans %>%
  as.list
drf_params$num.trees = round(drf_params$num.trees)
drf_params$honesty = FALSE
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


# GEFS benchmarking

# as long as we're using cross validation above for the ML benchmarking, we
# don't have to worry about matching the resampling for the NWP benchmarking.
# It's equivalent to testing the accuracy over the full dataset

# Get GEFS forecasts for the given day+network+lead time, as distr6. If no
# forecast is available, returns a forecast of -1 (because NA breaks things)
get_gefs_fcts = function(nc, network, days_ahead, valid_date) {
  nc_network = if (startsWith(network, 'system')) 'system' else network
  net_ind = which(attr(nc, 'dimensions')$network$values == nc_network)
  nc_valid_dates = as.Date(attr(nc, 'dimensions')$time$values) + days_ahead
  day_ind = which(nc_valid_dates == valid_date)
  # days ahead starts at 0 in the array
  tvs = nc[['TV']][net_ind,, days_ahead + 1, day_ind]
  if (days_ahead < 2) {
    # add observed (past) portion of TV to the forecast portion
    tvs = fill_incomplete_tv(tvs, valid_date, network, days_ahead)
  }
  if (all(is.na(tvs))) {
    warning('all GEFS members are missing: network ', network, ', valid_date ',
            valid_date, ', days_ahead ', days_ahead)
    return(distr6::Empirical$new(samples = -1))
  }
  if (any(is.na(tvs))) warning('some GEFS members are missing')
  # why are there some missing values?
  # convert to distr6 for mlr3
  distr6::Empirical$new(samples = na.omit(tvs))
}

# GEFS model, takes day+network+lead time as predictors, returns the GEFS
# forecasts
LearnerRegrGEFS = R6::R6Class(
  "LearnerRegrGEFS",
  inherit = LearnerRegr,
  public = list(
    initialize = function() {
      self$state$model = TRUE # no training required
      super$initialize(
        id = "regr.ident",
        feature_types = c("Date", "character", "numeric"),
        predict_types = 'distr',
        packages = c('magrittr', 'distr6'),
        properties = c("weights", "missings"),
        label = "GEFS"
      )
    }
  ),
  private = list(
      .predict = function(task) {
        distrs = task$data(cols = task$feature_names) %>%
          split(1:nrow(.)) %>%
          lapply(function(x) {
            get_gefs_fcts(gefs_tv_members, x$network, x$days_ahead, x$day)
          }) %>%
          distr6::VectorDistribution$new(distlist = .)
        list(distr = distrs)
      }
  )
)

# make a task where the predictors are day+network+leadtime, and the outcomes
# are observations
make_gefs_network_task = function(network) {
  subtasks = sapply(0:7, function(ahead) {
    id = paste('Forecast TV', network, ahead)
    # gefs = get_gefs_fcts(gefs_tv_members, network, ahead)
    if (startsWith(network, 'system')) {
      tv_col = network
    } else {
      tv_col = paste0('network.', network)
    }
    valid_dates = all_tvs$day
    obs_tv = all_tvs[, tv_col]
    data.frame(day = valid_dates, network = network, days_ahead = ahead,
               tv = obs_tv) %>%
      subset(!is.na(tv)) %>%
      as_task_regr('tv', id = id)
  })
  # combine lead times
  task = subtasks[[1]]
  task$cbind(data = data.frame(days_ahead = rep(0, task$nrow)))
  # append the other lead times
  for (i in 1:7) {
    task_i = subtasks[[i + 1]]
    task_i$cbind(data = data.frame(days_ahead = rep(i, task_i$nrow)))
    task$rbind(data = task_i$data())
  }
  task$id = paste('TV', network, sep = '_')
  task
}


# benchmarking approach: testing the random networks from before, plus the
# system TV

gefs_benchmark_tasks = rand_networks %>%
  c('system.orig', 'system.new') %>%
  lapply(make_gefs_network_task)

gefs_model = LearnerRegrGEFS$new()

gefs_preds = lapply(gefs_benchmark_tasks, function(x) gefs_model$predict(x))

# summarize GEFS benchmarking results for plotting
extract_gefs_losses = function(pred, task) {
  pred$obs_loss(msr('regr.mae')) %>%
    # have to calculate crps manually for now
    transform(regr.crps = mapply(function(x, y) {
      samples = x$getParameterValue('data')$samples
      scoringRules::crps_sample(y, samples)
    }, distr[[1]]$wrappedModels(), truth)) %>%
    subset(select = -distr) %>%
    # remove NA values, marked by -1 forecast value
    subset(response > 0) %>%
    # get the loss by lead time
    transform(days_ahead = task$data(rows = row_ids)$days_ahead,
              network = task$data(rows = row_ids)$network) %>%
    as.data.frame %>%
    by(.[c('days_ahead', 'network')], function(x) {
      data.frame(days_ahead = x$days_ahead[1],
                 network = x$network[1],
                 regr.mae = mean(x$regr.mae),
                 regr.crps = mean(x$regr.crps))
    }) %>%
    as.list %>%
    do.call(rbind, .)
}

gefs_losses = mapply(extract_gefs_losses, gefs_preds, gefs_benchmark_tasks,
                     SIMPLIFY = FALSE) %>%
  do.call(rbind, .)
write.csv(gefs_losses, file = 'results/forecast_tv/gefs_losses.csv',
          row.names = F)


# add DRF benchmarking for comparison with GEFS
benchmark_tasks = c('system.orig', 'system.new') %>%
  lapply(prepare_tv_task_combined, predict_error = TRUE)

bgrid = benchmark_grid(
    tasks = benchmark_tasks,
    learners = drf_at,
    resamplings = rsmp('block_cv', folds = 5)
)

system.time(bres4 <- benchmark(bgrid))

bres4_summ = bres4$aggregate() %>%
  transform(network = sub('.*_', '', task_id))
get_drf_bench_ahead = function(nr, row_ids) {
  n = length(nr)
  out = rep(integer(), times = n)
  for (i in 1:2) {
    ind_i = nr == i
    row_ids_i = row_ids[ind_i]
    out[ind_i] = benchmark_tasks[[i]]$data(rows = row_ids_i)$days_ahead
  }
  out
}
# summarize drf benchmarking results for plotting
drf_losses = bres4$obs_loss(msr('regr.mae')) %>%
  # have to calculate crps manually for now
  transform(regr.crps = scoringRules::crps_norm(truth, response, se)) %>%
  # get the loss by lead time
  transform(days_ahead = get_drf_bench_ahead(resample_result, row_ids)) %>%
  by(.[, c('days_ahead', 'resample_result')], function(x) {
    network = bres4_summ$network[bres4_summ$nr == x$resample_result[1]]
    data.frame(days_ahead = x$days_ahead[1],
               network = network,
               regr.mae = mean(x$regr.mae),
               regr.crps = mean(x$regr.crps))
  }) %>%
  as.list %>%
  do.call(rbind, .)
write.csv(drf_losses, file = 'results/forecast_tv/drf_losses.csv',
          row.names = F)
