# making TV forecasts for Brownsville (network 3B)

setwd('..')
library(psychrolib)
library(magrittr)
library(stars) # also requires ncmeta
# install.packages('ncmeta')
# library(crch)
# library(scoringRules)
# library(beepr)
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
# brownsville = networks[3, ] %>%
#   st_centroid %>%
#   st_transform(wgs84) %>%
#   merge(data.frame(time = get_valid_day(gefs_0p25_3day, 0)))
# gefs_3b = gefs_0p25_3day %>%
#     dplyr::slice('edt9pm_day', 3) %>%
#     st_extract(brownsville, bilinear = TRUE, time_column = 'time')

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
    out = out_dat %>%
      subset(select = -day) %>%
      as_task_regr('fct_err', id = id)
    list(task = out, dates = out_dat$day)
  } else {
    out %>%
      transform(tv_obs = system_tv$tv[match(day, system_tv$day)]) %>%
      subset(select = -day) %>%
      na.omit %>%
      as_task_regr('tv_obs', id = id)
  }
}
# st_extract for interpolation


# # set up the effective temp data
# data_wide = readRDS('results/load_vs_weather/tv_and_load.rds')

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


selected_sites = readRDS('results/select_stations/selected_sites.rds')
system_tv = get_combined_tv(selected_sites$'3B') %>%
  # not training on 2025 data
  subset(day < '2025-01-01')

# predict daily value (effective temp) instead of TV
# system_eff_tmp = get_combined_eff_tmp(system_tv_sites)



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
    term_evals = 10,
    store_tuning_instance = FALSE
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
    tuner = tnr('random_search'),
    learner = drf, 
    resampling = rsmp('block_cv'),
    measure = msr('regr.crps'),
    term_evals = 10,
    store_tuning_instance = FALSE
)

# run the inner loop in parallel and the outer loop sequentially
future::plan(list('sequential', 'multicore'), workers = 4)
# future::plan("sequential")

# benchmark methods, first predicting the result directly, then predicting the
# GEFS forecast error  
forecast_tv_tasks = lapply(0:7, prepare_tv_task, predict_error = FALSE)
bgrid = benchmark_grid(
    tasks = forecast_tv_tasks,
    learners = ngr,
    # learners = c(ngr, rangerdist_at, drf_at),
    resamplings = rsmp('forecast_holdout')
)
system.time(bres <- progressr::with_progress(benchmark(bgrid)))
bres$aggregate(c(msr('regr.crps'), msr('regr.mae'), msr('regr.rmse')))
saveRDS(bres, 'benchmarks.rds')

forecast_tv_tasks2 = lapply(0:7, prepare_tv_task, predict_error = TRUE)
forecast_tv_tasks3 = lapply(0:7, prepare_tv_task, predict_error = TRUE)
bgrid = benchmark_grid(
    tasks = forecast_tv_tasks2,
    # learners = ngr,
    learners = c(ngr, rangerdist_at, drf_at),
    resamplings = rsmp('forecast_holdout')
)
system.time(bres2 <- progressr::with_progress(benchmark(bgrid)))
bres2$aggregate(c(msr('regr.crps'), msr('regr.mae'), msr('regr.rmse')))

# does one stand out as best?
bres2$aggregate(c(msr('regr.crps'), msr('regr.mae'), msr('regr.rmse'))) %>%
  aggregate(regr.crps ~ learner_id, FUN = mean, data = .)
#              learner_id regr.crps
# 1        regr.drf.tuned  1.038836
# 2              regr.ngr  1.057980
# 3 regr.rangerdist.tuned  1.042399
# it's all so close

# train drf_at on the datasets, and make predictions for summer 2025
# run the tuning in parallel
future::plan('multicore', workers = 4)

models = lapply(0:7, function(days_ahead) {
  print(days_ahead)
  drf_at = auto_tuner(
      tuner = tnr('random_search'),
      learner = drf, 
      resampling = rsmp('block_cv'),
      measure = msr('regr.crps'),
      term_evals = 10,
      store_tuning_instance = FALSE
  )
  progressr::with_progress(drf_at$train(forecast_tv_tasks2[[days_ahead + 1]]))
  drf_at
})
saveRDS(models, 'results/brownsville/models.rds')

# now forecast 2025 values

gefs_0p25_3day = read_nwp_nc('results/process_nwp_data/2025_gefs_0p25_3day_wmean.nc')
gefs_0p5_3day = read_nwp_nc('results/process_nwp_data/2025_gefs_0p5_3day_wmean.nc')
gefs_tv = read_nwp_nc('results/process_nwp_data/2025_gefs_tv.nc')
system_tv = get_combined_tv(selected_sites$'3B') %>%
  # not training on 2025 data
  subset(day > '2025-01-01')

forecast_tv_2025 = lapply(0:7, prepare_tv_task, predict_error = TRUE)

pred_list = lapply(0:7, function(i) {
  models[[i + 1]]$predict(forecast_tv_2025[[i + 1]]$task)
})


# save forecasts
preds = 0:7 %>%
  lapply(function(i) {
    pred_err = unlist(pred_list[[i + 1]]$distr$getParameterValue('mean'))
    gefs_fct = forecast_tv_2025[[i + 1]]$task$data()$TV
    mean = unlist(pred_list[[i + 1]]$distr$getParameterValue('mean'))
    # the 12pm from the netcdf files shows up here as 0.5!
    dates = forecast_tv_2025[[i + 1]]$dates %>%
      as.integer %>%
      as.Date
    data.frame(date = dates - i,
               days_ahead = i,
               truth = gefs_fct - pred_list[[i + 1]]$truth,
               pred_mean = gefs_fct - pred_err,
               pred_sd = unlist(pred_list[[i + 1]]$distr$getParameterValue('sd'))) %>%
      # fill missing forecasts
      merge(data.frame(days_ahead = 0:7))
  }) %>%
  do.call(rbind, .) %>%
  transform(lower_95ci = pred_mean - 1.96 * pred_sd,
            upper_95ci = pred_mean + 1.96 * pred_sd) %>%
  transform(forecast_made = date, forecast_for = date + days_ahead) %>%
  subset(select = c('forecast_made', 'forecast_for', 'days_ahead', 'truth',
                    'pred_mean', 'pred_sd', 'lower_95ci', 'upper_95ci'))

preds = preds[order(preds$forecast_made, preds$days_ahead), ]
for (v in c('pred_mean', 'pred_sd', 'lower_95ci', 'upper_95ci')) preds[, v] = round(preds[, v], 1)

write.csv(preds, 'brownsville_forecasts.csv', row.names = FALSE)


# make plots of true and forecast TVs
library(ggplot2)

system_tv3 = as.list(preds$forecast_made) %>%
  lapply(function(x) {
    system_tv %>%
      subset(day >= x & day <= x + 7) %>%
      transform(forecast_made = x,
                days_ahead = as.integer(day - x))
  }) %>%
  do.call(rbind, .)

png('brownsville_forecasts.png', width = 1500, height = 1000, res = 120)
ggplot(preds, aes(x = forecast_for, y = pred_mean)) +
  geom_ribbon(aes(ymin = lower_95ci, ymax = upper_95ci), fill = "#66666666") +
  geom_line() +
  geom_line(aes(x = day, y = tv), system_tv3, color = 'red') +
  ylim(min(preds$lower_95ci), max(preds$upper_95ci)) +
  xlab('Date') + ylab('Brownsville-specific TV (network 3B)') +
  facet_wrap(~forecast_made, scales = 'free_x')
dev.off()
