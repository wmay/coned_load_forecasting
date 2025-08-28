# run some ML distributional regression models for forecasting TV

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
library(mlr3)
source('R/mlr3_distr.R')

# see notes in `?benchmark` and `?mlr_tuners_random_search`
lgr::get_logger("mlr3")$set_threshold("warn")
lgr::get_logger("bbotk")$set_threshold("warn")

all_eff_tmps = read.csv('results/process_station_data/eff_tmp.csv')
all_tvs = read.csv('results/process_station_data/tv.csv')

# if coordinates are not in the same order for every variable, `read_mdim` works
# fine while `read_ncdf` messes up the coordinates
gefs_3day = read_ncdf('results/process_nwp_data/gefs_3day_wmean.nc')


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

prepare_newdata = function(network, days_ahead = 0) {
  # Note: netcdf files use forecast day, while `system_tv` uses valid day. Must
  # be consistent!
  is_system = startsWith(network, 'system')
  extract_name = if (is_system) 'system' else network
  out = extract_values(extract_name, days_ahead) %>%
    subset(select = -eff_temp) # 3-day eff_temp is same as TV
  # what the hell is up with the different column names?
  param_dict = c(avg_ishf = 'msshf', avg_slhtf = 'mslhf')
  changed_names = which(names(out) %in% names(param_dict))
  if (length(changed_names)) {
    names(out)[changed_names] = param_dict[names(out)[changed_names]]
  }
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
  # Predict the GEFS forecast error instead of directly predicting TV
  out %>%
    subset(select = -day) %>%
    na.omit %>%
    transform(fct_err = as.numeric(NA)) %>%
    as_task_regr('fct_err', id = id)
}

prepare_newdata_combined = function(network) {
  task = prepare_newdata(network, 0)
  task$cbind(data = data.frame(days_ahead = rep(0, task$nrow)))
  # append the other lead times
  for (i in 1:7) {
    task_i = prepare_newdata(network, i)
    task_i$cbind(data = data.frame(days_ahead = rep(i, task_i$nrow)))
    task$rbind(data = task_i$data())
  }
  task$id = paste('TV', network, sep = '_')
  # make sure resampling is grouped by day_ahead
  task$set_col_roles('days_ahead', add_to = 'stratum')
  task
}

get_predictions = function(network) {
  task = prepare_newdata_combined(network)
  pred = tv_models[[network]]$predict(task)
  # because we predicted the forecast error, now re-add the forecast
  pred_tv = pred$response + task$data()$TV
  data.frame(network = network, days_ahead = task$data()$days_ahead,
             tv_mean = pred_tv, tv_sd = pred$se)
}


tv_models = readRDS('results/forecast_tv/tv_models.rds')
# this is a large file-- 10GB?

cur_date = as.Date(attr(gefs_3day, 'dimensions')$time$values)
out_file = format(cur_date, '%Y_%m%d') %>%
  paste0('forecast_tv_', ., '.csv') %>%
  file.path('scripts', 'forecasts', .)

preds = names(tv_models) %>%
  lapply(get_predictions) %>%
  do.call(rbind, .) %>%
  transform(forecast_made = cur_date) %>%
  transform(forecast_for = forecast_made + days_ahead) %>%
  subset(select = c(network, forecast_made, forecast_for, days_ahead, tv_mean,
                    tv_sd)) %>%
  transform(tv_mean = round(tv_mean, 2), tv_sd = round(tv_sd, 2))

write.csv(preds, file = out_file, row.names = FALSE)


# library(ggplot2)
# # library(tidyr)

# gefs_baseline = all_res %>%
#   subset(learner_id == 'GEFS') %>%
#   transform(days_ahead = as.integer(substr(task_id, 13, 13))) %>%
#   with(setNames(regr.crps, days_ahead))

# all_res %>%
#   transform(days_ahead = as.integer(substr(task_id, 13, 13))) %>%
#   ggplot(aes(x = days_ahead, y = regr.crps, color = learner_id)) +
#   geom_line()

# all_res %>%
#   transform(days_ahead = as.integer(substr(task_id, 13, 13))) %>%
#   transform(crps = regr.crps / (gefs_baseline[as.character(days_ahead)])) %>%
#   ggplot(aes(x = days_ahead, y = crps, color = learner_id)) +
#   geom_line()
