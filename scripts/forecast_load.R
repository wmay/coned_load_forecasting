# Evaluate daily peak load forecasting methods

# pkg_deps = c('ncmeta', 'scoringRules')
# install.packages(pkg_deps)
# remotes::install_github("mlr-org/mlr3temporal")
# distr6_repos =  c(CRAN = 'https://cloud.r-project.org',
#                   raphaels1 = 'https://raphaels1.r-universe.dev')
# install.packages('distr6', repos = distr6_repos)
# install.packages("mlr3proba", repos = "https://mlr-org.r-universe.dev")
# library(timeDate) # holidays
library(magrittr)
library(stars) # also requires ncmeta
library(mlr3)
library(torch)
# source('R/mlr3_distr.R')

# set directories and settings based on where this is running
if (nzchar(Sys.getenv('SUPERCRONIC'))) {
  # running in supercronic container
  out_dir = '/mnt/coe/web/coeweather/coned/forecasts'
} else {
  # running locally
  out_dir = 'scripts/forecasts'
}

# see notes in `?benchmark` and `?mlr_tuners_random_search`
lgr::get_logger("mlr3")$set_threshold("warn")
lgr::get_logger("bbotk")$set_threshold("warn")

all_eff_tmps = read.csv('scripts/data/eff_tmp.csv')

# if coordinates are not in the same order for every variable, `read_mdim` works
# fine while `read_ncdf` messes up the coordinates
gefs_3day = read_ncdf('scripts/data/gefs_3day_wmean.nc')

# get values from stars object at the centroid of a network
networks = readRDS('results/maps/coned_networks_cleaned.rds')
# add centroid coordinates
network_coords = networks %>%
  st_centroid %>%
  st_coordinates %>%
  as.data.frame
row.names(network_coords) = networks$id
names(network_coords) = c('nad83_x', 'nad83_y')
# stations = readRDS('results/station_data/stations.rds')

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
  eff_tmps_col =
    if (network == 'system') 'system.new' else paste0('network.', network)
  system_eff_tmp = all_eff_tmps[, c('day', eff_tmps_col)]
  names(system_eff_tmp)[2] = 'eff_temp'
  if (days_ahead == 1) {
    tv + .1 * system_eff_tmp$eff_temp[match(day - 2, system_eff_tmp$day)]
  } else if (days_ahead == 0) {
    tv + .2 * system_eff_tmp$eff_temp[match(day - 1, system_eff_tmp$day)] +
      .1 * system_eff_tmp$eff_temp[match(day - 2, system_eff_tmp$day)]
  }
}

prepare_load_task = function(network, days_ahead) {
  # Note: netcdf files use forecast day, while `system_tv` uses valid day. Must
  # be consistent!
  is_system = startsWith(network, 'system')
  extract_name = if (is_system) 'system' else network
  out = extract_values(extract_name, days_ahead) %>%
    subset(select = -eff_temp) %>% # 3-day eff_temp is same as TV
    # the time measure, to allow for a nonstationary model
    transform(fct_run = as.integer(valid_day - days_ahead))
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
    out$TV = with(out, fill_incomplete_tv(TV, valid_day, network, days_ahead))
  }
  out$doy = as.POSIXlt(out$valid_day)$yday
  # print(out)
  # get load data
  load_col = paste0('reading.', network)
  id = paste('load', network, days_ahead, sep = '_')
  out %>%
    # transform(load = all_loads[match(valid_day, all_loads$day), load_col]) %>%
    # na.omit %>% # does that make sense?
    # transform(load = as.numeric(NA)) %>%
    transform(load = 0) %>% # we aren't using the load number anyway
    subset(select = -valid_day) %>%
    as_task_regr('load', id = id)
}

# create a task that combines lead times, as in https://doi.org/10.1002/qj.4701
prepare_load_task_combined = function(network) {
  task = prepare_load_task(network, 0)
  task$cbind(data = data.frame(days_ahead = rep(0, task$nrow)))
  # append the other lead times
  for (i in 1:7) {
    task_i = prepare_load_task(network, i)
    task_i$cbind(data = data.frame(days_ahead = rep(i, task_i$nrow)))
    task$rbind(data = task_i$data())
  }
  task$id = paste('load', network, sep = '_')
  # make sure resampling is grouped by day_ahead
  task$set_col_roles('days_ahead', add_to = 'stratum')
  task
}

prepare_load_task_all_networks = function(network_ids = networks$id) {
  network_tasks = network_ids %>%
    lapply(prepare_load_task_combined) %>%
    lapply(function(x) {
      # add the network as a predictor
      network_id = sub('load_', '', x$id)
      # remove soilw since some networks don't have it
      if ('soilw' %in% x$feature_names) {
        x$select(setdiff(x$feature_names, 'soilw'))
      }
      x$cbind(data.frame(network = rep(network_id, x$nrow)))
    })
  out_task = network_tasks[[1]]
  for (i in 2:length(network_tasks)) {
    out_task$rbind(network_tasks[[i]]$data())
  }
  # add NAD83 coordinates as predictors
  out_task$cbind(network_coords[out_task$data()$network, ])
  # convert network column to factor
  network_fct = factor(out_task$data()$network, levels = network_ids)
  out_task$select(out_task$feature_names[out_task$feature_names != 'network'])
  out_task$cbind(data.frame(network = network_fct))
  out_task$id = paste('network_loads')
  out_task
}


network_model = readRDS('results/forecast_load/network_model.rds')
# network_model$unmarshal()

fct_task = prepare_load_task_all_networks()
pred = network_model$predict(fct_task)

cur_date = attr(gefs_3day, 'dimensions')$time$offset %>%
  as.Date
out_file = format(cur_date, '%Y_%m%d') %>%
  paste0('forecast_load_', ., '.csv') %>%
  file.path(out_dir, .)

out = data.frame(network = fct_task$data()$network,
                 forecast_for = cur_date + fct_task$data()$days_ahead,
                 days_ahead = fct_task$data()$days_ahead,
                 load = pred$response,
                 load_se = pred$se) %>%
  transform(load_lower95 = load - 1.96 * load_se,
            load_upper95 = load + 1.96 * load_se) %>%
  transform(load = round(load, 1), load_se = round(load_se, 1),
            load_lower95 = round(load_lower95, 1),
            load_upper95 = round(load_upper95, 1))

write.csv(out, file = out_file, row.names = FALSE)
message('Wrote load forecast file')
