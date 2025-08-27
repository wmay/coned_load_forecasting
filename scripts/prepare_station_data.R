# Organize weather station data for daily predictions

library(magrittr)
library(riem)
library(psychrolib)
source('R/load_data.R')
source('R/coned_tv.R')

cur_date = Sys.Date()
nysm_dir = '/run/user/1000/gvfs/sftp:host=swrcc/nysm/archive/nysm/netcdf/proc'
nycm_dir = '/run/user/1000/gvfs/sftp:host=swrcc/nysm/archive/nyc/netcdf/proc'

get_nysm_nc_files = function(network = c('nysm', 'nycm'), cur_date) {
  data_dir = if (network == 'nysm') nysm_dir else nycm_dir
  nysm_files = file.path(data_dir, format(cur_date, '%Y')) %>%
    list.files(full.names = T, recursive = T)
  file_dates =
    as.Date(nysm_files, format = paste0(data_dir, '/%Y/%m/%Y%m%d.nc'))
  out = nysm_files[file_dates >= cur_date - 7]
  if (!length(out)) warning(network, ' files not found')
  out
}

get_combined_nc_data = function(vars = NULL, times = NULL, cur_date) {
  nysm_stids = c("BKLN", "BRON", "MANH", "QUEE", "STAT")
  nycm = get_nysm_nc_files('nycm', cur_date) %>%
    load_nc_files(vars = vars, times = times)
  nysm = get_nysm_nc_files('nysm', cur_date) %>%
    load_nc_files(nysm_stids, vars = vars, times = times)
  rbind(nycm, nysm)
}

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

SetUnitSystem('SI')
tv_vars = c('tair', 'relh')
nysm_obs = get_combined_nc_data(tv_vars, hour_to_nc_index(7:21), cur_date) %>%
  subset(!(is.na(tair) | is.na(relh))) %>%
  transform(dwp = GetTDewPointFromRelHum(tair, relh / 100)) %>%
  transform(tmpf = as_fahrenheit(tair), dwpf = as_fahrenheit(dwp)) %>%
  subset(select = c(time, stid, tmpf, dwpf)) %>%
  reshape(direction = 'wide', idvar = 'time', timevar = 'stid')

stations = readRDS('results/station_data/stations.rds')
asos_obs = stations %>%
  subset(network == 'ASOS') %>%
  getElement('stid') %>%
  paste0('K', .) %>%
  lapply(riem_measures, date_start = cur_date - 7) %>%
  do.call(rbind, .) %>%
  asos_raw_to_hourly %>%
  subset(as.POSIXlt(valid_hour)$hour %in% 7:21) %>%
  transform(stid = station, time = valid_hour) %>%
  subset(!(is.na(tmpf) | is.na(dwpf))) %>%
  subset(select = c(time, stid, tmpf, dwpf)) %>%
  reshape(direction = 'wide', idvar = 'time', timevar = 'stid')

# print warning about missing ASOS station data
asos_desired_stations = stations %>%
  subset(network == 'ASOS') %>%
  getElement('stid')
asos_colnames = names(asos_obs)
missing_stations = asos_colnames[startsWith(asos_colnames, 'tmpf.')] %>%
  substr(6, 8) %>%
  setdiff(asos_desired_stations, .)
if (length(missing_stations)) {
  warning('Missing ASOS observations from ',
          paste(missing_stations, collapse = ', '))
}

station_obs = merge(nysm_obs, asos_obs, all = T)

selected_sites = readRDS('results/select_stations/selected_sites.rds')
system_results = readRDS('results/select_stations/system_results.rds')
system_sites = system_results$stids

all_sites = selected_sites
names(all_sites) = paste0('network.', names(selected_sites))
all_sites[['system.orig']] = c('NYC', 'LGA')
all_sites[['system.new']] = system_results$stids


# For each station collection, we want a file for TV, and a file for effective
# temperature

out_tv = data.frame(day = as.Date(character()))
out_eff_tmp = data.frame(day = as.Date(character()))

for (i in seq_along(all_sites)) {
  network_name = names(all_sites)[i]
  tvs = try(get_combined_tv(all_sites[[i]]))
  if (inherits(tvs, 'try-error')) {
    warning('Missing station data for ', network_name)
    tvs = data.frame(day = cur_date, tv = NA)
  }
  names(tvs)[2] = network_name
  out_tv = merge(out_tv, tvs, by = 'day', all = TRUE)
  eff_tmps = try(get_combined_eff_tmp(all_sites[[i]]))
  if (inherits(eff_tmps, 'try-error')) {
    eff_tmps = data.frame(day = cur_date, eff_tmp = NA)
  }
  names(eff_tmps)[2] = network_name
  out_eff_tmp = merge(out_eff_tmp, eff_tmps, by = 'day', all = TRUE)
}

write.csv(out_tv, file = 'scripts/data/tv.csv', row.names = F)
write.csv(out_eff_tmp, file = 'scripts/data/eff_tmp.csv', row.names = F)
