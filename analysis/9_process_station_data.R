# set up the station-based TV data for machine learning training

library(psychrolib)
library(magrittr)
source('R/load_data.R')
source('R/coned_tv.R')

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
  tvs = get_combined_tv(all_sites[[i]])
  names(tvs)[2] = names(all_sites)[i]
  out_tv = merge(out_tv, tvs, by = 'day', all = TRUE)
  eff_tmps = get_combined_eff_tmp(all_sites[[i]])
  names(eff_tmps)[2] = names(all_sites)[i]
  out_eff_tmp = merge(out_eff_tmp, eff_tmps, by = 'day', all = TRUE)
}

write.csv(out_tv, file = 'results/process_station_data/tv.csv', row.names = F)
write.csv(out_eff_tmp, file = 'results/process_station_data/eff_tmp.csv', row.names = F)
