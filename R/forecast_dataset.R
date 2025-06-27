# preparing the training datasets for forecasting

# get a dataet of mean TV from GEFS
make_nwp_tv_dataset = function(nc, lon, lat, days_ahead = 2) {
  # when read by `read_mdim` time values are read as dates despite having a time
  # as well
  dates = attr(nc, 'dimensions')$time$values + days_ahead
  # days ahead starts at 2 in the array
  if (days_ahead < 2) stop('TV only available starting day 2')
  tv = nc[['eff_temp']][lon, lat,, days_ahead - 1, ] %>%
    colMeans(na.rm = TRUE)
  data.frame(day = dates, tv = tv)
}

# get a matrix of GEFS TV forecasts
get_gefs_samples = function(nc, lon = 3, lat = 4, days_ahead = 2) {
  # when read by `read_mdim` time values are read as dates despite having a time
  # as well
  dates = attr(nc, 'dimensions')$time$values + days_ahead
  # days ahead starts at 2 in the array
  if (days_ahead < 2) stop('TV only available starting day 2')
  nc[['eff_temp']][lon, lat,, days_ahead - 1, ]
  tv = nc[['eff_temp']][lon, lat,, days_ahead - 1, ]
  list(day = dates, tv = t(tv))
}
