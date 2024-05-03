
library(ncdf4)
library(magrittr)

# functions to organize mesonet data from netcdf

get_nc_vars = function(x) {
  unname(sapply(x$var, getElement, name = 'name'))
}

# just the site info
get_nc_sites = function(nc_path, sites = NULL) {
  nc = nc_open(nc_path)
  on.exit(nc_close(nc))
  if (!is.null(sites)) stop('site selection not implemented yet')
  all_vars = get_nc_vars(nc)
  var_dims = sapply(nc$var, function(x) length(x$dim))
  out_vars = all_vars[var_dims == 1]
  stids = nc$dim[[1]]$vals
  out = data.frame(stid = stids)
  for (v in out_vars) out[, v] = ncvar_get(nc, v)
  out
}

# extract measurements
load_nc_file = function(nc_path, sites = NULL, vars = NULL, times = NULL) {
  nc = nc_open(nc_path)
  on.exit(nc_close(nc))
  nc_stids = nc$dim[[1]]$vals
  if (is.null(sites)) {
    out_sites = 1:length(nc_stids)
  } else {
    out_sites = which(nc_stids %in% sites)
    if (!length(out_sites)) stop('sites not found')
  }
  nc_times = as.POSIXct(nc$dim[[3]]$vals, origin = '1970-01-01')
  attr(nc_times, 'tzone') = 'EST5EDT'
  if (is.null(times)) {
    out_times = 1:length(nc_times)
  } else {
    out_times = times
  }
  if (is.null(vars)) {
    all_vars = get_nc_vars(nc)
    var_dims = sapply(nc$var, function(x) length(x$dim))
    out_vars = all_vars[var_dims == 2]
  } else {
    out_vars = vars
  }
  stids = nc_stids[out_sites]
  dat_list = lapply(out_vars, function(x) {
    res = as.data.frame(ncvar_get(nc, x)[out_times, out_sites, drop = FALSE])
    names(res) = paste(x, stids, sep = '.')
    res
  })
  # could also get the dot-separated names by naming the list
  # names(dat_list) = out_vars
  df1 = do.call(cbind, dat_list)
  reshape(df1, direction = 'long', varying = 1:ncol(df1), sep = '.',
          ids = nc_times[out_times], timevar = 'stid', idvar = 'time')
}

load_nc_files = function(nc_paths, sites = NULL, vars = NULL, times = NULL) {
  data_list = lapply(nc_paths, load_nc_file, sites = sites, vars = vars,
                     times = times)
  do.call(rbind, data_list)
}

# `x` is the hour in EDT
hour_to_nc_index = function(edt_hour) {
  utc_hour = sort((edt_hour + 4) %% 24)
  1 + utc_hour * 12 # indexing starts at 1, w/ 5-minute increments
}

# easy way to get the full dataset
get_combined_nc_data = function(vars = NULL, times = NULL) {
  nysm_stids = c("BKLN", "BRON", "MANH", "QUEE", "STAT")
  nycm = list.files('data/micronet', full.names = T, recursive = T) %>%
    load_nc_files(vars = vars, times = times)
  nysm = list.files('data/mesonet', full.names = T, recursive = T) %>%
    load_nc_files(nysm_stids, vars = vars, times = times)
  rbind(nycm, nysm)
}

read_asos_file = function(f) {
  out = read.csv(f)
  # times are omitted at midnight, which screws up parsing
  out$valid = sub('^([^:]*)$', '\\1 00:00:00', out$valid)
  out$valid = as.POSIXct(out$valid, tz = 'UTC')
  out
}

# get the fahrenheit measurements closest to the hour
asos_raw_to_hourly = function(dat) {
  dat = subset(dat, !is.na(tmpf)) %>%
    transform(valid_hour = as.POSIXct(round(valid, 'hours')),
              minute = as.POSIXlt(valid)$min) %>%
    # careful here-- this time zone change only works when `valid_hour` is
    # already POSIXct
    transform(valid_hour = as.POSIXct(valid_hour, tz = 'EST5EDT')) %>%
    # minutes from nearest hour
    transform(from_hour = ifelse(minute > 30, 60 - minute, minute))
  dat[with(dat, order(valid_hour, from_hour)), ] %>%
    subset(!duplicated(valid_hour), select = -c(minute, from_hour))
}

get_hourly_asos_data = function(edt_hours = NULL) {
  list.files('data/asos', full.names = T) %>%
    lapply(function(f) {
      out = asos_raw_to_hourly(read_asos_file(f))
      if (!is.null(edt_hours)) {
        out = subset(out, as.POSIXlt(valid_hour)$hour %in% edt_hours)
      }
      out
    }) %>%
    do.call(rbind, .)
}
