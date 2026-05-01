# Get NYISO load data to use as a predictor. This is a substitute for ConEd load
# data, since we can't get it in real time.

library(magrittr)

# download and organize the data

today = Sys.Date()
start_month = as.Date('2021-04-01')
end_month = trunc(today, 'month')

download_month = function(m) {
  m_dir = file.path('data', 'nyiso', format(m, '%Y-%m'))
  tmp_zip = tempfile()
  m_url = paste0('https://mis.nyiso.com/public/csv/pal/', format(m, '%Y%m'),
                 '01pal_csv.zip')
  on.exit(unlink(tmp_zip))
  download.file(m_url, tmp_zip)
  unzip(tmp_zip, exdir = m_dir)
}

months = seq(start_month, end_month, by = 'month')

for (m in as.list(months)) download_month(m)


# ok let's extract the daily peaks

get_daily_peak = function(day) {
  day_file = paste0(format(day, '%Y%m%d'), 'pal.csv')
  day_path = file.path('data', 'nyiso', format(day, '%Y-%m'), day_file)
  readings = try(read.csv(day_path, check.names = FALSE))
  if (inherits(readings, 'try-error')) {
    warning('Failed to read data for ', day)
    return()
  }
  # note that all dates are in daylight savings time (EDT)
  if (any(readings$`Time Zone` != 'EDT')) stop('Non-EDT time zone')
  readings$time = as.POSIXct(readings$`Time Stamp`, tz = 'EST5EDT',
                             format = '%m/%d/%Y %H:%M:%S')
  attr(readings$time, 'tzone') = 'EST'
  readings %>%
    subset(Name == 'N.Y.C.') %>%
    # must be between 9am and 9pm EST
    subset(as.POSIXlt(time)$hour >= 9 & as.POSIXlt(time)$hour < 21) %>%
    transform(date = as.Date(time, tz = 'EST')) %>%
    aggregate(Load ~ date, FUN = max, data = .)
}

get_dates = function(y) {
  start_date = as.Date(paste0(y, '-04-01'))
  end_date = as.Date(paste0(y, '-09-30'))
  seq(start_date, end_date, by='day')
}

load_dates = lapply(2021:2026, get_dates) %>%
  unlist %>%
  as.Date
load_dates = load_dates[load_dates < today]

peaks = lapply(load_dates, get_daily_peak) %>%
  do.call(rbind, .)

write.csv(peaks, file = 'results/get_nyiso_loads/loads.csv', row.names = FALSE)
