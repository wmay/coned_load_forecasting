# exploratory analysis of load data

library(magrittr)
library(sf)
library(mapboxapi)
library(cowplot)
source('R/plots.R')

wgs84 = 4326
web_merc = 3857

loads = list.files('data/coned', 'Network.*\\.csv', full.names = T) %>%
  lapply(read.csv) %>%
  do.call(rbind, .) %>%
  transform(DT = as.POSIXct(DT, tz = 'EST5EDT', '%m/%d/%Y %H:%M:%S'),
            # some inconsistency, but various forms of 'False' mean data is fine
            BAD = !startsWith(BAD, 'F')) %>%
  subset(!is.na(DT))
names(loads) = tolower(names(loads))
  
# B=Brooklyn, M=Manhattan, Q=Queens, R=Staten Island (Richmond County), X=Bronx

loads %>%
  transform(year = trunc(dt, 'years')) %>%
  subset(select = c(year, bad)) %>%
  table
#             bad
# year          FALSE   TRUE
#   2021-01-01 256940    100
#   2022-01-01 256408    632
#   2023-01-01 255440   1600

# with(load23, plot(DT, READING))

# get peak readings

loads %>%
  subset(!bad) %>%
  transform(day = as.Date(dt)) %>%
  subset(select = c(dt, day)) %>%
  head(30)

peaks = loads %>%
  subset(!bad) %>%
  transform(day = as.Date(dt, tz = 'EST5EDT')) %>%
  aggregate(reading ~ network + day, ., max)

# peak times of day
peak_times = loads %>%
  subset(!bad) %>%
  # aggregate loads
  aggregate(reading ~ dt, ., sum) %>%
  # important to get the time zone correct here
  transform(day = as.Date(dt, tz = 'EST5EDT')) %>%
  { .[with(., order(day, reading, decreasing = T)), ] } %>%
  # get first entry (highest load) for each day
  { .[!duplicated(.$day), ] } %>%
  transform(hour = as.POSIXlt(dt)$hour)
table(peak_times$hour) # mode at 4pm
times_hist = hist(peak_times$hour, 24, right = F)
saveRDS(times_hist, 'results/energy_loads/max_load_times.rds')



# peak by network

networks = readRDS('results/maps/coned_networks_cleaned.rds') %>%
  st_transform(web_merc)

network_peaks = loads %>%
  subset(!bad) %>%
  aggregate(reading ~ network, ., max)
networks$peak = network_peaks$reading[match(networks$id, network_peaks$network)]
# plot(networks[, 'peak'])

# save 2021 peak values, to use for weighting later
loads %>%
  subset(!bad & dt < '2022-01-01') %>%
  aggregate(reading ~ network, ., max) %>%
  write.csv(file = 'results/energy_loads/network_peaks_2021.csv',
            row.names = FALSE)

# plot_station_data(networks, aes(fill = peak), nyc_base, alpha = .9) +
#   scale_fill_viridis_c(option = 'magma')

# hourly peaks by borough
borough_peaks = loads %>%
  subset(!bad) %>%
  aggregate(reading ~ borough + dt, ., sum, na.rm = T) %>%
  aggregate(reading ~ borough, ., max, na.rm = T)
#         borough   reading
# 1         Bronx 1065.0876
# 2      Brooklyn 2569.9881
# 3     Manhattan 4184.9484
# 4        Queens 1830.9877
# 5 Staten Island  624.7099


# peak load densities

#networks$peak_density = with(networks, peak / Shape_Area)
networks$peak_density = with(networks, peak / (Shape_Area / 5280^2))
saveRDS(networks, 'results/energy_loads/network_peaks.rds')
# plot(networks[, 'peak_density'])

# plot_station_data(networks, aes(fill = peak_density), nyc_base, alpha = .9) +
#   scale_fill_viridis_c(trans = 'log10', option = 'plasma')

p3 = plot_station_data(networks, aes(fill = peak), nyc_base, alpha = .9) +
  scale_fill_viridis_c(option = 'magma')

p4 = plot_station_data(networks, aes(fill = peak_density), nyc_base, alpha = .9) +
  scale_fill_viridis_c(trans = 'log10', option = 'magma')

png('results/energy_loads/peak_loads.png', width = 1400, height = 700, res = 100)
plot_grid(p3, p4, ncol = 2, nrow = 1, align = 'v')
dev.off()


# peak by weekday/weekend/holiday



# mds

peaks_wide = reshape(peaks, direction = 'wide', timevar = 'network',
                     idvar = 'day')
peak_cors = cor(peaks_wide[, -1], use = 'pairwise.complete.obs')
colnames(peak_cors) = sub('.*\\.', '', colnames(peak_cors))
row.names(peak_cors) = sub('.*\\.', '', row.names(peak_cors))

peak_coords = cmdscale(sqrt(1 - peak_cors))
plot(peak_coords)
text(peak_coords, labels = row.names(peak_cors))

# compare with peak value

networks[, c('mds1', 'mds2')] = peak_coords[networks$id, ]

# with(networks, plot(mds1, peak))
# # nah it doesn't really look like coords represent size

# plot_station_data(networks, aes(fill = mds1), nyc_base, alpha = .9) +
#   scale_fill_viridis_c()
# # this looks almost exactly like the R^2 plot

with(networks, plot(mds1, log(peak_density)))
# it's reverse peak density!

# what's the second one?
plot_station_data(networks, aes(fill = mds2), nyc_base, alpha = .9) +
  scale_fill_viridis_c()
# just a few weird sites?
