# Organize weather station data

# create a bounding box based on the network map

library(sf)
library(riem)
source('R/load_data.R')

wgs84 = 4326
li_lcc = 2263

networks = st_read('data/CECONY_Network_Boundary.gdb') %>%
  # ignore westchester
  subset(Borough != 'Westchester')

# project_bbox = networks %>%
#   st_buffer(5280 * 10) %>%
#   st_bbox
# nah I don't like it, box shape doesn't fit the area well

project_area = networks %>%
  st_buffer(5280 * 10) %>%
  st_union
saveRDS(project_area, 'results/station_data/project_area.rds')


# get some ASOS stuff
# asos_nets = riem_networks()

ny_stations = riem_stations(network = 'NY_ASOS')

close_ny = ny_stations %>%
  subset(select = -archive_end) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = wgs84) %>%
  st_transform(li_lcc) %>%
  subset(st_within(geometry, project_area, sparse = F))

nj_stations = riem_stations(network = 'NJ_ASOS')

close_nj = nj_stations %>%
  st_as_sf(coords = c('lon', 'lat'), crs = wgs84) %>%
  st_transform(li_lcc) %>%
  subset(st_within(geometry, project_area, sparse = F))

plot(st_geometry(rbind(close_ny, close_nj)))

# let's show what this looks like

plot(project_area)
plot(st_geometry(networks), border = 'gray', add = T)
plot(st_geometry(rbind(close_ny, close_nj)), pch = 19, add = T)
# looks good!

# select stations within the bounding box

local_asos_stations = paste0('K', c(close_ny$id, close_nj$id))

download_asos_data = function(id) {
  for (year in 2021:2024) {
    out = file.path('data/asos', paste0(year, '_', id, '.csv'))
    # no need to re-download
    if (file.exists(out)) next()
    # start at April to get a head start on the load curve estimates
    start_str = paste0(year, '-04-01')
    end_str = paste0(year, '-10-01')
    obs = riem_measures(id, start_str, date_end = end_str)
    write.csv(obs, file = out, row.names = F)
  }
}

for (id in local_asos_stations) {
  message(id)
  download_asos_data(id)
}


# select NYS Mesonet sites

nysm_file = 'data/mesonet/2023/05/20230501.nc'
nysm_sites = get_nc_sites(nysm_file) %>%
  st_as_sf(coords = c('lon', 'lat'), crs = wgs84) %>%
  st_transform(li_lcc) %>%
  subset(st_within(geometry, project_area, sparse = F))
saveRDS(nysm_sites, 'results/station_data/nysm_sites.rds')

# maybe I should download the nysm data here?


# get micronet sites?


# get NJ mesonet data


# get buoy data, https://www.ndbc.noaa.gov/


# save the site locations

nycm_sites = get_nc_sites('data/micronet/2023/05/20230501.nc') %>%
  transform(network = 'NYC Micronet') %>%
  st_as_sf(coords = c('lon', 'lat'), crs = wgs84) %>%
  st_transform(li_lcc)
# remove MHLSQR, which has no data
nycm_sites = subset(nycm_sites, stid != 'MHLSQR')

nysm_sites$network = 'NYS Mesonet'

asos_sites = rbind(close_ny, close_nj) %>%
  transform(stid = id, network = 'ASOS')

all_sites = rbind(nycm_sites[, c('stid', 'network')],
                  nysm_sites[, c('stid', 'network')],
                  asos_sites[, c('stid', 'network')])
# make indexing easier
row.names(all_sites) = all_sites$stid

# get NYC borough by matching with networks
all_sites$borough = NA
for (i in 1:nrow(all_sites)) {
  idx = st_within(all_sites[i, ], networks, sparse = F)
  if (!any(idx)) next()
  all_sites$borough[i] = networks$Borough[which(idx)[1]]
}
# missed a couple
all_sites$borough[all_sites$stid == 'JRB'] = 'Manhattan'
all_sites$borough[all_sites$stid == 'NYC'] = 'Manhattan'

plot(st_geometry(networks), border = 'gray')
plot(st_geometry(all_sites), pch = 19, col = factor(all_sites$borough), add = T)
# good

saveRDS(all_sites, 'results/station_data/stations.rds')

# make a map?
