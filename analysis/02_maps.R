# Make reference maps of stations, networks, and study area

library(sf)
library(ggmap)
library(ggrepel)
library(cowplot)
library(mapboxapi)
source('R/load_data.R')

wgs84 = 4326
web_merc = 3857
mb_token = Sys.getenv('mapbox_key')

match_first_word = function(x, y) {
  # also remove a pesky apostrophe
  match(gsub(" .*|'", '', x), gsub(" .*|'", '', y))
}

# Get the corresponding network IDs from load data. `networks` is the map/sf
# collection, `loads` is the load csv data. It's a pain due to some naming
# inconsistencies
get_network_ids = function(networks, loads) {
  network_dict = loads %>%
    subset(network != '', select = c(network.name, network, borough)) %>%
    unique
  # quick check to make sure we don't have duplicated names, which would cause
  # issues
  name_counts = table(table(network_dict$network.name))
  if (length(name_counts) > 1) stop("Found duplicated names")
  ids = network_dict$network[match(networks$Network, network_dict$network.name)]
  # fix some values with slightly differing names by matching on only the first
  # word of the name
  missing_networks = subset(network_dict, !network %in% ids)
  missing_id = is.na(ids)
  ids[missing_id] =
    missing_networks$network[match_first_word(networks$Network[missing_id],
                                              missing_networks$network.name)]
  # remove false matches identified through incorrect borough
  false_match = !is.na(ids) &
    networks$Borough != network_dict$borough[match(ids, network_dict$network)]
  # this should only catch Richmond Hill in Queens
  if (sum(false_match) != 1)
    stop('Expected 1 borough mismatch, found ', sum(false_match))
  ids[false_match] = NA
  ids
}

# get nice basemap from mapbox
networks = st_read('data/CECONY_Network_Boundary.gdb') %>%
  # No westchester in the load data. SN (sub network?) networks aren't included
  # in the load data
  subset(Borough != 'Westchester' & !grepl('SN', Network))


nyc_base = get_static_tiles(st_union(st_make_valid(networks)), zoom = 10,
                            style_id = "light-v9", username = "mapbox",
                            access_token = mb_token)
saveRDS(nyc_base, 'results/maps/nyc_basemap.rds')

# and now the plots code works
source('R/plots.R')


# station reference map
stations = readRDS('results/station_data/stations.rds') %>%
  st_transform(web_merc)

png('results/maps/station_reference.png', width = 1300, height = 1000, res = 100)
plot_station_data(stations, aes(color = network, shape = network), size = 3,
                  alpha = .7) +
  scale_color_discrete() +
  scale_shape_manual(values = c(15, 19, 17)) +
  geom_text_repel(aes(label = stid, geometry = geometry), stations,
                  stat = "sf_coordinates")
dev.off()


# Clean up network data, and add network IDs to match the load data
loads = read.csv('data/coned/Network Data - SUNY 2023.csv') %>%
  transform(DT = as.POSIXct(DT, tz = 'EST5EDT', '%m/%d/%Y %H:%M:%S'),
            # some inconsistency, but various forms of 'False' mean data is fine
            BAD = !startsWith(BAD, 'F')) %>%
  subset(!is.na(DT))
names(loads) = tolower(names(loads))

networks$id = get_network_ids(networks, loads)
# fix Richmond Hill
richmond_hill_ind = which(startsWith(networks$Network, 'Richmond Hill'))
# plot(st_union(networks[richmond_hill_ind, ]))
# # borough boundary lines remain
richmond_hill_combined = st_union(networks$Shape[richmond_hill_ind])
richmond_hill_fixed_coords = st_coordinates(richmond_hill_combined)[1:73, 1:3]
networks$Shape[richmond_hill_ind[2]] = st_polygon(list(richmond_hill_fixed_coords))
# fix the area too (length is still wrong)
networks$Shape_Area[richmond_hill_ind[2]] = sum(networks$Shape_Area[richmond_hill_ind])
networks = networks[-richmond_hill_ind[1], ]
row.names(networks) = networks$id

# save this cleaned up data
saveRDS(networks, 'results/maps/coned_networks_cleaned.rds')

# save the network centroids to use for NWP interpolation later
networks %>%
  st_centroid %>%
  st_transform(wgs84) %>%
  transform(lon = st_coordinates(.)[, 'X'],
            lat = st_coordinates(.)[, 'Y']) %>%
  as.data.frame %>%
  subset(select = c(id, lon, lat)) %>%
  write.csv(file = 'results/maps/network_centroids.csv', row.names = FALSE)


# network reference maps
networks = st_transform(networks, web_merc)

# will have to be two plots-- get a bounding box for a separate map of southern
# manhattan
manh_bbox = st_bbox(subset(networks, Borough == 'Manhattan'))
manh_bbox[4] = mean(manh_bbox[c(2, 4)]) # remove northern manhattan
networks$south_manhattan =
  t(st_contains(st_as_sfc(manh_bbox), networks, sparse = F))
networks$south_manhattan =
  with(networks, south_manhattan & Borough == 'Manhattan')
# refine the bounding box again
manh_bbox = st_bbox(subset(networks, south_manhattan))

manh_basemap = subset(networks, south_manhattan) %>%
  st_make_valid %>%
  st_union %>%
  get_static_tiles(zoom = 12, style_id = "light-v9", username = "mapbox",
                   access_token = mb_token)
saveRDS(manh_basemap, 'results/maps/manh_basemap.rds')

p1 = plot_station_data(networks, aes(), alpha = 0, color = '#888888') +
  geom_text(aes(label = id, geometry = Shape),
            networks[!networks$south_manhattan, ], stat = "sf_coordinates")
p2 = subset(networks, south_manhattan) %>%
  plot_station_data(aes(), manh_basemap, alpha = 0, color = '#909090',
                    linewidth = .3) +
  geom_text(aes(label = id, geometry = Shape),
            networks[networks$south_manhattan, ], stat = "sf_coordinates") +
  xlim(manh_bbox[c(1, 3)]) +
  ylim(manh_bbox[c(2, 4)])

png('results/maps/network_reference.png', 1550, 750, res = 100)
plot_grid(p1, p2, align = 'v')
dev.off()


# Micronet+ConEd map

stadia_key = Sys.getenv('stadia_key')
register_stadiamaps(stadia_key)
has_stadiamaps_key() # yes

project_area = readRDS('results/station_data/project_area.rds') %>%
  st_transform(web_merc)

# simple version
plot(project_area)
plot(st_geometry(networks), border = 'gray', add = T)
plot(st_geometry(stations), pch = 19, add = T)

# for this part I want to use the original network boundaries, with networks
# split between boroughs
networks = st_read('data/CECONY_Network_Boundary.gdb') %>%
  # No westchester in the load data. SN (sub network?) networks aren't included
  # in the load data
  subset(Borough != 'Westchester' & !grepl('SN', Network)) %>%
  st_transform(web_merc)

# get bounding box in lon/lat
gg_bbox = st_transform(networks, wgs84) %>%
  st_bbox %>%
  { make_bbox(.[c(1, 3)], .[c(2, 4)], f = .1) }

stamen_basemap = get_stadiamap(gg_bbox, 11, 'stamen_terrain')
ggmap(stamen_basemap) # yay
saveRDS(stamen_basemap, 'results/maps/stamen_basemap.rds')

network_pch = c(ASOS = 22, `NYS Mesonet` = 24, `NYC Micronet` = 23)
network_col = c(ASOS = 'darkred', `NYS Mesonet` = 'darkred',
                `NYC Micronet` = 'darkorange3')
network_bg = c(ASOS = 'red', `NYS Mesonet` = 'red', `NYC Micronet` = 'yellow')

# NOTE: to avoid basemap pixelation, resize plot to original basemap size.
# Margins are removed to simplify the math

# Calculate x- and y-ranges. Then size of cropped basemap. Has to be in lon/lat
# for ggmap. 
plot_bbox = st_transform(networks, wgs84) %>%
  st_bbox %>%
  # expand by 4%, same as base R plot
  { make_bbox(.[c(1, 3)], .[c(2, 4)], f = .04) }

# rows are y, columns are x
bm_x = ncol(stamen_basemap) * ((plot_bbox[3] - plot_bbox[1]) / (gg_bbox[3] - gg_bbox[1]))
bm_y = nrow(stamen_basemap) * ((plot_bbox[4] - plot_bbox[2]) / (gg_bbox[4] - gg_bbox[2]))
# the aspect ratio will be constant, so it should be OK to vary one dimension as
# long as the constraining/binding dimension is correct

orig_mar = par('mar')
png('results/maps/micronet.png', width = bm_x, height = bm_y, res = 100)
par(mar = rep(0, 4))
palette('Set 2')
palette(paste0(head(palette(), -1), 'B0'))
plot(st_geometry(networks), border = '#4A4A4A', col = factor(networks$Borough),
     bgMap = stamen_basemap)
plot(st_geometry(stations), pch = network_pch[stations$network], cex = 1.5,
     col = network_col[stations$network], bg = network_bg[stations$network],
     add = T)
# add borough names
boroughs = networks %>%
  st_make_valid %>%
  aggregate(list(networks$Borough), head, n = 1) %>%
  st_centroid
borough_centers = st_coordinates(boroughs)
text(borough_centers[, 1], borough_centers[, 2], boroughs$Borough, cex = 1.8, col = '#333333B0')
# st_bbox(networks) # the extent of the plot
#     xmin     ymin     xmax     ymax 
# -8266186  4938199 -8194055  5063192
# st_bbox(subset(networks, Borough != 'Westchester'))
#     xmin     ymin     xmax     ymax 
# -8266186  4938199 -8204217  4999970 
legend(-8215217, 4947199, c('ASOS', 'NYS Mesonet', 'NYC Micronet'),
       col = c('darkred', 'darkred', 'darkorange3'), bg = '#FFFFFF',
       pch = c(22, 24, 23), pt.bg = c('red', 'red', 'yellow'), pt.cex = 1.5,
       title = 'Station Types', title.cex = 1.1)
dev.off()
par(mar = orig_mar)
