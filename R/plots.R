
library(sf)
library(ggplot2)
#library(mapboxapi)

# run maps.R to get this file
nyc_base = readRDS('results/maps/nyc_basemap.rds')

plot_station_data = function(sites, aes, basemap = nyc_base, size = 8, ...) {
  ggplot() +
    ggspatial::layer_spatial(basemap) +
    geom_sf(aes, sites, size = size, ...) +
    scale_color_viridis_c() +
    scale_shape_manual('Weather network', values = c(15, 19, 17)) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) + 
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), axis.title.x = element_blank(),
          axis.ticks.x = element_blank(), axis.text.x = element_blank(),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) + 
    labs(caption = "© Mapbox, © OpenStreetMap")
}

plot_network_data = function(networks, aes, basemap = nyc_base, alpha = .9,
                             ...) {
  ggplot() +
    ggspatial::layer_spatial(basemap) +
    geom_sf(aes, networks, alpha = alpha, ...) +
    scale_fill_viridis_c() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) + 
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), axis.title.x = element_blank(),
          axis.ticks.x = element_blank(), axis.text.x = element_blank(),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) + 
    labs(caption = "© Mapbox, © OpenStreetMap")
}
