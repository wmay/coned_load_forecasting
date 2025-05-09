# Run all analyses in order. Must be run from analysis folder

Rscript analysis/1_station_data.R
Rscript analysis/2_maps.R
Rscript analysis/3_weather.R
Rscript analysis/4_energy_loads.R
Rscript analysis/5_load_vs_weather.R
Rscript analysis/6_select_stations.R
