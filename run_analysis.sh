# Run all analyses in order. Must be run from analysis folder

Rscript 1_station_data.R
Rscript 2_maps.R
Rscript 3_weather.R
Rscript 4_energy_loads.R
Rscript 5_load_vs_weather.R
Rscript 6_select_stations.R
