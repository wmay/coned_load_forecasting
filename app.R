
library(shiny)
library(magrittr)
library(sf)
library(plotly)
library(leaflet)
library(mapboxapi)

mb_token = readLines('app/mapbox_key.txt')

# set directories and settings based on where this is running
is_local = is.na(Sys.getenv('SHINY_SERVER_VERSION', NA))
if (is_local) {
  forecast_path = 'scripts/forecasts/'
  tv_path = 'scripts/data/tv.csv'
} else {
  # running in shiny server container
  forecast_path = '/mnt/coe/web/coeweather/coned/forecasts'
  tv_path = '/mnt/coe/Will/scheduled/coned_load_forecasting/scripts/data/tv.csv'
}

networks = readRDS('results/maps/coned_networks_cleaned.rds') %>%
  st_transform(4326) %>%
  st_make_valid %>%
  transform(idnum = as.integer(sub('[A-Z]', '', id)),
            idlet = sub('[0-9]*', '', id)) %>%
  transform(label = paste0(id, ' (', Network, ')'))
networks = networks[with(networks, order(idlet, idnum)), ]

network_choices = c('System (official ConEd TV)' = 'system.orig',
                    'System (new UAlbany TV)' = 'system.new',
                    setNames(networks$id, networks$label))

tv_dat = read.csv('scripts/data/tv.csv') %>%
  transform(day = as.Date(day))

forecast_file = list.files(forecast_path, full.names = T) %>%
  max
cur_day = forecast_file %>%
  basename %>%
  as.Date(format = 'forecast_tv_%Y_%m%d.csv')
preds = read.csv(forecast_file) %>%
  transform(forecast_for = as.Date(forecast_for))
load_preds = cur_day %>%
  format('forecast_load_%Y_%m%d.csv') %>%
  file.path(forecast_path, .) %>%
  read.csv %>%
  transform(forecast_for = as.Date(forecast_for))

# map color scale
map_pal = colorNumeric('viridis', range(preds$tv_mean, na.rm = T))
# following https://stackoverflow.com/a/56334156/5548959
map_pal_rev = colorNumeric('viridis', range(-preds$tv_mean, na.rm = T), reverse = T)

# source: https://github.com/rstudio/leaflet/issues/496#issuecomment-650122985
setShapeStyle <- function(map, data = getMapData(map), layerId, stroke = NULL,
                          color = NULL, weight = NULL, opacity = NULL,
                          fill = NULL, fillColor = NULL, fillOpacity = NULL,
                          dashArray = NULL, smoothFactor = NULL, noClip = NULL,
                          options = NULL) {
  options <- c(
      list(layerId = layerId),
      options,
      filterNULL(
          list(stroke = stroke, color = color, weight = weight,
               opacity = opacity, fill = fill, fillColor = fillColor,
               fillOpacity = fillOpacity, dashArray = dashArray,
               smoothFactor = smoothFactor, noClip = noClip)
      )
  )
  # evaluate all options
  options <- evalFormula(options, data = data)
  # make them the same length (by building a data.frame)
  options <- do.call(data.frame, c(options, list(stringsAsFactors=FALSE)))
  layerId <- options[[1]]
  style <- options[-1] # drop layer column
  leaflet::invokeMethod(map, data, "setStyle", "shape", layerId, style)
}

# based on https://github.com/rstudio/leaflet/issues/496#issuecomment-651625559
setShapeLabel <- function(map, data = getMapData(map), layerId, label = NULL,
                          options = labelOptions()) {
  stopifnot(length(layerId) == length(label))
  leaflet::invokeMethod(map, data, "setLabel", "shape", layerId, label, options)
}

date_note = paste('Forecasts made', cur_day)
maintainence_msg = 'Note: Individual network forecasts may not be accurate during plant maintenance or if a network has been recently altered.'
month_warning = 'Warning: Forecasts are only valid May through September.'
cur_month = as.POSIXlt(cur_day)$mon + 1
table_discl = 'Note that individual ConEd networks may have different design criteria thresholds than the thresholds for the ConEd system as a whole (82, 84, 86).'

ui = fluidPage(
    tags$head(tags$link(rel='shortcut icon', href='COE_16px.png'),
              tags$script(src = 'app.js')),
    includeCSS('app/style.css'),
    div(
        titlePanel('NYC TV and Load forecasts | UAlbany Center of Excellence'),
        class = 'page-header'
    ),
    h4(date_note),
    if (cur_month < 5 || cur_month > 9)
      p(month_warning, style = 'color: #E80202;'),
    p(maintainence_msg),
    hr(),
    fluidRow(
        column(6,
               h3('TV forecasts'),
               p('Click the map to select a network.'),
               sliderInput('day', 'Map day', cur_day, cur_day + 7, cur_day,
                           step = 1, ticks = FALSE, animate = TRUE),
               leafletOutput('map', height = 600)
               ),
        column(6,
               h3('Network details'),
               selectInput('network', 'Network', network_choices, width = 400),
               plotlyOutput('plot_ts'),
               tableOutput('table_ts'),
               p(table_discl),
               style = 'border-left: #EEE solid 1px;'
               )
    ),
    hr(),
    includeHTML('app/footer.html')
)

server <- function(input, output) {
  output$plot_ts <- renderPlotly({
    preds_n = preds[preds$network == input$network, ] %>%
      subset(select = c(forecast_for, tv_mean, tv_lower95, tv_upper95))
    load_preds_n = load_preds[load_preds$network == input$network, ]
    if (startsWith(input$network, 'system')) {
      if (input$network == 'system.orig') {
        network_name = 'System (official ConEd TV)'
      } else {
        network_name = 'System (new UAlbany TV)'
      }
      tv_col = input$network
    } else {
      network_name = networks$label[networks$id == input$network]
      tv_col = paste0('network.', input$network)
    }
    ts_title = paste('TV Forecast for', network_name)
    y_label = paste(input$network, 'TV')
    tv_obs = tv_dat[, c('day', tv_col)] %>%
      subset(day < min(preds_n$forecast_for))
    names(tv_obs)[2] = 'tv'

    # connect the forecast to the most recent observation
    preds_n = tv_obs %>%
      subset(day == max(day)) %>%
      transform(forecast_for = day, tv_mean = tv, tv_lower95 = tv, tv_upper95 = tv) %>%
      subset(select = -c(day, tv)) %>%
      rbind(preds_n) %>%
      transform(tv_mean = round(tv_mean, 1),
                tv_lower95 = round(tv_lower95, 1),
                tv_upper95 = round(tv_upper95, 1)) %>%
      transform(hover = paste0(tv_mean, ' (', tv_lower95, ' – ', tv_upper95, ')'))
    preds_n = preds_n[order(preds_n$forecast_for), ]

    plot_bgcolor = '#E6E6E6'
    unwanted_plotly_buttons = c('zoom', 'pan', 'select', 'zoomIn', 'zoomOut',
                                'autoScale', 'resetScale', 'lasso2d',
                                'hoverClosestCartesian', 'hoverCompareCartesian')
    fig1 = plot_ly(preds_n, x = ~forecast_for) %>%
      add_ribbons(ymin = ~tv_lower95, ymax = ~tv_upper95, #hoverinfo = "none",
                  name = '95% Prediction Interval', line = list(width = 0),
                  opacity = 0.6) %>%
      add_trace(y = ~tv_mean,
                # hovertext = ~hover,
                # text = ~hover, hovertemplate = "TV: %{text}",
                type = 'scatter', mode = 'lines',
                name = 'Forecast', line = list(color = '#1f77b4', dash = 'dash')) %>%
      add_trace(x = ~day, y = ~tv, data = tv_obs, type = 'scatter',
                mode = 'lines', name = 'Observation',
                line = list(color = '#1f77b4', dash = 'solid')) %>%
      layout(#title = ts_title,
          xaxis = list(title = 'Day', gridcolor = 'white'),
          yaxis = list(title = 'TV', gridcolor = 'white'),
          plot_bgcolor = plot_bgcolor)

    fig2 = plot_ly(load_preds_n, x = ~forecast_for) %>%
      add_ribbons(ymin = ~load_lower95, ymax = ~load_upper95, #hoverinfo = "none",
                  name = '95% Prediction Interval', line = list(width = 0),
                  # need opacity in the fillcolor to match the default
                  fillcolor = 'rgba(31, 119, 180, 0.5)',
                  opacity = 0.6, showlegend = FALSE) %>%
      add_trace(y = ~load,
                # hovertext = ~hover,
                # text = ~hover, hovertemplate = "TV: %{text}",
                type = 'scatter', mode = 'lines',
                name = 'Forecast', line = list(color = '#1f77b4', dash = 'dash'),
                showlegend = FALSE) %>%
      # add_trace(x = ~day, y = ~tv, data = tv_obs, type = 'scatter',
      #           mode = 'lines', name = 'Observation',
      #           line = list(color = '#1f77b4', dash = 'solid')) %>%
      layout(xaxis = list(title = 'Day', gridcolor = 'white'),
             yaxis = list(title = 'Peak Load (MW)', gridcolor = 'white'),
             plot_bgcolor = plot_bgcolor)

    subplot(fig1, fig2, nrows = 2, shareX = TRUE, titleY = TRUE) %>%
      layout(title = list(text = 'TV and Load Forecasts')) %>%
      config(modeBarButtonsToRemove = unwanted_plotly_buttons,
             displaylogo = F)
  })
  output$table_ts <- renderTable({
    preds_n = preds[preds$network == input$network, ]
    # ts_title = paste('TV Forecast for', input$network)
    # y_label = paste(input$network, 'TV')
    out = preds_n[, -(1:2)] %>%
      transform(ge82 = 100 * pnorm(82, tv_mean, tv_sd, lower.tail = FALSE),
                ge84 = 100 * pnorm(84, tv_mean, tv_sd, lower.tail = FALSE),
                ge86 = 100 * pnorm(86, tv_mean, tv_sd, lower.tail = FALSE)) %>%
      transform(forecast_for = as.character(forecast_for),
                tv_mean = as.character(round(tv_mean, 1)),
                tv_lower95 = as.character(round(tv_lower95, 1)),
                tv_upper95 = as.character(round(tv_upper95, 1))) %>%
      transform(tv_95int = paste0(tv_lower95, ' – ', tv_upper95)) %>%
      subset(select = -c(tv_sd, tv_lower95, tv_upper95))
    # replace tiny percentages
    for (thresh in paste0('ge', c(82, 84, 86))) {
      out[, thresh] = out[, thresh] %>%
        replace(out[, thresh] >= 1, round(out[, thresh])) %>%
        replace(out[, thresh] < 1, '<1')
    }
    names(out) = c('Forecast for', 'Days ahead', 'TV', '>82 (%)', '>84 (%)',
                   '>86 (%)', '95% interval')
    out = out[, c(1:3, 7, 4:6)]
    out
  })
  output$map = renderLeaflet({
    leaflet(networks) %>%
      addMapboxTiles(style_id = 'light-v11', username = 'mapbox',
                     access_token = mb_token,
                     options = tileOptions(tileSize = 512, zoomOffset = -1)) %>%
      # addTiles(
      #     urlTemplate = paste0(
      #         "https://api.mapbox.com/styles/v1/mapbox/light-v11/tiles/{z}/{x}/{y}?access_token=",
      #         mb_token
      #     ),
      #     options = tileOptions(
      #         tileSize   = 512,
      #         zoomOffset = -1
      #     ),
      #     attribution = '© Mapbox © OpenStreetMap'
      # ) %>%
      addPolygons(weight = 1, opacity = 1, color = 'gray',
                  fillColor = '#A0A0A0', fillOpacity = 0.9,
                  layerId = networks$id) %>%
      # following https://stackoverflow.com/a/56334156/5548959
      addLegend("bottomright", pal = map_pal_rev, values = -preds$tv_mean,
                labFormat = labelFormat(transform = function(x) -x),
                title = 'TV')
  })
  # Update the choropleth and date label without recreating the entire leaflet
  # map
  observe({
    req(input$day)
    map_day = input$day
    map_title = paste('Forecast for', map_day)
    title_div = div(map_title, style = 'color: #444444')

    preds_n = preds[preds$forecast_for == input$day, ]
    networks_n = networks %>%
      transform(TV = preds_n$tv_mean[match(id, preds_n$network)])
    network_labels = paste0('Network: ', networks_n$label,
                            '<br>TV forecast: ',
                            round(networks_n$TV, 1)) %>%
      lapply(htmltools::HTML)
    leafletProxy('map') %>%
      removeControl('mapTitle') %>%
      addControl(title_div, 'topleft', layerId = 'mapTitle') %>%
      setShapeStyle(layerId = networks_n$id, fillColor = map_pal(networks_n$TV)) %>%
      setShapeLabel(layerId = networks_n$id, label = network_labels)
  })
  # update the network selection when the map is clicked
  observeEvent(input$map_shape_click, {
    updateSelectInput(inputId = 'network', selected = input$map_shape_click$id)
  })
  # Highlight the selected network on the map
  observe({
    req(input$network)
    leafletProxy('map') %>%
      setShapeStyle(layerId = setdiff(networks$id, input$network),
                    weight = 1, color = 'gray') %>%
      setShapeStyle(layerId = input$network, weight = 3, color = 'white')
  })
}

shinyApp(ui = ui, server = server)
