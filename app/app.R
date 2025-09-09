
library(shiny)
library(magrittr)
library(ggplot2)
library(sf)

networks = readRDS('../results/maps/coned_networks_cleaned.rds') %>%
  transform(idnum = as.integer(sub('[A-Z]', '', id)),
            idlet = sub('[0-9]*', '', id)) %>%
  transform(label = paste0(id, ' (', Network, ')'))
networks = networks[with(networks, order(idlet, idnum)), ]

network_choices = c('System (original)' = 'system.orig',
                    'System (UAlbany)' = 'system.new',
                    setNames(networks$id, networks$label))
  
# day_choices = setNames(0:7, networks$label)

forecast_file = list.files('../scripts/forecasts/', full.names = T) %>%
  max
cur_day = forecast_file %>%
  basename %>%
  as.Date(format = 'forecast_tv_%Y_%m%d.csv')
date_note = paste('Forecasts made', cur_day)
preds = read.csv(forecast_file) %>%
  transform(forecast_for = as.Date(forecast_for))

ui = fluidPage(
    titlePanel('TV+Load forecasts | UAlbany Center of Excellence'),
    sidebarLayout(
        sidebarPanel(
            p(date_note),
            selectInput('network', 'Network', network_choices),
            selectInput('day', 'Map day', 0:7)
        ),
        mainPanel(
            fluidRow(plotOutput('plot_ts')),
            fluidRow(plotOutput('map', height = '800px'))
        )
    )
)

server <- function(input, output) {
  output$plot_ts <- renderPlot({
    preds_n = preds[preds$network == input$network, ]
    ts_title = paste('TV Forecast for', input$network)
    y_label = paste(input$network, 'TV')
    ggplot(preds_n, aes(x = forecast_for, y = tv_mean)) +
      geom_ribbon(aes(ymin = tv_lower95, ymax = tv_upper95), fill = "#66666666") +
      geom_line() +
      ylim(min(preds_n$tv_lower95, na.rm = T),
           max(preds_n$tv_upper95, na.rm = T)) +
      xlab('Date') + ylab('TV') +
      ggtitle(ts_title)
  }, res = 120)
  output$map <- renderPlot({
    map_day = cur_day + as.integer(input$day)
    map_title = paste('TV Forecast for', map_day)
    preds_n = preds[preds$days_ahead == input$day, ]
    networks_n = networks %>%
      transform(forecast = preds_n$tv_mean[match(id, preds_n$network)])
    # plot(networks_n[, 'forecast'], main = map_title)
    ggplot(data = networks_n) +
      geom_sf(aes(fill = forecast)) +
      scale_fill_viridis_c() +
      ggtitle(map_title) +
      theme(panel.grid.major = element_blank(), axis.title.x = element_blank(),
            axis.ticks.x = element_blank(), axis.text.x = element_blank(),
            axis.title.y = element_blank(), axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    # for a multi-day version
    # networks_n = merge(networks, preds, by.x = 'id', by.y = 'network')
    # ggplot(data = networks_n) +
    #   geom_sf(aes(fill = tv_mean)) +
    #   scale_fill_viridis_c() +
    #   ggtitle(map_title) +
    #   theme(panel.grid.major = element_blank(), axis.title.x = element_blank(),
    #         axis.ticks.x = element_blank(), axis.text.x = element_blank(),
    #         axis.title.y = element_blank(), axis.text.y = element_blank(),
    #         axis.ticks.y = element_blank()) +
    #   facet_wrap(~ forecast_for, nrow = 2)
  }, res = 120)
}

shinyApp(ui = ui, server = server)
