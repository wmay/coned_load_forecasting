
library(shiny)
library(magrittr)
library(ggplot2)
library(sf)
library(plotly)

networks = readRDS('../results/maps/coned_networks_cleaned.rds') %>%
  transform(idnum = as.integer(sub('[A-Z]', '', id)),
            idlet = sub('[0-9]*', '', id)) %>%
  transform(label = paste0(id, ' (', Network, ')'))
networks = networks[with(networks, order(idlet, idnum)), ]

network_choices = c('System (official ConEd TV)' = 'system.orig',
                    'System (new UAlbany TV)' = 'system.new',
                    setNames(networks$id, networks$label))

tv_dat = read.csv('../scripts/data/tv.csv') %>%
  transform(day = as.Date(day))

forecast_file = list.files('../scripts/forecasts/', full.names = T) %>%
  max
cur_day = forecast_file %>%
  basename %>%
  as.Date(format = 'forecast_tv_%Y_%m%d.csv')
date_note = paste('Forecasts made', cur_day)
preds = read.csv(forecast_file) %>%
  transform(forecast_for = as.Date(forecast_for))

ui = fluidPage(
    # tags$head(tags$title(page_title),
    #           tags$link(rel='shortcut icon', href='COE_16px.png')),
    tags$head(tags$link(rel='shortcut icon', href='COE_16px.png')),
    includeCSS('style.css'),
    div(
        titlePanel('NYC TV and Load forecasts | UAlbany Center of Excellence'),
        # h3(date_note),
        class = 'page-header'
    ),
    h3(date_note, style = 'font-style: italic;'),
    hr(),
    # wellPanel(
    #     selectInput('network', 'Network:', network_choices, width = 400),
    #     selectInput('day', 'Map day', 0:7),
    #     # uiOutput('categorySelect', inline = T),
    #     # div(textOutput('updateTime'), actionLink('about', 'App documentation.'),
    #     #     style = 'display: inline-block;'),
    #     id = 'headerBar'
    # ),
    # fluidRow(
    #     selectInput('network', 'Network:', network_choices, width = 400),
    #     selectInput('day', 'Map day', 0:7),
    #     # uiOutput('categorySelect', inline = T),
    #     # div(textOutput('updateTime'), actionLink('about', 'App documentation.'),
    #     #     style = 'display: inline-block;'),
    #     id = 'headerBar'
    # ),
    fluidRow(
        column(6,
               # h3('TV forecasts'),
               # sliderInput('day', 'Map day', 0, 7, 0, step = 1, ticks = FALSE, animate = TRUE),
               plotOutput('map', height = '1000px')
               ),
        column(6,
               h3('Network details'),
               selectInput('network', 'Network', network_choices, width = 400),
               plotlyOutput('plot_ts'),
               tableOutput('table_ts'),
               style = 'border-left: #EEE solid 1px;'
               )
    ),
    # fluidRow(
    #     column(6,
    #            plotOutput('plot_ts'),
    #            tableOutput('table_ts')
    #            ),
    #     column(6, plotOutput('map', height = '800px'))
    # ),
    # fluidRow(tableOutput('table_ts')),
        # fluidRow(plotOutput('map', height = '800px'))
    # sidebarLayout(
    #     sidebarPanel(
    #         p(date_note),
    #         # selectInput('network', 'Network', network_choices),
    #         # selectInput('day', 'Map day', 0:7)
    #         fluidRow(tableOutput('table_ts'))
    #     ),
    #     mainPanel(
    #         fluidRow(plotOutput('plot_ts')),
    #         fluidRow(tableOutput('table_ts')),
    #         fluidRow(plotOutput('map', height = '800px'))
    #     )
    # ),
    hr(),
    includeHTML('footer.html')
)

server <- function(input, output) {
  output$plot_ts <- renderPlotly({
    preds_n = preds[preds$network == input$network, ] %>%
      subset(select = c(forecast_for, tv_mean, tv_lower95, tv_upper95))
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

    unwanted_plotly_buttons = c('zoom', 'pan', 'select', 'zoomIn', 'zoomOut',
                                'autoScale', 'resetScale', 'lasso2d',
                                'hoverClosestCartesian', 'hoverCompareCartesian')
    plot_ly(preds_n, x = ~forecast_for) %>%
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
      layout(title = ts_title,
             xaxis = list(title = 'Day'),
             yaxis = list(title = 'TV')) %>%
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
    # print(names(out))
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
  output$map <- renderPlot({
    map_day = cur_day + as.integer(input$day)
    map_title = paste('TV Forecasts for', map_day)
    preds_n = preds[preds$days_ahead == input$day, ]
    networks_n = networks %>%
      transform(TV = preds_n$tv_mean[match(id, preds_n$network)])
    # plot(networks_n[, 'forecast'], main = map_title)
    # ggplot(data = networks_n) +
    #   geom_sf(aes(fill = TV)) +
    #   scale_fill_viridis_c() +
    #   ggtitle(map_title) +
    #   theme(panel.grid.major = element_blank(), axis.title.x = element_blank(),
    #         axis.ticks.x = element_blank(), axis.text.x = element_blank(),
    #         axis.title.y = element_blank(), axis.text.y = element_blank(),
    #         axis.ticks.y = element_blank())
    # for a multi-day version
    networks_n = merge(networks, preds, by.x = 'id', by.y = 'network') %>%
      transform(TV = tv_mean)
    ggplot(data = networks_n) +
      geom_sf(aes(fill = TV)) +
      scale_fill_viridis_c() +
      ggtitle('TV Forecasts') +
      theme(panel.grid.major = element_blank(), axis.title.x = element_blank(),
            axis.ticks.x = element_blank(), axis.text.x = element_blank(),
            axis.title.y = element_blank(), axis.text.y = element_blank(),
            axis.ticks.y = element_blank()) +
      facet_wrap(~ forecast_for, nrow = 3)
  }, res = 120)
}

shinyApp(ui = ui, server = server)
