#' @export
control_explorer_ui <- function(id, measurement_model, process_model) {
  ns <- shiny::NS(id)

  soundings <- measurement_model$soundings
  xco2_range <- range(soundings$xco2)
  xco2_range[1] <- floor(xco2_range[1])
  xco2_range[2] <- ceiling(xco2_range[2])

  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shiny::sliderInput(
        ns('start_date'),
        'Start date:',
        value = date(min(soundings$observation_time)),
        min = date(min(soundings$observation_time)),
        max = date(max(soundings$observation_time)),
        step = 1
      ),
      shiny::sliderInput(
        ns('days_to_include'),
        'Days to include:',
        value = 30,
        min = 1,
        max = 60,
        step = 1
      ),
      shiny::sliderInput(
        ns('xco2_range'),
        'XCO2 range:',
        min = xco2_range[1],
        max = xco2_range[2],
        value = xco2_range,
        step = 0.1
      ),
      width = 3
    ),
    shiny::mainPanel(
      shiny::plotOutput(ns('controlMap')),
      width = 9
    )
  )
}

#' @export
control_explorer <- function(input, output, session, measurement_model, process_model) {
  reactive <- shiny::reactive

  control <- process_model$control

  end_date <- reactive(input$start_date + days(input$days_to_include))
  control_window <- reactive({
    control %>%
      filter(
        time >= input$start_date,
        time <= end_date()
      )
  })

  output$controlMap <- shiny::renderPlot({
    df_average <- control_window() %>%
      group_by(longitude, latitude) %>%
      summarise(xco2_mean = mean(xco2))

    ggplot(df_average, aes(longitude, latitude, fill = xco2_mean)) +
      geom_world() +
      geom_tile() +
      coord_quickmap() +
      scale_fill_wes_palette_c(limits = input$xco2_range) +
      labs(x = 'Longitude', y = 'Latitude', fill = 'XCO2') +
      ggtitle(sprintf(
        'Grid cell averaged XCO2 between %s and %s',
        input$start_date,
        end_date()
      ))
  })
}
