.binary_search_bounds <- function(A, T, start_L = 1, start_R = length(A)) {
  L <- start_L - 1
  R <- start_R
  while (L < R) {
    m <- floor((L + R) / 2)
    if (A[m + 1] < T) {
      L <- m + 1
    } else {
      R <- m
    }
  }
  lower <- L

  L <- start_L - 1
  R <- start_R
  while (L < R) {
    m <- floor((L + R) / 2)
    if (A[m + 1] > T) {
      R <- m
    } else {
      L <- m + 1
    }
  }
  upper <- L - 1

  c(lower, upper) + 1
}

#' @export
sensitivity_explorer_ui <- function(id, process_model) {
  sliderInput <- shiny::sliderInput
  plotOutput <- shiny::plotOutput

  ns <- shiny::NS(id)

  transcoms <- process_model$regions
  names(transcoms) <- sprintf('Transcom %02d', transcoms)

  emissions <- process_model$emissions
  control <- process_model$control
  sensitivities <- process_model$sensitivities

  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shinyWidgets::sliderTextInput(
        ns('month'),
        'Month:',
        choices = format(sort(unique(emissions$month_start)), '%Y-%m'),
        animate = TRUE
      ),
      sliderInput(
        ns('start_date'),
        'Start date:',
        value = date(min(control$time)),
        min = date(min(control$time)),
        max = date(max(control$time)),
        step = 1
      ),
      sliderInput(
        ns('days_to_include'),
        'Days to include:',
        value = 30,
        min = 1,
        max = 60,
        step = 1
      ),
      shiny::selectInput(
        ns('transcom_str'),
        'Transcom:',
        transcoms
      ),
      sliderInput(
        ns('sensitivity_max_abs'),
        'Max absolute sensitivity:',
        value = 1,
        min = 0.01,
        max = max(abs(sensitivities$xco2_sensitivity))
      ),
      sliderInput(
        ns('flux_max_abs'),
        'Max absolute flux:',
        value = 2,
        min = 0.01,
        max = ceiling(max(abs(emissions$flux)))
      ),
      colourInput(ns('low_colour'), 'Low colour:', value = '#35978f'),
      colourInput(ns('high_colour'), 'High colour:', value = '#bf812d'),
      width = 3
    ),
    shiny::mainPanel(
      plotOutput(ns('emissionsMap')),
      plotOutput(ns('sensitivityMap')),
      width = 9
    )
  )
}

#' @export
sensitivity_explorer <- function(
  input, output, session,
  process_model
) {
  reactive <- shiny::reactive
  renderPlot <- shiny::renderPlot

  emissions <- process_model$emissions
  control <- process_model$control
  sensitivities <- process_model$sensitivities

  end_date <- reactive(input$start_date + days(input$days_to_include))
  transcom <- reactive(as.integer(input$transcom_str))

  shiny::observe({
    shiny::updateSliderInput(
      session,
      'start_date',
      min = date(sprintf('%s-01', input$month))
    )
  })

  control_window <- reactive({
    control %>%
      select(
        model_id,
        time,
        longitude,
        latitude
      ) %>%
      filter(
        time >= input$start_date,
        time <= end_date()
      )
  })

  sensitivities_window <- reactive({
    month_date <- date(sprintf('%s-02', input$month))

    transcom_range <- .binary_search_bounds(
      sensitivities$region,
      transcom()
    )
    from_month_start_range <- .binary_search_bounds(
      sensitivities$from_month_start,
      month_date,
      transcom_range[1],
      transcom_range[2]
    )

    control_window() %>%
      inner_join(
        sensitivities[from_month_start_range[1] : from_month_start_range[2], ],
        by = 'model_id'
      )
  })

  output$emissionsMap <- renderPlot({
    ggplot() +
      geom_tile(
        data = emissions %>%
          filter(
            region == transcom(),
            strftime(month_start, '%Y-%m') == input$month
          ),
        mapping = aes(longitude, latitude, fill = flux)
      ) +
      geom_sf(
        data = transcom_boundary,
        fill = NA,
        colour = '#555555',
        size = 0.1
      ) +
      scale_fill_gradient2(
        low = input$low_colour,
        high = input$high_colour,
        limits = c(-1, 1) * input$flux_max_abs,
        oob = scales::squish
      ) +
      labs(
        x = 'Longitude',
        y = 'Latitude',
        fill = expression('Flux [kg/' * m ^ 2 * '/year]')
      ) +
      xlim(-180, 180) +
      ylim(-90, 90)
  })

  output$sensitivityMap <- renderPlot({
    ggplot(
      sensitivities_window(),
      aes(longitude, latitude, colour = xco2_sensitivity)
    ) +
      geom_world() +
      geom_point(
        position = position_jitter(width = 0.5, height = 0.5),
        size = 0.75
      ) +
      coord_quickmap() +
      scale_colour_gradient2(
        low = input$low_colour,
        high = input$high_colour,
        limits = c(-1, 1) * input$sensitivity_max_abs,
        oob = scales::squish
      ) +
      labs(x = 'Longitude', y = 'Latitude', colour = 'XCO2')
  })
}
