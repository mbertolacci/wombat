#' @export
emission_inventory_explorer_ui <- function(id, process_model) {
  ns <- shiny::NS(id)

  emissions <- process_model$emissions

  shiny::sidebarLayout(
    shiny::sidebarPanel(
      shinyWidgets::sliderTextInput(
        ns('month'),
        'Month:',
        choices = format(sort(unique(emissions$month_start)), '%Y-%m'),
        animate = TRUE
      ),
      shiny::sliderInput(
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
      shiny::plotOutput(ns('emissionsMap')),
      width = 9
    )
  )
}

#' @export
emission_inventory_explorer <- function(
  input,
  output,
  session,
  process_model,
  transcom_boundary
) {
  emissions <- process_model$emissions
  output$emissionsMap <- shiny::renderPlot({
    ggplot() +
      geom_tile(
        data = emissions %>%
          filter(format(month_start, '%Y-%m') == input$month),
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
        limits = c(-input$flux_max_abs, input$flux_max_abs),
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
}
