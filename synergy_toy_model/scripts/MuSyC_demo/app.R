#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Old Faithful Geyser Data"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            sliderInput("d2",
                        "Drug 2 dose",
                        min = 0,
                        max = 400,
                        value = 50),
            sliderInput("E0",
                        "E0",
                        min = 0,
                        max = 1,
                        value = .4),
            sliderInput("E1",
                        "E1",
                        min = 0,
                        max = 1,
                        value = .01),
            sliderInput("C1",
                        "Drug 1 IC50",
                        min = 0,
                        max = 400,
                        value = 33),
            sliderInput("s1",
                        "s1",
                        min = -2,
                        max = 5,
                        step = .1,
                        value = 1),
            sliderInput("E2",
                        "E2",
                        min = 0,
                        max = 1,
                        value = .4),
            sliderInput("C2",
                        "Drug 2 IC50",
                        min = 0,
                        max = 400,
                        value = 50),
            sliderInput("s2",
                        "s2",
                        min = -2,
                        max = 5,
                        step = .1,
                        value = 1),
            sliderInput("E3",
                        "E3",
                        min = 0,
                        max = 1,
                        value = .06),
            sliderInput("alpha",
                        "alpha",
                        min = 0,
                        max = 5,
                        value = 1)
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("distPlot"),
           textOutput("describe_h1")
        )
    )
)

infectivity_y_scales <- list(
    ggplot2::scale_y_continuous(
        name = "% Infeected cells per well",
        labels = scales::percent_format(accuracy = 1)))

lf_single_agent_scales <- list(
    ggplot2::scale_x_continuous(
        name = "Log[Drug 1]",
        breaks = c(1.52, 3.12, 6.25, 12.5, 25, 50, 100, 200, 400) %>% log10(),
        labels = c("0", "3.12", "6.25", "12.5", "25", "50", "100", "200", "400")))

generate_MuSyC_effects <- function(
    d1,
    d2,
    E0,
    h1, C1, E1,
    h2, C2, E2,
    alpha,
    E3 = NULL,
    beta = NULL) {
    
    if(!is.null(beta)){
        E3 <- min(E1, E2) - beta * min(E1, E2)
    } else if(is.null(E3)){
        stop("either E3 or beta must be non-null")
    }
    
    if(d1 == -Inf & d2 == -Inf){
        response <- E0
    } else if(d2 == -Inf){
        if(d1 == Inf){
            response <- E1
        } else {
            response <-
                (C1^h1 * E0 + d1^h1 * E1) / 
                (C1^h1      + d1^h1)
        }
    } else if(d1 == -Inf){
        if(d2 == Inf){
            response <- E2
        } else {
            response <-
                (C2^h2 * E0 + d2^h2 * E2) /
                (C2^h2      + d2^h2)
        }
    } else {
        if(d1 == Inf & d2 == Inf) {
            response <- E3
        } else {
            response <-
                (
                    C1^h1 * C2^h2 * E0 +
                    d1^h1 * C2^h2 * E1 +
                    C1^h1 * d2^h2 * E2 +
                    d1^h1 * d2^h2 * E3 * alpha) /
                (
                    C1^h1 * C2^h2 +
                    d1^h1 * C2^h2 +
                    C1^h1 * d2^h2 +
                    d1^h1 * d2^h2 * alpha)
        }
    }
    return(response)
}


# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        data_background <- expand.grid(
            d1 = seq(from = 0, 400, length.out = 200),
            d2 = 10^seq(from = log10(input$C2)-3, log10(input$C2) + 3, length.out = 7)) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(
                Ed = generate_MuSyC_effects(
                    d1 = log10(d1),
                    d2 = log10(d2),
                    E0 = input$E0,
                    h1 = input$s1 * (4 * log10(input$C1) / (input$E0 + input$E1)),
                    C1 = log10(input$C1),
                    E1 = input$E1,
                    h2 = input$s2 * (4 * input$C2) / (input$E0 + input$E2),
                    C2 = log10(input$C2),
                    E2 = input$E2,
                    E3 = input$E3,
                    alpha = input$alpha)) %>%
            dplyr::ungroup()
        
        data_forground <- data.frame(
            d1 = seq(from = 0, 400, length.out = 200),
            d2 = input$d2) %>%
            dplyr::rowwise() %>%
            dplyr::mutate(
                Ed = generate_MuSyC_effects(
                    d1 = log10(d1),
                    d2 = log10(d2),
                    E0 = input$E0,
                    h1 = input$s1 * (4 * log10(input$C1) / (input$E0 + input$E1)),
                    C1 = log10(input$C1),
                    E1 = input$E1,
                    h2 = input$s2 * (4 * input$C2) / (input$E0 + input$E2),
                    C2 = log10(input$C2),
                    E2 = input$E2,
                    E3 = input$E3,
                    alpha = input$alpha)) %>%
            dplyr::ungroup()

        ggplot2::ggplot() +
            ggplot2::theme_bw() +
            ggplot2::geom_line(
                data = data_background,
                mapping = ggplot2::aes(
                    x=log10(d1),
                    y=Ed,
                    group = d2,
                    color = d2),
                size = .8,
                color = "grey80") +
            ggplot2::geom_line(
                data = data_forground,
                mapping = ggplot2::aes(
                    x=log10(d1),
                    y=Ed,
                    group = d2,
                    color = d2),
                size = 1.5) +
            ggplot2::geom_hline(yintercept = input$E0) +
            ggplot2::geom_hline(yintercept = input$E1) +
            ggplot2::geom_vline(xintercept = log10(input$C1)) +
            ggplot2::scale_x_continuous("log[Drug 1]") +
            ggplot2::scale_y_continuous(
                "Percent infected cells per well",
                limits = c(0, .5),
                labels = scales::percent_format(accuracy = 1))
    })
    
    output$describe_h1 <- renderText({ 
        paste0(
            "h1 = s1 * (4 * C1) / (E0 + E1) = ",
            input$s1, " * (4 * ", input$C1, ") / (", input$E0, " + ", input$E1, ") = ",
            input$s1 * (4 * input$C1) / (input$E0 + input$E1), "\n")
    })
}

# Run the application 
shinyApp(ui = ui, server = server)
