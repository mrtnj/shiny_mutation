library(shiny)
library(ggplot2)
library(plyr)

## Package and compiler flags for the C++ simulation function
library(Rcpp)
Sys.setenv( "PKG_CXXFLAGS"="-std=c++11" )
sourceCpp("simulation.cpp")

ui <- fluidPage(
  tags$h1("Simulation of selection, mutation and drift at a single locus"),
  column(4,
         sliderInput(inputId = "log_N",
                     label = "log10(Population size)",
                     value = 2, min = 1, max = 4),
         sliderInput(inputId = "log_mu",
                     label = "log10(Per locus mutation rate)",
                     value = -6, min = -10, max = -3),
         sliderInput(inputId = "log_s",
                     label = "log10(Selection coefficient)",
                     value = -1, min = -4, max = 0),
         sliderInput(inputId = "h",
                     label = "Dominance coefficient",
                     value = 0.5, min = 0, max = 1),
         sliderInput(inputId = "q0",
                     label = "Starting frequency",
                     value = 0.1, min = 0, max = 1),
         actionButton(inputId = "run_button",
                      label = "Run"),
         tags$p("See the code on ",
                tags$a(href = "https://github.com/mrtnj/shiny_mutation", "GitHub"))),
  column(8,
         plotOutput(outputId = "plot"),
         tableOutput(outputId = "endpoint_table"))
)


## Simulation in R

sim_variation <- function(N, mu, s, h, gen, q0) {
  q <- numeric(gen)
  q[1] <- q0
  for (i in 2:gen) {
    geno <- rbinom(n = N, size = 2, prob = q[i - 1])
    fitness <- rep(1, N)
    fitness <- ifelse(geno == 1, 1 - h * s, fitness)
    fitness <- ifelse(geno == 2, 1 - s, fitness)
    survival_roll <- runif(N, 0, 1)
    survival <- ifelse(survival_roll < fitness, TRUE, FALSE)
    surviving <- geno[survival]
    mutation_roll <- runif(length(surviving), 0, 1)
    mutated <- surviving
    mutated[which(mutation_roll < mu & surviving == 0)] <- 1
    mutated[which(mutation_roll < mu & surviving == 1)] <- 2
    mutation_roll2 <- runif(length(surviving), 0, 1)
    mutated[which(mutation_roll2 < mu & mutated == 0)] <- 1
    mutated[which(mutation_roll2 < mu & mutated == 1)] <- 2
    n_P = length(which(surviving == 0))
    n_H = length(which(surviving == 1))
    n_Q = length(which(surviving == 2))
    q[i] <- (n_H + n_Q * 2) / (2 * (n_P + n_H + n_Q))
  }
  q
}

## Haldane's expression for equilibrium frequency

q_eq <- function(mu, k1, k2) {
  (k1 + mu - sqrt((k1 - mu)^2 + 4 * k2 * mu)) /
    (2 * k1 - 2 * k2)
}


## The plot

plot_simulations <- function(sim) {
  combined_sim <- ldply(1:length(sim), function(i) {
    transform(sim[[i]], replicate = i)
  })
  qplot(x = gen, y = q, group = replicate, data = combined_sim,
        geom = "line", alpha = I(0.5)) + theme_bw() + ylim(0, 1) +
        xlab("Generation") + ylab("Variant frequency")
}

server <- function(input, output) {
  simulated_data <- reactiveValues(data = NULL)
  
  observeEvent(input$run_button, {
    N <- 10^input$log_N
    mu <- 10^input$log_mu
    s <- 10^input$log_s
    h <- input$h
    q0 <- input$q0
    
    q_eq <- q_eq(mu, s * input$h, s)
    
    ## If you'd prefer the simulation function in R, this is where you would
    ## remove the "_cpp" on the function call.
    sim <- replicate(10,
      data.frame(gen = 1:200,
      q = sim_variation_cpp(N = N,
                           mu = mu,
                           s = s,
                           h = h,
                           q0 = q0,
                           gen = 200)),
      simplify = FALSE)
      
    endpoints <- ldply(1:10, function(i) data.frame(final_frequency = sim[[i]]$q[200],
                                                    replicate = i))
    
    frequency_plot <- plot_simulations(sim) + 
      geom_hline(yintercept = q_eq, colour = "red") +
      annotate("text", x = 200, y = q_eq + 0.1, label = signif(q_eq, 2), colour = "red")  
    
    simulated_data$data <- list(frequency_plot, endpoints)
  })
  output$plot <- renderPlot({
    if (! is.null(simulated_data$data))
      print(simulated_data$data[[1]])
  })
  output$endpoint_table <- renderTable({
    if (! is.null(simulated_data$data))
      simulated_data$data[[2]]
  })
}

shinyApp(ui = ui, server = server)
