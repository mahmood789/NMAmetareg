#  Network Meta-Analysis and Meta-Regression Shiny App

# Load required packages
library(shiny)
library(shinydashboard)
library(shinyBS)
library(shinyjs)
library(netmeta)
library(meta)
library(DT)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(readxl)
library(readr)
library(metafor)

# Create sample data files for different effect measures
create_sample_data <- function(type = "binary") {
  set.seed(42)  # For reproducibility
  
  if(type == "binary") {
    # Binary outcome data (OR, RR)
    binary_data <- data.frame(
      Study = paste0("Study", sprintf("%02d", 1:15)),
      Treatment1 = c(rep("Placebo", 5), "TreatmentA", "TreatmentA", "TreatmentB", 
                     "TreatmentB", "TreatmentC", rep("Placebo", 5)),
      Treatment2 = c("TreatmentA", "TreatmentB", "TreatmentC", "TreatmentD", "TreatmentE",
                     "TreatmentB", "TreatmentC", "TreatmentC", "TreatmentD", "TreatmentD",
                     "TreatmentA", "TreatmentB", "TreatmentC", "TreatmentD", "TreatmentE"),
      Events1 = round(runif(15, 10, 50)),
      Total1 = round(runif(15, 80, 150)),
      Events2 = round(runif(15, 15, 60)),
      Total2 = round(runif(15, 80, 150)),
      Year = sample(2000:2020, 15, replace = TRUE),
      Age = round(runif(15, 40, 75), 1),
      Duration = round(runif(15, 6, 24), 1),
      Quality = sample(1:5, 15, replace = TRUE)
    )
    
    return(binary_data)
  }
  else if(type == "continuous") {
    # Continuous outcome data (MD, SMD)
    continuous_data <- data.frame(
      Study = paste0("Study", sprintf("%02d", 1:15)),
      Treatment1 = c(rep("Placebo", 5), "TreatmentA", "TreatmentA", "TreatmentB", 
                     "TreatmentB", "TreatmentC", rep("Placebo", 5)),
      Treatment2 = c("TreatmentA", "TreatmentB", "TreatmentC", "TreatmentD", "TreatmentE",
                     "TreatmentB", "TreatmentC", "TreatmentC", "TreatmentD", "TreatmentD",
                     "TreatmentA", "TreatmentB", "TreatmentC", "TreatmentD", "TreatmentE"),
      Mean1 = round(runif(15, 20, 40), 1),
      SD1 = round(runif(15, 5, 15), 1),
      N1 = round(runif(15, 80, 150)),
      Mean2 = round(runif(15, 15, 50), 1),
      SD2 = round(runif(15, 5, 15), 1),
      N2 = round(runif(15, 80, 150)),
      Year = sample(2000:2020, 15, replace = TRUE),
      Age = round(runif(15, 40, 75), 1),
      Duration = round(runif(15, 6, 24), 1),
      Quality = sample(1:5, 15, replace = TRUE)
    )
    
    return(continuous_data)
  }
  else if(type == "hazard") {
    # Hazard ratio data (HR)
    log_hrs <- rnorm(15, mean = -0.2, sd = 0.3)
    hrs <- exp(log_hrs)
    
    sample_sizes <- round(runif(15, 80, 250))
    se_log_hrs <- sqrt(4/sample_sizes)
    
    hazard_data <- data.frame(
      Study = paste0("Study", sprintf("%02d", 1:15)),
      Treatment1 = c(rep("Placebo", 5), "TreatmentA", "TreatmentA", "TreatmentB", 
                     "TreatmentB", "TreatmentC", rep("Placebo", 5)),
      Treatment2 = c("TreatmentA", "TreatmentB", "TreatmentC", "TreatmentD", "TreatmentE",
                     "TreatmentB", "TreatmentC", "TreatmentC", "TreatmentD", "TreatmentD",
                     "TreatmentA", "TreatmentB", "TreatmentC", "TreatmentD", "TreatmentE"),
      HR = round(hrs, 2),
      LowerCI = round(exp(log_hrs - 1.96 * se_log_hrs), 2),
      UpperCI = round(exp(log_hrs + 1.96 * se_log_hrs), 2),
      SampleSize = sample_sizes,
      Year = sample(2000:2020, 15, replace = TRUE),
      Age = round(runif(15, 40, 75), 1),
      Duration = round(runif(15, 6, 24), 1),
      Quality = sample(1:5, 15, replace = TRUE)
    )
    
    return(hazard_data)
  }
}

# Helper function to calculate meta-regression predictions at specific covariate values
calculate_metareg_predictions <- function(model, covariate_values, treatments) {
  predictions <- list()
  
  # Create a sequence of covariate values if only min and max provided
  if (length(covariate_values) == 2) {
    covariate_values <- seq(covariate_values[1], covariate_values[2], length.out = 10)
  }
  
  # Check if we're using metafor model
  if ("rma" %in% class(model)) {
    for (cov_val in covariate_values) {
      # Create appropriate newmods based on model structure
      newmods <- matrix(cov_val, ncol = 1)
      
      # Get predictions with confidence intervals
      tryCatch({
        pred <- predict(model, newmods = newmods, level = 0.95)
        
        # Store predictions
        predictions[[as.character(cov_val)]] <- data.frame(
          covariate_value = cov_val,
          effect = pred$pred,
          lower_ci = pred$ci.lb,
          upper_ci = pred$ci.ub
        )
      }, error = function(e) {
        # Handle prediction errors
        predictions[[as.character(cov_val)]] <- data.frame(
          covariate_value = cov_val,
          effect = NA,
          lower_ci = NA,
          upper_ci = NA
        )
      })
    }
  } else {
    # Simple approach for lm models
    for (cov_val in covariate_values) {
      # Calculate predicted value
      pred_val <- coef(model)[1] + coef(model)[2] * cov_val
      
      # Calculate confidence interval using predict
      newdata <- data.frame(matched_covs = cov_val)
      pred_interval <- predict(model, newdata = newdata, interval = "confidence", level = 0.95)
      
      # Store prediction
      predictions[[as.character(cov_val)]] <- data.frame(
        covariate_value = cov_val,
        effect = pred_interval[1, "fit"],
        lower_ci = pred_interval[1, "lwr"],
        upper_ci = pred_interval[1, "upr"]
      )
    }
  }
  
  return(predictions)
}

# Helper function to create a bubble plot for meta-regression
create_metareg_bubble_plot <- function(data, xvar, yvar, sizevar, title, xlab, ylab) {
  # Check for required variables
  req_vars <- c(xvar, yvar, sizevar)
  if (!all(req_vars %in% names(data))) {
    return(ggplot() + 
             ggtitle("Required variables not found in data") +
             theme_minimal())
  }
  
  # Check if we have enough data points
  if (nrow(data) < 2) {
    return(ggplot() + 
             ggtitle("Not enough data points for bubble plot") +
             theme_minimal())
  }
  
  # Calculate study weights (inverse of standard error squared)
  weights <- 1 / (data[[sizevar]] ^ 2)
  data$weights <- weights
  
  # Create bubble plot
  p <- ggplot(data, aes_string(x = xvar, y = yvar)) +
    geom_point(aes(size = weights), alpha = 0.6, color = "blue") +
    geom_smooth(aes(weight = weights), method = "lm", se = TRUE, color = "red") +
    labs(title = title, x = xlab, y = ylab, size = "Weight") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    geom_hline(yintercept = 0, linetype = "dashed")
  
  return(p)
}

# Function to create a forest plot with meta-regression adjustments
create_metareg_forest_plot <- function(nma, metareg_model, covariate_value, 
                                       reference_treatment, sm, cov_name) {
  # Get treatments from NMA
  treatments <- nma$trts
  treatments <- treatments[treatments != reference_treatment]
  
  # Extract NMA estimates for comparison with reference
  effects <- c()
  lower_ci <- c()
  upper_ci <- c()
  
  for (treat in treatments) {
    # Get the network estimate
    idx <- which(nma$treat1 == reference_treatment & nma$treat2 == treat)
    if (length(idx) == 0) {
      idx <- which(nma$treat2 == reference_treatment & nma$treat1 == treat)
      if (length(idx) > 0) {
        effects <- c(effects, -nma$TE.nma.random[idx])
        lower_ci <- c(lower_ci, -nma$upper.nma.random[idx])
        upper_ci <- c(upper_ci, -nma$lower.nma.random[idx])
      }
    } else {
      effects <- c(effects, nma$TE.nma.random[idx])
      lower_ci <- c(lower_ci, nma$lower.nma.random[idx])
      upper_ci <- c(upper_ci, nma$upper.nma.random[idx])
    }
  }
  
  # Adjust for covariate if meta-regression model is provided
  if (!is.null(metareg_model) && "rma" %in% class(metareg_model)) {
    # Simple adjustment based on covariate value
    adjustment <- metareg_model$beta[2] * (covariate_value - mean(metareg_model$mods[,1], na.rm = TRUE))
    effects <- effects + adjustment
    # Note: CI adjustment would be more complex in practice
  }
  
  # Create data frame for forest plot
  forest_data <- data.frame(
    treatment = treatments,
    effect = effects,
    lower = lower_ci,
    upper = upper_ci
  )
  
  # Sort by effect size
  forest_data <- forest_data[order(forest_data$effect),]
  
  # Create forest plot
  par(mar = c(5, 8, 4, 2))
  plot(forest_data$effect, 1:nrow(forest_data), 
       xlim = range(c(forest_data$lower, forest_data$upper), na.rm = TRUE),
       ylim = c(0.5, nrow(forest_data) + 0.5),
       xlab = paste("Effect (", sm, ")", sep = ""),
       ylab = "",
       yaxt = "n",
       pch = 15,
       cex = 1.5,
       main = paste("Forest Plot at", cov_name, "=", round(covariate_value, 2),
                    "\nCompared with", reference_treatment))
  
  # Add confidence intervals
  for (i in 1:nrow(forest_data)) {
    lines(c(forest_data$lower[i], forest_data$upper[i]), 
          c(i, i), lwd = 2)
  }
  
  # Add treatment labels
  axis(2, at = 1:nrow(forest_data), labels = forest_data$treatment, las = 1)
  
  # Add vertical line at null effect
  abline(v = 0, lty = 2)
  
  # Add grid
  grid(col = "lightgray", lty = "dotted")
}

# Function to generate treatment rankings adjusted for covariate values
generate_adjusted_rankings <- function(nma, metareg_model, covariate_value, cov_name) {
  # Get treatments from NMA
  treatments <- nma$trts
  n_treatments <- length(treatments)
  
  # Calculate P-scores from network estimates
  # Get all pairwise comparisons
  comparison_matrix <- matrix(NA, n_treatments, n_treatments,
                              dimnames = list(treatments, treatments))
  
  for (i in 1:n_treatments) {
    for (j in 1:n_treatments) {
      if (i != j) {
        # Find the comparison
        idx <- which((nma$treat1 == treatments[i] & nma$treat2 == treatments[j]) |
                     (nma$treat1 == treatments[j] & nma$treat2 == treatments[i]))
        
        if (length(idx) > 0) {
          if (nma$treat1[idx] == treatments[i]) {
            comparison_matrix[i, j] <- nma$TE.nma.random[idx]
          } else {
            comparison_matrix[i, j] <- -nma$TE.nma.random[idx]
          }
        }
      }
    }
  }
  
  # Adjust comparisons based on meta-regression if available
  if (!is.null(metareg_model) && "rma" %in% class(metareg_model)) {
    # Apply covariate adjustment
    mean_cov <- mean(metareg_model$mods[,1], na.rm = TRUE)
    adjustment <- metareg_model$beta[2] * (covariate_value - mean_cov)
    comparison_matrix <- comparison_matrix + adjustment
  }
  
  # Calculate P-scores
  p_scores <- rep(0, n_treatments)
  names(p_scores) <- treatments
  
  for (i in 1:n_treatments) {
    # Count how many treatments this treatment beats
    wins <- 0
    for (j in 1:n_treatments) {
      if (i != j && !is.na(comparison_matrix[i, j])) {
        # For smaller values being better
        if (nma$small.values == "good") {
          wins <- wins + pnorm(0, mean = comparison_matrix[i, j], sd = 1)
        } else {
          wins <- wins + pnorm(0, mean = -comparison_matrix[i, j], sd = 1)
        }
      }
    }
    p_scores[i] <- wins / (n_treatments - 1)
  }
  
  # Sort treatments by P-score
  sorted_idx <- order(p_scores, decreasing = TRUE)
  
  return(list(
    treatments = treatments[sorted_idx],
    p_scores = p_scores[sorted_idx]
  ))
}

# UI
ui <- dashboardPage(
  dashboardHeader(title = "786-MIII NMA Metaregression"),
  
  dashboardSidebar(
    sidebarMenu(id = "sidebar",
                menuItem("Home", tabName = "home", icon = icon("home")),
                menuItem("1. Data", tabName = "data", icon = icon("table")),
                menuItem("2. Analysis", tabName = "analysis", icon = icon("cogs")),
                menuItem("3. NMA Results", tabName = "results", icon = icon("chart-bar")),
                menuItem("4. Meta-regression", icon = icon("chart-line"),
                         menuSubItem("4a. Summary", tabName = "metareg_summary"),
                         menuSubItem("4b. Baseline Risk", tabName = "baseline_risk"),
                         menuSubItem("4c. Covariate Analysis", tabName = "covariate_analysis")
                )
    ),
    
    # Treatment selection preference
    conditionalPanel(
      condition = "input.sidebar == 'analysis' || input.sidebar == 'results' || input.sidebar.indexOf('metareg') !== -1",
      div(
        style = "padding: 15px; background-color: #f8f9fa; margin-top: 15px; border-radius: 5px;",
        p("For treatment rankings, smaller outcome values are:"),
        radioButtons("values_direction", NULL,
                     choices = c("Desirable" = "good", 
                                 "Undesirable" = "bad"),
                     selected = "good", inline = TRUE),
        checkboxInput("random_effects", "Random Effects Model", TRUE)
      )
    )
  ),
  
  dashboardBody(
    useShinyjs(),
    
    # Custom CSS for MetaInsight-like styling
    tags$style(HTML("
      .nav-tabs {
        background-color: #f8f9fa;
        border-radius: 4px 4px 0 0;
      }
      .nav-tabs-custom .nav-tabs li.active {
        border-top-color: #3c8dbc;
      }
      .box-primary {
        border-top-color: #3c8dbc;
      }
      .meta-section-label {
        font-weight: bold;
        color: #3c8dbc;
        margin-bottom: 15px;
      }
      .beta-label {
        font-size: 0.8em;
        color: white;
        background-color: #f39c12;
        padding: 2px 5px;
        border-radius: 3px;
        margin-left: 5px;
      }
      .sample-data-buttons {
        display: flex;
        justify-content: center;
        gap: 10px;
        margin-top: 15px;
        margin-bottom: 15px;
      }
      .download-button {
        margin-top: 10px;
      }
      .help-icon {
        color: #3c8dbc;
        margin-left: 5px;
        cursor: pointer;
      }
      .info-box {
        background-color: #f8f9fa;
        border-left: 4px solid #3c8dbc;
        padding: 10px;
        margin-bottom: 15px;
        border-radius: 3px;
      }
    ")),
    
    tabItems(
      # Home Tab
      tabItem(
        tabName = "home",
        fluidRow(
          box(
            width = 12,
            title = "Network Meta-Analysis with Meta-Regression",
            status = "primary",
            solidHeader = TRUE,
            p("This application provides an interface for conducting network meta-analysis with meta-regression."),
            br(),
            # Sample data options
            div(
              class = "text-center",
              h4("Load Sample Data"),
              div(
                class = "sample-data-buttons",
                actionButton("load_binary_sample", "Binary Outcomes (OR/RR)", 
                             icon = icon("table"), class = "btn btn-success"),
                actionButton("load_continuous_sample", "Continuous Outcomes (MD/SMD)", 
                             icon = icon("table"), class = "btn btn-success"),
                actionButton("load_hazard_sample", "Hazard Ratio (HR)", 
                             icon = icon("table"), class = "btn btn-success")
              ),
              p("Choose the type of sample data to load based on your effect measure.")
            ),
            # Download sample data
            div(
              class = "text-center",
              h4("Download Sample Data"),
              div(
                class = "sample-data-buttons",
                downloadButton("download_binary", "Binary Data (CSV)", class = "btn-info"),
                downloadButton("download_continuous", "Continuous Data (CSV)", class = "btn-info"),
                downloadButton("download_hazard", "Hazard Data (CSV)", class = "btn-info")
              ),
              p("Download sample datasets to use in your own analyses.")
            ),
            br(),
            p("Please follow these steps to perform your analysis:"),
            tags$ol(
              tags$li("Upload your data or use the sample dataset"),
              tags$li("Configure your analysis settings"),
              tags$li("View NMA results"),
              tags$li("Perform meta-regression to examine study-level covariates")
            ),
            
            # Add citation and version information
            div(
              class = "info-box",
              p(strong("Version:"), "1.2.0"),
              p(strong("Last Updated:"), "March 15, 2025"),
              p(strong("Based on:"), "netmeta and metafor R packages")
            )
          )
        ),
        fluidRow(
          box(
            width = 6,
            title = "Effect Measures",
            status = "info",
            p("The app supports the following effect measures:"),
            tags$ul(
              tags$li("Odds Ratio (OR)"),
              tags$li("Risk Ratio (RR)"),
              tags$li("Mean Difference (MD)"),
              tags$li("Standardized Mean Difference (SMD)"),
              tags$li("Hazard Ratio (HR)")
            )
          ),
          box(
            width = 6,
            title = "Meta-Regression Features",
            status = "info",
            p("Meta-regression allows you to explore:"),
            tags$ul(
              tags$li("Effect of baseline risk on treatment effectiveness"),
              tags$li("Impact of study-level covariates on treatment effects"),
              tags$li("Treatment rankings based on covariate values")
            )
          )
        )
      ),
      
      # Data Tab
      tabItem(
        tabName = "data",
        fluidRow(
          box(
            width = 12,
            title = "1. Data Input",
            status = "primary",
            fluidRow(
              column(6,
                     fileInput("file", "Upload CSV or Excel File",
                               accept = c(".csv", ".xls", ".xlsx")),
                     checkboxInput("header", "File has header", TRUE),
                     selectInput("sep", "Separator (for CSV)", 
                                 choices = c(Comma = ",", Semicolon = ";", Tab = "\t"),
                                 selected = ",")
              ),
              column(6,
                     # Sample data buttons
                     div(
                       style = "margin-bottom: 15px;",
                       actionButton("load_binary_sample_data", "Load Binary Data (OR/RR)", 
                                    icon = icon("table"), class = "btn btn-info btn-block"),
                       actionButton("load_continuous_sample_data", "Load Continuous Data (MD/SMD)", 
                                    icon = icon("table"), class = "btn btn-info btn-block"),
                       actionButton("load_hazard_sample_data", "Load Hazard Ratio Data (HR)", 
                                    icon = icon("table"), class = "btn btn-info btn-block")
                     ),
                     actionButton("proceed_to_analysis", "Proceed to Analysis →", 
                                  class = "btn btn-success btn-block")
              )
            )
          )
        ),
        fluidRow(
          box(
            width = 12,
            title = "Data Preview",
            status = "primary",
            DTOutput("data_preview"),
            div(
              class = "download-button text-center",
              downloadButton("download_current_data", "Download Current Data", class = "btn-primary")
            )
          )
        ),
        # Add data validation section
        fluidRow(
          box(
            width = 12,
            title = "Data Validation",
            status = "primary",
            collapsible = TRUE,
            collapsed = TRUE,
            verbatimTextOutput("data_validation_output")
          )
        )
      ),
      
      # Analysis Tab
      tabItem(
        tabName = "analysis",
        fluidRow(
          box(
            width = 12,
            title = "2. Analysis Settings",
            status = "primary"
          )
        ),
        fluidRow(
          box(
            width = 6,
            title = "Effect Measure",
            status = "primary",
            selectInput("effect_measure", "Effect Measure",
                        choices = c("Odds Ratio" = "OR", 
                                    "Risk Ratio" = "RR",
                                    "Mean Difference" = "MD",
                                    "Standardized Mean Difference" = "SMD",
                                    "Hazard Ratio" = "HR"),
                        selected = "OR")
          ),
          box(
            width = 6,
            title = "Study Selection",
            status = "primary",
            uiOutput("study_selection_ui")
          )
        ),
        fluidRow(
          box(
            width = 12,
            title = "Column Mapping",
            status = "primary",
            fluidRow(
              column(4,
                     selectInput("study_col", "Study Column", choices = NULL),
                     selectInput("treat1_col", "Treatment 1 Column", choices = NULL),
                     selectInput("treat2_col", "Treatment 2 Column", choices = NULL)
              ),
              column(4,
                     # Conditional UI for binary outcomes
                     conditionalPanel(
                       condition = "input.effect_measure == 'OR' || input.effect_measure == 'RR'",
                       selectInput("event1_col", "Events Treatment 1", choices = NULL),
                       selectInput("total1_col", "Total Treatment 1", choices = NULL),
                       selectInput("event2_col", "Events Treatment 2", choices = NULL),
                       selectInput("total2_col", "Total Treatment 2", choices = NULL)
                     ),
                     # Conditional UI for continuous outcomes
                     conditionalPanel(
                       condition = "input.effect_measure == 'MD' || input.effect_measure == 'SMD'",
                       selectInput("mean1_col", "Mean Treatment 1", choices = NULL),
                       selectInput("sd1_col", "SD Treatment 1", choices = NULL),
                       selectInput("n1_col", "N Treatment 1", choices = NULL),
                       selectInput("mean2_col", "Mean Treatment 2", choices = NULL),
                       selectInput("sd2_col", "SD Treatment 2", choices = NULL),
                       selectInput("n2_col", "N Treatment 2", choices = NULL)
                     ),
                     # Conditional UI for hazard ratios
                     conditionalPanel(
                       condition = "input.effect_measure == 'HR'",
                       selectInput("hr_col", "Hazard Ratio Column", choices = NULL),
                       selectInput("lowerci_col", "Lower CI Column", choices = NULL),
                       selectInput("upperci_col", "Upper CI Column", choices = NULL)
                     )
              ),
              column(4,
                     selectInput("metareg_vars", "Covariates for Meta-Regression",
                                 choices = NULL, multiple = TRUE),
                     selectInput("metareg_type", "Covariate Type",
                                 choices = c("Continuous" = "continuous",
                                             "Categorical" = "categorical"),
                                 selected = "continuous")
              )
            ),
            div(
              style = "text-align: center;",
              actionButton("run_nma", "Run Network Meta-Analysis", 
                           class = "btn-lg btn-primary", icon = icon("play"))
            )
          )
        ),
        # Add advanced settings section
        fluidRow(
          box(
            width = 12,
            title = "Advanced Settings",
            status = "primary",
            collapsible = TRUE,
            collapsed = TRUE,
            fluidRow(
              column(6,
                     numericInput("level", "Confidence Level", 
                                  value = 0.95, min = 0.5, max = 0.99, step = 0.01),
                     selectInput("method", "Statistical Method for Meta-Analysis",
                                 choices = c("Inverse Variance" = "inverse", 
                                             "Mantel-Haenszel" = "MH",
                                             "Peto" = "Peto"),
                                 selected = "inverse")
              ),
              column(6,
                     selectInput("tau2_estimator", "Tau² Estimator (for Random Effects)",
                                 choices = c("DerSimonian-Laird" = "DL",
                                             "Paule-Mandel" = "PM",
                                             "Restricted Maximum Likelihood" = "REML"),
                                 selected = "DL"),
                     numericInput("n_simulation", "Number of Simulations for Rankings",
                                  value = 1000, min = 100, max = 10000, step = 100)
              )
            )
          )
        )
      ),
      
      # Results Tab
      tabItem(
        tabName = "results",
        fluidRow(
          box(
            width = 12,
            title = "3. Network Meta-Analysis Results",
            status = "primary"
          )
        ),
        tabBox(
          id = "nma_results_tabs",
          width = 12,
          tabPanel("Summary", 
                   fluidRow(
                     column(12, DTOutput("nma_summary"))
                   )),
          tabPanel("League Table", 
                   fluidRow(
                     column(12, DTOutput("league_table"))
                   )),
          tabPanel("Network Plot", 
                   fluidRow(
                     column(8, plotOutput("network_plot", height = "500px")),
                     column(4, 
                            selectInput("network_style", "Plot Style",
                                        choices = c("Default" = "default", 
                                                    "3D" = "3d",
                                                    "Circle" = "circle"),
                                        selected = "default"),
                            checkboxInput("show_weights", "Show line weights based on precision", TRUE),
                            actionButton("redraw_network", "Redraw Network", class = "btn-info"),
                            # Add network plot export
                            downloadButton("download_network_plot", "Download Plot", class = "btn-success download-button")
                     )
                   )),
          tabPanel("Forest Plot", 
                   fluidRow(
                     column(8, plotOutput("forest_plot", height = "500px")),
                     column(4, 
                            selectInput("forest_reference", "Reference Treatment", choices = NULL),
                            actionButton("redraw_forest", "Redraw Forest Plot", class = "btn-info"),
                            downloadButton("download_forest_plot", "Download Plot", class = "btn-success download-button")
                     )
                   )),
          tabPanel("Ranking", 
                   fluidRow(
                     column(12, DTOutput("treatment_ranking")),
                     column(12, plotOutput("ranking_plot", height = "400px")),
                     column(12, 
                            div(
                              class = "download-button text-center",
                              downloadButton("download_rankings", "Download Rankings", class = "btn-success")
                            )
                     )
                   )),
          tabPanel("Heterogeneity", 
                   fluidRow(
                     column(12, verbatimTextOutput("heterogeneity"))
                   )),
          tabPanel("Inconsistency", 
                   fluidRow(
                     column(12, verbatimTextOutput("inconsistency"))
                   )),
          # Add contribution matrix tab
          tabPanel("Contribution Matrix", 
                   fluidRow(
                     column(12, 
                            h4("Direct Evidence Contribution"),
                            p("This matrix shows the percentage contribution of each direct comparison to each network estimate."),
                            DTOutput("contribution_matrix"))
                   ))
        )
      ),
      
      # Meta-regression Summary Tab
      tabItem(
        tabName = "metareg_summary",
        fluidRow(
          box(
            width = 12,
            title = span("4. Meta-regression", span(class = "beta-label", "beta")),
            status = "primary"
          )
        ),
        fluidRow(
          box(
            width = 12,
            title = "4a. Summary",
            status = "primary",
            tabsetPanel(
              id = "summary_tabs",
              tabPanel("Characteristic Plot", 
                       fluidRow(
                         column(8, plotOutput("summary_char_plot", height = "500px")),
                         column(4, 
                                radioButtons("summary_type", "Plot Type",
                                             choices = c("Covariate" = "covariate", 
                                                         "Baseline risk" = "baseline"),
                                             selected = "covariate"),
                                uiOutput("summary_controls")
                         )
                       )),
              # Add study characteristics summary tab
              tabPanel("Study Characteristics Table", 
                       fluidRow(
                         column(12, 
                                h4("Summary of Study Characteristics"),
                                DTOutput("study_characteristics_table"))
                       ))
            )
          )
        )
      ),
      
      # Baseline Risk Analysis Tab
      tabItem(
        tabName = "baseline_risk",
        fluidRow(
          box(
            width = 12,
            title = span("4. Meta-regression", span(class = "beta-label", "beta")),
            status = "primary"
          )
        ),
        fluidRow(
          box(
            width = 12,
            title = "4b. Baseline Risk Analysis",
            status = "primary",
            tabsetPanel(
              id = "baseline_tabs",
              tabPanel("Regression plot", 
                       fluidRow(
                         column(8, plotOutput("baseline_regression_plot", height = "500px")),
                         column(4, 
                                selectInput("baseline_reference", "Reference Treatment", choices = NULL),
                                actionButton("run_baseline", "Run Baseline Risk Analysis", 
                                             class = "btn-info", icon = icon("calculator")),
                                # Add a download button for baseline risk plot
                                downloadButton("download_baseline_plot", "Download Plot", 
                                               class = "btn-success download-button")
                         )
                       )),
              tabPanel("Forest plot", 
                       fluidRow(
                         column(8, plotOutput("baseline_forest_plot", height = "500px")),
                         column(4, 
                                numericInput("baseline_value", "Baseline Risk Value", 
                                             value = 0.3, min = 0, max = 1, step = 0.05),
                                actionButton("redraw_baseline_forest", "Redraw Forest Plot", 
                                             class = "btn-info"),
                                # Add a download button for baseline forest plot
                                downloadButton("download_baseline_forest", "Download Plot", 
                                               class = "btn-success download-button")
                         )
                       )),
              tabPanel("Comparison of all treatment pairs", 
                       fluidRow(
                         column(12, DTOutput("baseline_comparison_table")),
                         column(12, 
                                div(
                                  class = "download-button text-center",
                                  downloadButton("download_baseline_comparison", "Download Table", 
                                                 class = "btn-success")
                                )
                         )
                       )),
              tabPanel("Ranking", 
                       fluidRow(
                         column(8, plotOutput("baseline_ranking_plot", height = "500px")),
                         column(4, 
                                sliderInput("baseline_ranking_value", "Baseline Risk", 
                                            min = 0, max = 1, value = 0.3, step = 0.05),
                                downloadButton("download_baseline_ranking", "Download Rankings", 
                                               class = "btn-success download-button")
                         )
                       )),
              tabPanel("Result details", 
                       fluidRow(
                         column(12, verbatimTextOutput("baseline_details"))
                       ))
            )
          )
        )
      ),
      
      # Covariate Analysis Tab - Enhanced
      tabItem(
        tabName = "covariate_analysis",
        fluidRow(
          box(
            width = 12,
            title = span("4. Meta-regression", span(class = "beta-label", "beta")),
            status = "primary"
          )
        ),
        fluidRow(
          box(
            width = 12,
            title = "4c. Covariate Analysis",
            status = "primary",
            fluidRow(
              column(4,
                     selectInput("covariate_var", "Select Covariate", choices = NULL),
                     selectInput("covariate_type", "Covariate Type", 
                                 choices = c("Continuous" = "continuous", 
                                             "Categorical" = "categorical"),
                                 selected = "continuous"),
                     selectInput("reg_method", "Regression Method",
                                 choices = c("Meta-regression" = "metareg", 
                                             "Random effects meta-regression" = "rma"),
                                 selected = "metareg"),
                     # Add interaction term options
                     checkboxInput("include_interaction", "Include treatment-by-covariate interaction", TRUE),
                     conditionalPanel(
                       condition = "input.include_interaction == true",
                       selectInput("interaction_reference", "Reference Treatment for Interaction", choices = NULL)
                     ),
                     actionButton("run_metareg", "Run Meta-Regression", 
                                  class = "btn-primary", icon = icon("calculator"))
              ),
              column(8,
                     conditionalPanel(
                       condition = "input.covariate_type == 'continuous'",
                       # Add range information for the selected covariate
                       uiOutput("covariate_range_info"),
                       sliderInput("covariate_value", "Covariate Value for Analysis", 
                                   min = 0, max = 100, value = 50, step = 1,
                                   width = "100%")
                     ),
                     conditionalPanel(
                       condition = "input.covariate_type == 'categorical'",
                       selectInput("category_value", "Category Value", choices = NULL)
                     )
              )
            ),
            tabsetPanel(
              id = "covariate_tabs",
              tabPanel("Regression plot", 
                       fluidRow(
                         column(12, plotOutput("covariate_regression_plot", height = "500px")),
                         column(12, 
                                div(
                                  class = "download-button text-right",
                                  downloadButton("download_covariate_regression", "Download Plot", 
                                                 class = "btn-success")
                                )
                         )
                       )),
              tabPanel("Forest plot", 
                       fluidRow(
                         column(12, plotOutput("covariate_forest_plot", height = "500px")),
                         column(12, 
                                div(
                                  class = "download-button text-right",
                                  downloadButton("download_covariate_forest", "Download Plot", 
                                                 class = "btn-success")
                                )
                         )
                       )),
              tabPanel("Comparison of all treatment pairs", 
                       fluidRow(
                         column(12, DTOutput("covariate_comparison_table")),
                         column(12, 
                                div(
                                  class = "download-button text-center",
                                  downloadButton("download_covariate_comparison", "Download Table", 
                                                 class = "btn-success")
                                )
                         )
                       )),
              tabPanel("Ranking", 
                       fluidRow(
                         column(8, plotOutput("covariate_ranking_plot", height = "500px")),
                         column(4,
                                conditionalPanel(
                                  condition = "input.covariate_type == 'continuous'",
                                  sliderInput("ranking_covariate_value", "Covariate Value", 
                                              min = 0, max = 100, value = 50, step = 1)
                                ),
                                conditionalPanel(
                                  condition = "input.covariate_type == 'categorical'",
                                  selectInput("ranking_category_value", "Category Value", choices = NULL)
                                ),
                                downloadButton("download_covariate_ranking", "Download Rankings", 
                                               class = "btn-success download-button")
                         )
                       )),
              tabPanel("Bubble Plot", 
                       fluidRow(
                         column(8, plotOutput("metareg_bubble_plot", height = "500px")),
                         column(4, 
                                selectInput("bubble_treatment", "Select Treatment", choices = NULL),
                                checkboxInput("bubble_weighted", "Weight by Precision", TRUE),
                                downloadButton("download_bubble_plot", "Download Plot", 
                                               class = "btn-success download-button")
                         )
                       )),
              tabPanel("Meta-Regression Model", 
                       fluidRow(
                         column(12, verbatimTextOutput("metareg_model_summary"))
                       )),
              tabPanel("Predicted Effects", 
                       fluidRow(
                         column(12, 
                                h4("Predicted Treatment Effects at Different Covariate Values"),
                                DTOutput("predicted_effects_table"),
                                div(
                                  class = "download-button text-center",
                                  downloadButton("download_predicted_effects", "Download Predictions", 
                                                 class = "btn-success")
                                )
                         )
                       )),
              tabPanel("Model Diagnostics", 
                       fluidRow(
                         column(12, 
                                h4("Meta-Regression Diagnostics"),
                                tabsetPanel(
                                  id = "diagnostics_tabs",
                                  tabPanel("Goodness of Fit", verbatimTextOutput("goodness_of_fit")),
                                  tabPanel("Residual Plot", plotOutput("residual_plot", height = "400px")),
                                  tabPanel("Q-Q Plot", plotOutput("qq_plot", height = "400px"))
                                )
                         )
                       ))
            )
          )
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Reactive values
  data <- reactiveVal(NULL)
  nma_result <- reactiveVal(NULL)
  pairwise_data <- reactiveVal(NULL)
  baseline_result <- reactiveVal(NULL)
  metareg_result <- reactiveVal(NULL)
  
  # Sample data storage
  binary_data <- reactiveVal(NULL)
  continuous_data <- reactiveVal(NULL)
  hazard_data <- reactiveVal(NULL)
  
  # Create sample datasets once at app startup
  observe({
    binary_data(create_sample_data("binary"))
    continuous_data(create_sample_data("continuous"))
    hazard_data(create_sample_data("hazard"))
  }, priority = 1000)
  
  # Navigate to next section
  observeEvent(input$proceed_to_analysis, {
    updateTabItems(session, "sidebar", "analysis")
  })
  
  # Download handlers for sample data
  output$download_binary <- downloadHandler(
    filename = function() {
      "sample_binary_data.csv"
    },
    content = function(file) {
      write.csv(binary_data(), file, row.names = FALSE)
    }
  )
  
  output$download_continuous <- downloadHandler(
    filename = function() {
      "sample_continuous_data.csv"
    },
    content = function(file) {
      write.csv(continuous_data(), file, row.names = FALSE)
    }
  )
  
  output$download_hazard <- downloadHandler(
    filename = function() {
      "sample_hazard_data.csv"
    },
    content = function(file) {
      write.csv(hazard_data(), file, row.names = FALSE)
    }
  )
  
  # Download current data
  output$download_current_data <- downloadHandler(
    filename = function() {
      "nma_data.csv"
    },
    content = function(file) {
      req(data())
      write.csv(data(), file, row.names = FALSE)
    }
  )
  
  # Download forest plot
  output$download_forest_plot <- downloadHandler(
    filename = function() {
      "forest_plot.png"
    },
    content = function(file) {
      req(nma_result(), input$forest_reference)
      
      # Save the plot
      png(file, width = 800, height = 600, res = 100)
      
      # Safe sorting - check if effect.random exists and is numeric
      if (!is.null(nma_result()$effect.random) && is.numeric(nma_result()$effect.random)) {
        # Use the original sorting approach
        forest(nma_result(), 
               reference = input$forest_reference,
               sortvar = nma_result()$effect.random,
               smlab = paste("Compared with", input$forest_reference))
      } else {
        # Fall back to default sorting (alphabetical)
        forest(nma_result(), 
               reference = input$forest_reference,
               smlab = paste("Compared with", input$forest_reference))
      }
      
      dev.off()
    }
  )
  
  # Download network plot
  output$download_network_plot <- downloadHandler(
    filename = function() {
      "network_plot.png"
    },
    content = function(file) {
      req(nma_result())
      
      # Save the plot
      png(file, width = 800, height = 600, res = 100)
      
      netgraph(nma_result(), 
               plastic = FALSE,
               thickness = if (input$show_weights) "se.fixed" else "equal",
               points = TRUE, col = "blue",
               start = "circle")
      
      dev.off()
    }
  )
  
  # Download rankings
  output$download_rankings <- downloadHandler(
    filename = function() {
      "treatment_rankings.csv"
    },
    content = function(file) {
      req(nma_result())
      
      ranking <- netrank(nma_result(), small.values = input$values_direction)
      
      # Create a clean data frame for display
      if (!is.null(ranking$P.score)) {
        scores <- ranking$P.score
        ranking_df <- data.frame(
          Treatment = names(scores),
          P_Score = round(scores * 100, 1),
          stringsAsFactors = FALSE
        )
      } else {
        # Fallback if P-scores not available
        ranking_df <- data.frame(
          Treatment = nma_result()$trts,
          Message = "Ranking not available",
          stringsAsFactors = FALSE
        )
      }
      
      write.csv(ranking_df, file, row.names = FALSE)
    }
  )
  
  # Download baseline plot
  output$download_baseline_plot <- downloadHandler(
    filename = function() {
      "baseline_regression_plot.png"
    },
    content = function(file) {
      req(baseline_result())
      
      # Save the plot
      png(file, width = 800, height = 600, res = 100)
      
      # Recreation of the baseline plot (same code as in the renderPlot section)
      br <- baseline_result()
      
      if (!is.null(br$model)) {
        # Create a scatter plot with regression line
        plot(br$data$baseline_risk, br$data$TE, 
             xlab = "Baseline Risk", 
             ylab = paste("Effect (", input$effect_measure, ")", sep = ""),
             main = paste("Meta-regression with Baseline Risk\nReference:", br$reference),
             pch = 19, col = "blue")
        
        # Add regression line
        abline(br$model$beta[1], br$model$beta[2], col = "red", lwd = 2)
        
        # Add horizontal line at y = 0
        abline(h = 0, lty = 2)
        
        # Add confidence bands if available
        if (requireNamespace("metafor", quietly = TRUE)) {
          # Generate predicted values for plotting
          newdata <- data.frame(baseline_risk = seq(min(br$data$baseline_risk), 
                                                    max(br$data$baseline_risk), 
                                                    length.out = 100))
          pred <- predict(br$model, newmods = cbind(newdata$baseline_risk), level = 0.95)
          
          # Add the prediction interval
          lines(newdata$baseline_risk, pred$pred, col = "red", lwd = 2)
          lines(newdata$baseline_risk, pred$ci.lb, col = "red", lty = 2)
          lines(newdata$baseline_risk, pred$ci.ub, col = "red", lty = 2)
        }
      } else {
        # Create placeholder plot
        plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(-2, 2),
             xlab = "Baseline Risk", ylab = "Effect",
             main = "Baseline Risk Meta-Regression\n(Run analysis first)")
        text(0.5, 0, "Run baseline risk analysis", cex = 1.5)
      }
      
      dev.off()
    }
  )
  
  # Download covariate regression plot
  output$download_covariate_regression <- downloadHandler(
    filename = function() {
      paste0("covariate_regression_", gsub(" ", "_", input$covariate_var), ".png")
    },
    content = function(file) {
      req(metareg_result())
      
      # Save the plot
      png(file, width = 800, height = 600, res = 100)
      
      # Recreation of the covariate regression plot
      mr <- metareg_result()
      
      if (!is.null(mr$model)) {
        # Create a scatter plot with regression line
        plot(mr$data$matched_covs, mr$data$TE, 
             xlab = mr$covariate, 
             ylab = paste("Effect (", input$effect_measure, ")", sep = ""),
             main = paste("Meta-regression with", mr$covariate),
             pch = 19, col = "blue")
        
        # Add regression line
        if (mr$method %in% c("metafor", "rma")) {
          abline(mr$model$beta[1], mr$model$beta[2], col = "red", lwd = 2)
        } else if (mr$method == "lm") {
          abline(mr$model, col = "red", lwd = 2)
        }
        
        # Add horizontal line at y = 0
        abline(h = 0, lty = 2)
        
        # Add confidence bands if available and using metafor
        if (mr$method %in% c("metafor", "rma") && requireNamespace("metafor", quietly = TRUE)) {
          # Generate predicted values for plotting
          x_range <- range(mr$data$matched_covs, na.rm = TRUE)
          newdata <- data.frame(matched_covs = seq(x_range[1], x_range[2], length.out = 100))
          pred <- predict(mr$model, newmods = cbind(newdata$matched_covs), level = 0.95)
          
          # Add the prediction interval
          lines(newdata$matched_covs, pred$pred, col = "red", lwd = 2)
          lines(newdata$matched_covs, pred$ci.lb, col = "red", lty = 2)
          lines(newdata$matched_covs, pred$ci.ub, col = "red", lty = 2)
        }
      } else {
        # Create placeholder plot
        plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(-2, 2),
             xlab = input$covariate_var, ylab = "Effect",
             main = paste("Meta-Regression with", input$covariate_var, "\n(Run analysis first)"))
        text(0.5, 0, "Run meta-regression analysis", cex = 1.5)
      }
      
      dev.off()
    }
  )
  
  # Download bubble plot
  output$download_bubble_plot <- downloadHandler(
    filename = function() {
      paste0("bubble_plot_", gsub(" ", "_", input$bubble_treatment), ".png")
    },
    content = function(file) {
      req(metareg_result(), input$bubble_treatment)
      
      # Save the plot using ggplot's ggsave
      ggsave(file, plot = create_metareg_bubble_plot(
        data = metareg_result()$data,
        xvar = "matched_covs",
        yvar = "TE",
        sizevar = "seTE",
        title = paste("Meta-Regression Bubble Plot for", input$bubble_treatment),
        xlab = metareg_result()$covariate,
        ylab = paste("Effect (", input$effect_measure, ")")
      ), width = 8, height = 6, dpi = 100)
    }
  )
  
  # Load sample data from home page
  observeEvent(input$load_binary_sample, {
    data(binary_data())
    
    # Update column selections
    updateSelectInput(session, "study_col", choices = names(binary_data()), selected = "Study")
    updateSelectInput(session, "treat1_col", choices = names(binary_data()), selected = "Treatment1")
    updateSelectInput(session, "treat2_col", choices = names(binary_data()), selected = "Treatment2")
    updateSelectInput(session, "event1_col", choices = names(binary_data()), selected = "Events1")
    updateSelectInput(session, "total1_col", choices = names(binary_data()), selected = "Total1")
    updateSelectInput(session, "event2_col", choices = names(binary_data()), selected = "Events2")
    updateSelectInput(session, "total2_col", choices = names(binary_data()), selected = "Total2")
    
    # Update effect measure
    updateSelectInput(session, "effect_measure", selected = "OR")
    
    # Update meta-regression variables
    metareg_cols <- setdiff(names(binary_data()), c("Study", "Treatment1", "Treatment2", 
                                                    "Events1", "Total1", "Events2", "Total2"))
    updateSelectInput(session, "metareg_vars", choices = metareg_cols, selected = "Year")
    updateSelectInput(session, "covariate_var", choices = metareg_cols, selected = "Year")
    
    # Navigate to data tab
    updateTabItems(session, "sidebar", "data")
  })
  
  # Load continuous sample data
  observeEvent(input$load_continuous_sample, {
    data(continuous_data())
    
    # Update column selections
    updateSelectInput(session, "study_col", choices = names(continuous_data()), selected = "Study")
    updateSelectInput(session, "treat1_col", choices = names(continuous_data()), selected = "Treatment1")
    updateSelectInput(session, "treat2_col", choices = names(continuous_data()), selected = "Treatment2")
    updateSelectInput(session, "mean1_col", choices = names(continuous_data()), selected = "Mean1")
    updateSelectInput(session, "sd1_col", choices = names(continuous_data()), selected = "SD1")
    updateSelectInput(session, "n1_col", choices = names(continuous_data()), selected = "N1")
    updateSelectInput(session, "mean2_col", choices = names(continuous_data()), selected = "Mean2")
    updateSelectInput(session, "sd2_col", choices = names(continuous_data()), selected = "SD2")
    updateSelectInput(session, "n2_col", choices = names(continuous_data()), selected = "N2")
    
    # Update effect measure
    updateSelectInput(session, "effect_measure", selected = "MD")
    
    # Update meta-regression variables
    metareg_cols <- setdiff(names(continuous_data()), c("Study", "Treatment1", "Treatment2", 
                                                        "Mean1", "SD1", "N1", "Mean2", "SD2", "N2"))
    updateSelectInput(session, "metareg_vars", choices = metareg_cols, selected = "Year")
    updateSelectInput(session, "covariate_var", choices = metareg_cols, selected = "Year")
    
    # Navigate to data tab
    updateTabItems(session, "sidebar", "data")
  })
  
  # Load hazard ratio sample data
  observeEvent(input$load_hazard_sample, {
    data(hazard_data())
    
    # Update column selections
    updateSelectInput(session, "study_col", choices = names(hazard_data()), selected = "Study")
    updateSelectInput(session, "treat1_col", choices = names(hazard_data()), selected = "Treatment1")
    updateSelectInput(session, "treat2_col", choices = names(hazard_data()), selected = "Treatment2")
    updateSelectInput(session, "hr_col", choices = names(hazard_data()), selected = "HR")
    updateSelectInput(session, "lowerci_col", choices = names(hazard_data()), selected = "LowerCI")
    updateSelectInput(session, "upperci_col", choices = names(hazard_data()), selected = "UpperCI")
    
    # Update effect measure
    updateSelectInput(session, "effect_measure", selected = "HR")
    
    # Update meta-regression variables
    metareg_cols <- setdiff(names(hazard_data()), c("Study", "Treatment1", "Treatment2", 
                                                    "HR", "LowerCI", "UpperCI", "SampleSize"))
    updateSelectInput(session, "metareg_vars", choices = metareg_cols, selected = "Year")
    updateSelectInput(session, "covariate_var", choices = metareg_cols, selected = "Year")
    
    # Navigate to data tab
    updateTabItems(session, "sidebar", "data")
  })
  
  # Load sample data from data tab (duplicated buttons for convenience)
  observeEvent(input$load_binary_sample_data, {
    data(binary_data())
    
    updateSelectInput(session, "study_col", choices = names(binary_data()), selected = "Study")
    updateSelectInput(session, "treat1_col", choices = names(binary_data()), selected = "Treatment1")
    updateSelectInput(session, "treat2_col", choices = names(binary_data()), selected = "Treatment2")
    updateSelectInput(session, "event1_col", choices = names(binary_data()), selected = "Events1")
    updateSelectInput(session, "total1_col", choices = names(binary_data()), selected = "Total1")
    updateSelectInput(session, "event2_col", choices = names(binary_data()), selected = "Events2")
    updateSelectInput(session, "total2_col", choices = names(binary_data()), selected = "Total2")
    
    # Update effect measure
    updateSelectInput(session, "effect_measure", selected = "OR")
    
    # Update meta-regression variables
    metareg_cols <- setdiff(names(binary_data()), c("Study", "Treatment1", "Treatment2", 
                                                    "Events1", "Total1", "Events2", "Total2"))
    updateSelectInput(session, "metareg_vars", choices = metareg_cols, selected = "Year")
    updateSelectInput(session, "covariate_var", choices = metareg_cols, selected = "Year")
  })
  
  observeEvent(input$load_continuous_sample_data, {
    data(continuous_data())
    
    updateSelectInput(session, "study_col", choices = names(continuous_data()), selected = "Study")
    updateSelectInput(session, "treat1_col", choices = names(continuous_data()), selected = "Treatment1")
    updateSelectInput(session, "treat2_col", choices = names(continuous_data()), selected = "Treatment2")
    updateSelectInput(session, "mean1_col", choices = names(continuous_data()), selected = "Mean1")
    updateSelectInput(session, "sd1_col", choices = names(continuous_data()), selected = "SD1")
    updateSelectInput(session, "n1_col", choices = names(continuous_data()), selected = "N1")
    updateSelectInput(session, "mean2_col", choices = names(continuous_data()), selected = "Mean2")
    updateSelectInput(session, "sd2_col", choices = names(continuous_data()), selected = "SD2")
    updateSelectInput(session, "n2_col", choices = names(continuous_data()), selected = "N2")
    
    # Update effect measure
    updateSelectInput(session, "effect_measure", selected = "MD")
    
    # Update meta-regression variables
    metareg_cols <- setdiff(names(continuous_data()), c("Study", "Treatment1", "Treatment2", 
                                                        "Mean1", "SD1", "N1", "Mean2", "SD2", "N2"))
    updateSelectInput(session, "metareg_vars", choices = metareg_cols, selected = "Year")
    updateSelectInput(session, "covariate_var", choices = metareg_cols, selected = "Year")
  })
  
  observeEvent(input$load_hazard_sample_data, {
    data(hazard_data())
    
    updateSelectInput(session, "study_col", choices = names(hazard_data()), selected = "Study")
    updateSelectInput(session, "treat1_col", choices = names(hazard_data()), selected = "Treatment1")
    updateSelectInput(session, "treat2_col", choices = names(hazard_data()), selected = "Treatment2")
    updateSelectInput(session, "hr_col", choices = names(hazard_data()), selected = "HR")
    updateSelectInput(session, "lowerci_col", choices = names(hazard_data()), selected = "LowerCI")
    updateSelectInput(session, "upperci_col", choices = names(hazard_data()), selected = "UpperCI")
    
    # Update effect measure
    updateSelectInput(session, "effect_measure", selected = "HR")
    
    # Update meta-regression variables
    metareg_cols <- setdiff(names(hazard_data()), c("Study", "Treatment1", "Treatment2", 
                                                    "HR", "LowerCI", "UpperCI", "SampleSize"))
    updateSelectInput(session, "metareg_vars", choices = metareg_cols, selected = "Year")
    updateSelectInput(session, "covariate_var", choices = metareg_cols, selected = "Year")
  })
  
  # Read data file
  observeEvent(input$file, {
    req(input$file)
    
    tryCatch({
      ext <- tools::file_ext(input$file$name)
      
      if (ext == "csv") {
        df <- read.csv(input$file$datapath, header = input$header, sep = input$sep)
      } else if (ext %in% c("xls", "xlsx")) {
        df <- readxl::read_excel(input$file$datapath)
      } else {
        showNotification("Unsupported file type. Please upload CSV or Excel files.", type = "error")
        return(NULL)
      }
      
      data(df)
      
      # Update column choices for all inputs
      updateSelectInput(session, "study_col", choices = names(df))
      updateSelectInput(session, "treat1_col", choices = names(df))
      updateSelectInput(session, "treat2_col", choices = names(df))
      updateSelectInput(session, "event1_col", choices = names(df))
      updateSelectInput(session, "total1_col", choices = names(df))
      updateSelectInput(session, "event2_col", choices = names(df))
      updateSelectInput(session, "total2_col", choices = names(df))
      updateSelectInput(session, "mean1_col", choices = names(df))
      updateSelectInput(session, "sd1_col", choices = names(df))
      updateSelectInput(session, "n1_col", choices = names(df))
      updateSelectInput(session, "mean2_col", choices = names(df))
      updateSelectInput(session, "sd2_col", choices = names(df))
      updateSelectInput(session, "n2_col", choices = names(df))
      updateSelectInput(session, "hr_col", choices = names(df))
      updateSelectInput(session, "lowerci_col", choices = names(df))
      updateSelectInput(session, "upperci_col", choices = names(df))
      
      # Update meta-regression variables
      metareg_cols <- setdiff(names(df), c(input$study_col, input$treat1_col, input$treat2_col))
      updateSelectInput(session, "metareg_vars", choices = metareg_cols)
      updateSelectInput(session, "covariate_var", choices = metareg_cols)
      
      showNotification("File loaded successfully", type = "message")
    }, error = function(e) {
      showNotification(paste("Error loading file:", e$message), type = "error")
    })
  })
  
  # Data preview
  output$data_preview <- renderDT({
    req(data())
    datatable(data(), options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # Data validation output
  output$data_validation_output <- renderPrint({
    req(data())
    
    cat("Data Validation Results:\n\n")
    
    # Check for missing values
    missing_count <- sum(is.na(data()))
    cat("Missing values:", missing_count, "\n")
    
    # Check for duplicate rows
    dup_count <- sum(duplicated(data()))
    cat("Duplicate rows:", dup_count, "\n")
    
    # Check data types for required columns based on effect measure
    if (input$effect_measure %in% c("OR", "RR")) {
      # Check numeric columns for binary outcomes
      if (!is.null(input$event1_col) && !is.null(data()[[input$event1_col]])) {
        cat("Events1 column (", input$event1_col, "): ", 
            ifelse(is.numeric(data()[[input$event1_col]]), "Numeric ✓", "NOT Numeric ✗"), "\n")
      }
      if (!is.null(input$total1_col) && !is.null(data()[[input$total1_col]])) {
        cat("Total1 column (", input$total1_col, "): ", 
            ifelse(is.numeric(data()[[input$total1_col]]), "Numeric ✓", "NOT Numeric ✗"), "\n")
      }
      if (!is.null(input$event2_col) && !is.null(data()[[input$event2_col]])) {
        cat("Events2 column (", input$event2_col, "): ", 
            ifelse(is.numeric(data()[[input$event2_col]]), "Numeric ✓", "NOT Numeric ✗"), "\n")
      }
      if (!is.null(input$total2_col) && !is.null(data()[[input$total2_col]])) {
        cat("Total2 column (", input$total2_col, "): ", 
            ifelse(is.numeric(data()[[input$total2_col]]), "Numeric ✓", "NOT Numeric ✗"), "\n")
      }
    } else if (input$effect_measure %in% c("MD", "SMD")) {
      # Check numeric columns for continuous outcomes
      if (!is.null(input$mean1_col) && !is.null(data()[[input$mean1_col]])) {
        cat("Mean1 column (", input$mean1_col, "): ", 
            ifelse(is.numeric(data()[[input$mean1_col]]), "Numeric ✓", "NOT Numeric ✗"), "\n")
      }
      if (!is.null(input$sd1_col) && !is.null(data()[[input$sd1_col]])) {
        cat("SD1 column (", input$sd1_col, "): ", 
            ifelse(is.numeric(data()[[input$sd1_col]]), "Numeric ✓", "NOT Numeric ✗"), "\n")
      }
      if (!is.null(input$n1_col) && !is.null(data()[[input$n1_col]])) {
        cat("N1 column (", input$n1_col, "): ", 
            ifelse(is.numeric(data()[[input$n1_col]]), "Numeric ✓", "NOT Numeric ✗"), "\n")
      }
      if (!is.null(input$mean2_col) && !is.null(data()[[input$mean2_col]])) {
        cat("Mean2 column (", input$mean2_col, "): ", 
            ifelse(is.numeric(data()[[input$mean2_col]]), "Numeric ✓", "NOT Numeric ✗"), "\n")
      }
      if (!is.null(input$sd2_col) && !is.null(data()[[input$sd2_col]])) {
        cat("SD2 column (", input$sd2_col, "): ", 
            ifelse(is.numeric(data()[[input$sd2_col]]), "Numeric ✓", "NOT Numeric ✗"), "\n")
      }
      if (!is.null(input$n2_col) && !is.null(data()[[input$n2_col]])) {
        cat("N2 column (", input$n2_col, "): ", 
            ifelse(is.numeric(data()[[input$n2_col]]), "Numeric ✓", "NOT Numeric ✗"), "\n")
      }
    } else if (input$effect_measure == "HR") {
      # Check numeric columns for hazard ratio outcomes
      if (!is.null(input$hr_col) && !is.null(data()[[input$hr_col]])) {
        cat("HR column (", input$hr_col, "): ", 
            ifelse(is.numeric(data()[[input$hr_col]]), "Numeric ✓", "NOT Numeric ✗"), "\n")
      }
      if (!is.null(input$lowerci_col) && !is.null(data()[[input$lowerci_col]])) {
        cat("Lower CI column (", input$lowerci_col, "): ", 
            ifelse(is.numeric(data()[[input$lowerci_col]]), "Numeric ✓", "NOT Numeric ✗"), "\n")
      }
      if (!is.null(input$upperci_col) && !is.null(data()[[input$upperci_col]])) {
        cat("Upper CI column (", input$upperci_col, "): ", 
            ifelse(is.numeric(data()[[input$upperci_col]]), "Numeric ✓", "NOT Numeric ✗"), "\n")
      }
    }
    
    # Check treatment connectivity
    if (!is.null(input$treat1_col) && !is.null(input$treat2_col) && 
        !is.null(data()[[input$treat1_col]]) && !is.null(data()[[input$treat2_col]])) {
      
      unique_treatments <- unique(c(data()[[input$treat1_col]], data()[[input$treat2_col]]))
      cat("\nUnique treatments found:", length(unique_treatments), "\n")
      cat("Treatments:", paste(sort(unique_treatments), collapse = ", "), "\n\n")
      
      # Calculate network connectivity (simple check)
      treatments_in_t1 <- unique(data()[[input$treat1_col]])
      treatments_in_t2 <- unique(data()[[input$treat2_col]])
      
      treatments_in_both <- intersect(treatments_in_t1, treatments_in_t2)
      cat("Treatments that appear in both T1 and T2:", length(treatments_in_both), "\n")
      cat("This is a rough indicator of network connectivity.\n")
    }
  })
  
  # Study selection UI
  output$study_selection_ui <- renderUI({
    req(data(), input$study_col)
    
    # Validate that the study column exists in the data
    if (!input$study_col %in% names(data())) {
      return(div(
        class = "alert alert-warning",
        "Please select a valid study column from your data."
      ))
    }
    
    studies <- unique(data()[[input$study_col]])
    checkboxGroupInput("included_studies", "Select studies to include:", 
                       choices = studies, selected = studies)
  })
  
  # Covariate range information
  output$covariate_range_info <- renderUI({
    req(data(), input$covariate_var)
    
    # Get the selected covariate
    cov_data <- data()[[input$covariate_var]]
    
    # Only proceed if the covariate is numeric
    if (is.numeric(cov_data)) {
      min_val <- min(cov_data, na.rm = TRUE)
      max_val <- max(cov_data, na.rm = TRUE)
      mean_val <- mean(cov_data, na.rm = TRUE)
      
      # Update the slider input with the actual range
      updateSliderInput(session, "covariate_value", 
                        min = min_val, max = max_val, value = mean_val,
                        step = (max_val - min_val) / 100)
      
      # Also update the ranking slider
      updateSliderInput(session, "ranking_covariate_value", 
                        min = min_val, max = max_val, value = mean_val,
                        step = (max_val - min_val) / 100)
      
      # Return the range information
      div(
        class = "info-box",
        p(strong("Covariate Range:"), 
          paste("Min =", round(min_val, 2), ", Max =", round(max_val, 2), 
                ", Mean =", round(mean_val, 2)))
      )
    } else {
      div(
        class = "alert alert-warning",
        "Selected covariate is not numeric. Please select a numeric covariate for continuous analysis."
      )
    }
  })
  
  # Summary controls UI
  output$summary_controls <- renderUI({
    if (input$summary_type == "covariate") {
      req(input$metareg_vars)
      selectInput("summary_covariate", "Select Covariate", 
                  choices = input$metareg_vars, 
                  selected = if(length(input$metareg_vars) > 0) input$metareg_vars[1] else NULL)
    } else {
      req(input$baseline_reference)
      selectInput("summary_baseline", "Reference Treatment", 
                  choices = if(!is.null(nma_result())) nma_result()$trts else NULL, 
                  selected = input$baseline_reference)
    }
  })
  
  # Run Network Meta-Analysis
  observeEvent(input$run_nma, {
    req(data())
    
    # Validate required inputs
    if (is.null(input$study_col) || input$study_col == "" || 
        is.null(input$treat1_col) || input$treat1_col == "" || 
        is.null(input$treat2_col) || input$treat2_col == "") {
      showNotification("Please select valid study and treatment columns", type = "error")
      return()
    }
    
    # Create a progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Running network meta-analysis...", value = 0.1)
    on.exit(progress$close())
    
    # Get data and filter by selected studies if applicable
    df <- data()
    if (!is.null(input$included_studies)) {
      df <- df[df[[input$study_col]] %in% input$included_studies, ]
    }
    
    # Create pairwise dataset
    pw_data <- NULL
    
    progress$set(message = "Preparing data...", value = 0.3)
    
    tryCatch({
      if (input$effect_measure %in% c("OR", "RR")) {
        # Binary outcome
        req(input$treat1_col, input$treat2_col, input$event1_col, input$total1_col, 
            input$event2_col, input$total2_col, input$study_col)
        
        # Prepare data for netmeta
        pw_data <- df %>%
          select(
            studlab = !!sym(input$study_col),
            treat1 = !!sym(input$treat1_col),
            treat2 = !!sym(input$treat2_col),
            event1 = !!sym(input$event1_col),
            n1 = !!sym(input$total1_col),
            event2 = !!sym(input$event2_col),
            n2 = !!sym(input$total2_col)
          )
        
        # Add selected covariates
        if (!is.null(input$metareg_vars) && length(input$metareg_vars) > 0) {
          for (cov in input$metareg_vars) {
            if (cov %in% names(df)) {
              pw_data[[cov]] <- df[[cov]]
            }
          }
        }
        
        # Calculate effect sizes with metabin
        meta_bin <- metabin(event1, n1, event2, n2, 
                            studlab = studlab, 
                            data = pw_data,
                            sm = input$effect_measure,
                            method = input$method,
                            level = input$level)
        
        # Extract TE and seTE for netmeta
        pw_data$TE <- meta_bin$TE
        pw_data$seTE <- meta_bin$seTE
        
      } else if (input$effect_measure %in% c("MD", "SMD")) {
        # Continuous outcome
        req(input$treat1_col, input$treat2_col, input$mean1_col, input$sd1_col, 
            input$n1_col, input$mean2_col, input$sd2_col, input$n2_col, input$study_col)
        
        # Prepare data for netmeta
        pw_data <- df %>%
          select(
            studlab = !!sym(input$study_col),
            treat1 = !!sym(input$treat1_col),
            treat2 = !!sym(input$treat2_col),
            n1 = !!sym(input$n1_col),
            mean1 = !!sym(input$mean1_col),
            sd1 = !!sym(input$sd1_col),
            n2 = !!sym(input$n2_col),
            mean2 = !!sym(input$mean2_col),
            sd2 = !!sym(input$sd2_col)
          )
        
        # Add selected covariates
        if (!is.null(input$metareg_vars) && length(input$metareg_vars) > 0) {
          for (cov in input$metareg_vars) {
            if (cov %in% names(df)) {
              pw_data[[cov]] <- df[[cov]]
            }
          }
        }
        
        # Calculate effect sizes with metacont
        meta_cont <- metacont(n1, mean1, sd1, n2, mean2, sd2, 
                              studlab = studlab, 
                              data = pw_data,
                              sm = input$effect_measure,
                              method = input$method,
                              level = input$level)
        
        # Extract TE and seTE for netmeta
        pw_data$TE <- meta_cont$TE
        pw_data$seTE <- meta_cont$seTE
        
      } else if (input$effect_measure == "HR") {
        # Hazard ratio
        req(input$treat1_col, input$treat2_col, input$study_col, 
            input$hr_col, input$lowerci_col, input$upperci_col)
        
        # Prepare data for netmeta
        pw_data <- df %>%
          select(
            studlab = !!sym(input$study_col),
            treat1 = !!sym(input$treat1_col),
            treat2 = !!sym(input$treat2_col),
            HR = !!sym(input$hr_col),
            lower = !!sym(input$lowerci_col),
            upper = !!sym(input$upperci_col)
          )
        
        # Add selected covariates
        if (!is.null(input$metareg_vars) && length(input$metareg_vars) > 0) {
          for (cov in input$metareg_vars) {
            if (cov %in% names(df)) {
              pw_data[[cov]] <- df[[cov]]
            }
          }
        }
        
        # Calculate log HR and SE for netmeta
        pw_data$TE <- log(pw_data$HR)
        pw_data$seTE <- (log(pw_data$upper) - log(pw_data$lower)) / (2 * 1.96)
      }
      
      # Store pairwise data
      pairwise_data(pw_data)
      
      progress$set(message = "Running network meta-analysis...", value = 0.6)
      
      # Run network meta-analysis
      nma <- netmeta(TE, seTE, treat1, treat2, studlab, 
                     data = pw_data, 
                     sm = input$effect_measure,
                     random = input$random_effects,
                     small.values = input$values_direction,
                     level = input$level,
                     method.tau = input$tau2_estimator)
      
      nma_result(nma)
      
      # Update reference treatment selections
      updateSelectInput(session, "forest_reference", 
                        choices = nma$trts, selected = nma$trts[1])
      updateSelectInput(session, "baseline_reference", 
                        choices = nma$trts, selected = nma$trts[1])
      updateSelectInput(session, "summary_baseline", 
                        choices = nma$trts, selected = nma$trts[1])
      updateSelectInput(session, "interaction_reference", 
                        choices = nma$trts, selected = nma$trts[1])
      updateSelectInput(session, "bubble_treatment",
                        choices = nma$trts, selected = nma$trts[1])
      
      # Switch to results tab
      updateTabItems(session, "sidebar", "results")
      
      showNotification("Network meta-analysis completed successfully", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error in NMA:", e$message), type = "error")
    })
    
    progress$set(message = "Analysis complete!", value = 1)
  })
  
  # NMA Summary Table
  output$nma_summary <- renderDT({
    req(nma_result())
    
    result_table <- summary(nma_result())$comparison
    datatable(result_table, options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # League Table
  output$league_table <- renderDT({
    req(nma_result())
    
    league_table <- netleague(nma_result(), digits = 2)
    if (input$random_effects) {
      datatable(league_table$random, options = list(scrollX = TRUE))
    } else {
      datatable(league_table$fixed, options = list(scrollX = TRUE))
    }
  })
  
  # Contribution Matrix - Fixed implementation
  output$contribution_matrix <- renderDT({
    req(nma_result())
    
    # Calculate actual contribution matrix using netcontrib
    contrib <- netcontrib(nma_result())
    
    # Extract the contribution matrix
    contrib_matrix <- contrib$contribution.matrix.random
    
    # Convert to percentage
    contrib_matrix <- contrib_matrix * 100
    
    # Convert to data frame for display
    contrib_df <- as.data.frame(contrib_matrix)
    contrib_df$Comparison <- rownames(contrib_df)
    contrib_df <- contrib_df[, c("Comparison", setdiff(names(contrib_df), "Comparison"))]
    
    datatable(contrib_df, 
              caption = "Percentage Contribution of Direct Comparisons to Network Estimates",
              options = list(scrollX = TRUE, pageLength = 10)) %>%
      formatRound(columns = 2:ncol(contrib_df), digits = 1)
  })
  
  # Study characteristics table
  output$study_characteristics_table <- renderDT({
    req(data(), input$study_col)
    
    # Create a summary table of study characteristics
    study_data <- data()
    
    # Get all potential covariates
    if (!is.null(input$metareg_vars) && length(input$metareg_vars) > 0) {
      covariate_cols <- input$metareg_vars
    } else {
      # Fallback to identify likely covariates (exclude treatment and outcome columns)
      outcome_cols <- c(
        input$study_col, input$treat1_col, input$treat2_col,
        input$event1_col, input$total1_col, input$event2_col, input$total2_col,
        input$mean1_col, input$sd1_col, input$n1_col, input$mean2_col, input$sd2_col, input$n2_col,
        input$hr_col, input$lowerci_col, input$upperci_col
      )
      covariate_cols <- setdiff(names(study_data), outcome_cols)
    }
    
    # Create a summary data frame
    summary_df <- data.frame(Study = unique(study_data[[input$study_col]]))
    
    # Add summary for each covariate
    for (cov in covariate_cols) {
      if (cov %in% names(study_data)) {
        # Calculate study-level averages for numeric covariates
        if (is.numeric(study_data[[cov]])) {
          # Group by study and calculate mean
          study_means <- aggregate(study_data[[cov]], 
                                   by = list(Study = study_data[[input$study_col]]), 
                                   FUN = mean, na.rm = TRUE)
          
          # Rename the column
          names(study_means)[2] <- cov
          
          # Merge with summary_df
          summary_df <- merge(summary_df, study_means, by = "Study", all.x = TRUE)
        } else {
          # For categorical covariates, take the most frequent value
          study_modes <- aggregate(as.character(study_data[[cov]]),
                                   by = list(Study = study_data[[input$study_col]]),
                                   FUN = function(x) {
                                     tab <- table(x)
                                     names(tab)[which.max(tab)]
                                   })
          
          # Rename the column
          names(study_modes)[2] <- cov
          
          # Merge with summary_df
          summary_df <- merge(summary_df, study_modes, by = "Study", all.x = TRUE)
        }
      }
    }
    
    datatable(summary_df, options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # Network Plot
  output$network_plot <- renderPlot({
    req(nma_result())
    
    netgraph(nma_result(), 
             plastic = FALSE,
             thickness = if (input$show_weights) "se.fixed" else "equal",
             points = TRUE, col = "blue",
             start = "circle")
  })
  
  # Redraw Network Plot
  observeEvent(input$redraw_network, {
    req(nma_result())
    # This will trigger the renderPlot for network_plot
  })
  
  # Forest Plot - FIXED VERSION
  output$forest_plot <- renderPlot({
    req(nma_result(), input$forest_reference)
    
    # Safe sorting - check if effect.random exists and is numeric
    if (!is.null(nma_result()$effect.random) && is.numeric(nma_result()$effect.random)) {
      # Use the original sorting approach
      forest(nma_result(), 
             reference = input$forest_reference,
             sortvar = nma_result()$effect.random,
             smlab = paste("Compared with", input$forest_reference))
    } else {
      # Fall back to default sorting (alphabetical)
      forest(nma_result(), 
             reference = input$forest_reference,
             smlab = paste("Compared with", input$forest_reference))
    }
  })
  
  # Redraw Forest Plot
  observeEvent(input$redraw_forest, {
    req(nma_result(), input$forest_reference)
    # This will trigger the renderPlot for forest_plot
  })
  
  # Treatment Ranking
  output$treatment_ranking <- renderDT({
    req(nma_result())
    
    # Create ranking table
    ranking <- netrank(nma_result(), small.values = input$values_direction)
    
    # Create a clean data frame for display
    if (!is.null(ranking$P.score)) {
      scores <- ranking$P.score
      ranking_df <- data.frame(
        Treatment = names(scores),
        P_Score = round(scores * 100, 1),
        stringsAsFactors = FALSE
      )
    } else {
      # Fallback if P-scores not available
      ranking_df <- data.frame(
        Treatment = nma_result()$trts,
        Message = "Ranking not available",
        stringsAsFactors = FALSE
      )
    }
    
    datatable(ranking_df, options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # Ranking Plot
  output$ranking_plot <- renderPlot({
    req(nma_result())
    
    ranking <- netrank(nma_result(), small.values = input$values_direction)
    
    # Create P-score barplot
    if (!is.null(ranking$P.score)) {
      p_scores <- ranking$P.score * 100
      barplot(p_scores, 
              main = paste("P-Score Rankings (", 
                           ifelse(input$values_direction == "good", 
                                  "smaller values are better", 
                                  "larger values are better"), 
                           ")", sep = ""),
              xlab = "Treatment", ylab = "P-Score (%)",
              col = "skyblue",
              ylim = c(0, 100),
              names.arg = names(p_scores))
    } else {
      # Placeholder if P-scores are not available
      plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "",
           main = "Ranking information not available")
      text(0, 0, "Ranking plot cannot be generated with current data")
    }
  })
  
  # Heterogeneity
  output$heterogeneity <- renderPrint({
    req(nma_result())
    
    # Calculate and print heterogeneity statistics
    cat("Quantifying heterogeneity and inconsistency:\n\n")
    cat("I² = ", round(nma_result()$I2 * 100), "% (between studies)\n", sep = "")
    cat("tau² = ", round(nma_result()$tau^2, 5), "\n", sep = "")
    
    # Add interpretation of I² values
    cat("\nInterpretation guideline for I²:\n")
    cat("0-25%: Might not be important\n")
    cat("25-50%: May represent moderate heterogeneity\n")
    cat("50-75%: May represent substantial heterogeneity\n")
    cat(">75%: Considerable heterogeneity\n\n")
    
    # Check for closed loops in the network
    nma <- nma_result()
    
    # Calculate the number of treatments
    n_treats <- length(nma$trts)
    n_comparisons <- nrow(nma$data)
    
    # A simple check: for a fully connected network, we need at least n_treats - 1 comparisons
    # For closed loops, we need at least n_treats comparisons
    has_loops <- n_comparisons >= n_treats
    
    if (has_loops) {
      cat("\nTests for inconsistency:\n")
      tryCatch({
        ns <- netsplit(nma_result())
        print(ns)
      }, error = function(e) {
        cat("Error in inconsistency assessment:", e$message, "\n")
      })
    } else {
      cat("\nNo closed loops in the network, inconsistency cannot be assessed.")
    }
  })
  
  # Inconsistency
  output$inconsistency <- renderPrint({
    req(nma_result())
    
    nma <- nma_result()
    
    # Calculate the number of treatments and comparisons
    n_treats <- length(nma$trts)
    n_comparisons <- nrow(nma$data)
    
    # Check if we have enough data for inconsistency assessment
    has_loops <- n_comparisons >= n_treats
    
    if (has_loops) {
      tryCatch({
        decomp <- decomp.design(nma_result())
        print(decomp)
        cat("\nQ statistics to assess inconsistency:\n")
        cat("Within-design heterogeneity Q = ", round(decomp$Q.heterogeneity, 2), 
            ", df = ", decomp$df.Q.heterogeneity, 
            ", p = ", round(decomp$pval.Q.heterogeneity, 4), "\n", sep = "")
        cat("Between-design inconsistency Q = ", round(decomp$Q.inconsistency, 2), 
            ", df = ", decomp$df.Q.inconsistency, 
            ", p = ", round(decomp$pval.Q.inconsistency, 4), "\n", sep = "")
      }, error = function(e) {
        cat("Error in inconsistency assessment:", e$message, "\n")
        cat("This may occur when the network structure doesn't allow for inconsistency assessment.\n")
      })
    } else {
      cat("No closed loops in the network, inconsistency cannot be assessed.")
    }
  })
  
  # Summary Plot
  output$summary_char_plot <- renderPlot({
    req(nma_result())
    
    if (input$summary_type == "covariate") {
      req(input$summary_covariate, data())
      
      # Get covariate data
      covariate_data <- data()[[input$summary_covariate]]
      
      # Get pairwise data if available
      if (!is.null(pairwise_data())) {
        pw <- pairwise_data()
        
        # Match covariate data to pairwise data
        pw$covariate <- data()[[input$summary_covariate]][match(pw$studlab, data()[[input$study_col]])]
        
        # Create scatter plot
        plot(pw$covariate, pw$TE,
             xlab = input$summary_covariate, 
             ylab = paste("Effect (", input$effect_measure, ")", sep = ""),
             main = paste("Effect by", input$summary_covariate),
             pch = 19, col = "blue")
        
        # Add a simple regression line
        fit <- lm(TE ~ covariate, data = pw)
        abline(fit, col = "red", lwd = 2)
        abline(h = 0, lty = 2)
        
        # Add R-squared value
        r2 <- summary(fit)$r.squared
        legend("topright", 
               legend = paste("R² =", round(r2, 3)),
               bty = "n")
      } else {
        # Create placeholder if no pairwise data
        plot(1, 1, type = "n", xlab = input$summary_covariate, ylab = "Effect",
             main = paste("Effect by", input$summary_covariate),
             xlim = c(min(covariate_data, na.rm = TRUE), max(covariate_data, na.rm = TRUE)),
             ylim = c(-2, 2))
        abline(h = 0, lty = 2)
        text(mean(range(covariate_data, na.rm = TRUE)), 0, 
             "Run network meta-analysis first", cex = 1.2)
      }
      
    } else {
      req(input$summary_baseline)
      
      # Create baseline risk plot if baseline analysis has been run
      if (!is.null(baseline_result())) {
        br <- baseline_result()
        
        plot(br$data$baseline_risk, br$data$TE, 
             xlab = "Baseline Risk", 
             ylab = paste("Effect (", input$effect_measure, ")", sep = ""),
             main = paste("Effect by Baseline Risk (Reference:", input$summary_baseline, ")"),
             pch = 19, col = "blue")
        
        # Add regression line
        abline(br$model$beta[1], br$model$beta[2], col = "red", lwd = 2)
        abline(h = 0, lty = 2)
      } else {
        # Create placeholder
        plot(1, 1, type = "n", xlab = "Baseline Risk", ylab = "Effect",
             main = paste("Effect by Baseline Risk (Reference:", input$summary_baseline, ")"),
             xlim = c(0, 1), ylim = c(-2, 2))
        abline(h = 0, lty = 2)
        text(0.5, 0, "Run baseline risk analysis first", cex = 1.2)
      }
    }
  })
  
  # Baseline Risk Analysis
  observeEvent(input$run_baseline, {
    req(nma_result(), pairwise_data(), input$baseline_reference)
    
    # Create progress indicator
    progress <- shiny::Progress$new()
    progress$set(message = "Running baseline risk analysis...", value = 0.1)
    on.exit(progress$close())
    
    # Extract the data needed for baseline risk analysis
    tryCatch({
      # Get pairwise data
      pw <- pairwise_data()
      
      # Calculate baseline risk (for binary outcomes)
      if (input$effect_measure %in% c("OR", "RR")) {
        # For binary outcomes, baseline risk is event rate in control group
        baseline_data <- pw %>%
          mutate(
            baseline_risk = ifelse(treat1 == input$baseline_reference, 
                                   event1 / n1,
                                   ifelse(treat2 == input$baseline_reference, 
                                          event2 / n2, 
                                          NA))
          ) %>%
          filter(!is.na(baseline_risk))
      } else {
        # For continuous/HR outcomes, calculate a pseudo baseline risk
        # based on the control group mean or median effect
        baseline_data <- pw %>%
          filter(treat1 == input$baseline_reference | treat2 == input$baseline_reference)
        
        if (nrow(baseline_data) > 0) {
          # Use a normalized version of the control group response
          if (input$effect_measure %in% c("MD", "SMD")) {
            baseline_data$baseline_risk <- ifelse(
              baseline_data$treat1 == input$baseline_reference,
              baseline_data$mean1 / max(abs(c(baseline_data$mean1, baseline_data$mean2)), na.rm = TRUE),
              baseline_data$mean2 / max(abs(c(baseline_data$mean1, baseline_data$mean2)), na.rm = TRUE)
            )
          } else {
            # For HR, use a simulated baseline risk
            baseline_data$baseline_risk <- runif(nrow(baseline_data), 0.1, 0.4)
          }
        }
      }
      
      progress$set(message = "Fitting meta-regression model...", value = 0.5)
      
      # Create a basic meta-regression model using metafor
      if (requireNamespace("metafor", quietly = TRUE)) {
        # Handle possible empty data
        if (nrow(baseline_data) == 0) {
          showNotification("No studies with selected reference treatment. Please select another reference.", type = "error")
          return(NULL)
        }
        
        # Simple meta-regression model
        mr_model <- rma(yi = TE, sei = seTE, mods = ~ baseline_risk, 
                        data = baseline_data, method = "REML")
        
        # Store results
        baseline_result(list(
          reference = input$baseline_reference,
          data = baseline_data,
          model = mr_model,
          covariate = "baseline_risk",
          type = "continuous"
        ))
        
        # Show success notification
        showNotification("Baseline risk analysis completed", type = "message")
      } else {
        showNotification("Package 'metafor' is required for meta-regression", type = "error")
      }
      
    }, error = function(e) {
      showNotification(paste("Error in baseline risk analysis:", e$message), type = "error")
    })
    
    progress$set(message = "Baseline risk analysis complete!", value = 1)
  })
  
  # Baseline regression plot
  output$baseline_regression_plot <- renderPlot({
    req(baseline_result())
    
    # Create the plot
    br <- baseline_result()
    
    if (!is.null(br$model)) {
      # Create a scatter plot with regression line
      plot(br$data$baseline_risk, br$data$TE, 
           xlab = "Baseline Risk", 
           ylab = paste("Effect (", input$effect_measure, ")", sep = ""),
           main = paste("Meta-regression with Baseline Risk\nReference:", br$reference),
           pch = 19, col = "blue")
      
      # Add regression line
      abline(br$model$beta[1], br$model$beta[2], col = "red", lwd = 2)
      
      # Add horizontal line at y = 0
      abline(h = 0, lty = 2)
      
      # Add confidence bands
      newdata <- data.frame(baseline_risk = seq(min(br$data$baseline_risk), 
                                                max(br$data$baseline_risk), 
                                                length.out = 100))
      pred <- predict(br$model, newmods = cbind(newdata$baseline_risk), level = 0.95)
      
      # Add the prediction interval
      lines(newdata$baseline_risk, pred$pred, col = "red", lwd = 2)
      lines(newdata$baseline_risk, pred$ci.lb, col = "red", lty = 2)
      lines(newdata$baseline_risk, pred$ci.ub, col = "red", lty = 2)
      
      # Add legend
      legend("topright", 
             legend = c("Observed", "Fitted", "95% CI"),
             col = c("blue", "red", "red"),
             lty = c(NA, 1, 2),
             pch = c(19, NA, NA),
             bty = "n")
    } else {
      # Create placeholder plot
      plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(-2, 2),
           xlab = "Baseline Risk", ylab = "Effect",
           main = "Baseline Risk Meta-Regression\n(Run analysis first)")
      text(0.5, 0, "Run baseline risk analysis", cex = 1.5)
    }
  })
  
  # Baseline forest plot
  output$baseline_forest_plot <- renderPlot({
    req(baseline_result(), nma_result(), input$baseline_value)
    
    # Create adjusted forest plot
    create_metareg_forest_plot(
      nma = nma_result(),
      metareg_model = baseline_result()$model,
      covariate_value = input$baseline_value,
      reference_treatment = baseline_result()$reference,
      sm = input$effect_measure,
      cov_name = "Baseline Risk"
    )
  })
  
  # Baseline comparison table
  output$baseline_comparison_table <- renderDT({
    req(baseline_result(), nma_result())
    
    # Get all treatments
    treatments <- nma_result()$trts
    n_treats <- length(treatments)
    
    # Create comparison matrix
    comparison_matrix <- matrix(NA, nrow = n_treats, ncol = n_treats,
                                dimnames = list(treatments, treatments))
    
    # Get baseline risk value
    baseline_value <- input$baseline_value
    
    # Calculate adjusted effects
    br <- baseline_result()
    adjustment <- 0
    
    if (!is.null(br$model)) {
      # Calculate adjustment based on baseline risk
      mean_baseline <- mean(br$data$baseline_risk, na.rm = TRUE)
      adjustment <- br$model$beta[2] * (baseline_value - mean_baseline)
    }
    
    # Fill comparison matrix with adjusted NMA estimates
    for (i in 1:n_treats) {
      for (j in 1:n_treats) {
        if (i != j) {
          # Find the comparison in NMA results
          idx <- which((nma_result()$treat1 == treatments[i] & nma_result()$treat2 == treatments[j]) |
                       (nma_result()$treat1 == treatments[j] & nma_result()$treat2 == treatments[i]))
          
          if (length(idx) > 0) {
            if (nma_result()$treat1[idx] == treatments[i]) {
              comparison_matrix[i, j] <- round(nma_result()$TE.nma.random[idx] + adjustment, 3)
            } else {
              comparison_matrix[i, j] <- round(-nma_result()$TE.nma.random[idx] + adjustment, 3)
            }
          }
        }
      }
    }
    
    # Convert to data frame
    comparison_df <- as.data.frame(comparison_matrix)
    comparison_df$Treatment <- rownames(comparison_df)
    comparison_df <- comparison_df[, c("Treatment", treatments)]
    
    datatable(comparison_df, 
              caption = paste("Treatment Comparisons at Baseline Risk =", baseline_value),
              options = list(scrollX = TRUE))
  })
  
  # Baseline ranking plot
  output$baseline_ranking_plot <- renderPlot({
    req(baseline_result(), nma_result(), input$baseline_ranking_value)
    
    # Generate adjusted rankings
    adjusted_rankings <- generate_adjusted_rankings(
      nma = nma_result(),
      metareg_model = baseline_result()$model,
      covariate_value = input$baseline_ranking_value,
      cov_name = "Baseline Risk"
    )
    
    # Create barplot
    barplot(adjusted_rankings$p_scores * 100, 
            main = paste("P-Scores at Baseline Risk =", input$baseline_ranking_value),
            xlab = "Treatment", ylab = "P-Score (%)",
            col = "lightgreen",
            ylim = c(0, 100),
            names.arg = adjusted_rankings$treatments)
    
    # Add a note
    mtext(paste("Rankings adjusted for baseline risk (", 
                ifelse(input$values_direction == "good", 
                       "smaller values are better", 
                       "larger values are better"), ")", sep = ""),
          side = 3, line = 0.5, cex = 0.8)
  })
  
  # Baseline details
  output$baseline_details <- renderPrint({
    req(baseline_result())
    
    br <- baseline_result()
    
    cat("Baseline Risk Meta-Regression Results:\n")
    cat("=====================================\n\n")
    
    cat("Reference treatment:", br$reference, "\n")
    cat("Number of studies included:", nrow(br$data), "\n\n")
    
    cat("Model Summary:\n")
    print(summary(br$model))
    
    cat("\n\nInterpretation:\n")
    cat("- The intercept represents the effect at baseline risk = 0\n")
    cat("- The slope coefficient indicates how the treatment effect changes\n")
    cat("  with each unit increase in baseline risk\n")
    
    if (br$model$beta[2] > 0) {
      cat("- Positive slope: treatment effect increases with higher baseline risk\n")
    } else if (br$model$beta[2] < 0) {
      cat("- Negative slope: treatment effect decreases with higher baseline risk\n")
    } else {
      cat("- No significant relationship between baseline risk and treatment effect\n")
    }
  })
  
  # Covariate analysis (meta-regression)
  observeEvent(input$run_metareg, {
    req(nma_result(), pairwise_data(), input$covariate_var)
    
    # Input validation
    if (is.null(data()) || nrow(data()) == 0) {
      showNotification("No data available. Please upload or select a dataset first.", type = "error")
      return(NULL)
    }
    
    if (is.null(pairwise_data()) || nrow(pairwise_data()) == 0) {
      showNotification("No pairwise data available. Please run the network meta-analysis first.", type = "error")
      return(NULL)
    }
    
    if (is.null(input$covariate_var) || input$covariate_var == "") {
      showNotification("Please select a covariate for meta-regression.", type = "error")
      return(NULL)
    }
    
    # Create progress indicator
    progress <- shiny::Progress$new()
    progress$set(message = "Running meta-regression...", value = 0.1)
    on.exit(progress$close())
    
    progress$set(message = "Extracting covariates...", value = 0.3)
    
    # Extract the data needed for covariate meta-regression
    tryCatch({
      # Get pairwise data and merge with covariate data
      pw <- pairwise_data()
      
      # Match covariate information to studies
      covariate_name <- input$covariate_var
      pw$matched_covs <- data()[[covariate_name]][match(pw$studlab, data()[[input$study_col]])]
      
      # Check for missing data
      valid_data <- !is.na(pw$TE) & !is.na(pw$seTE) & !is.na(pw$matched_covs)
      if (sum(valid_data) == 0) {
        showNotification("No valid studies for meta-regression. Please check your data and covariate selection.", 
                         type = "error")
        return(NULL)
      }
      
      # Filter valid data
      pw <- pw[valid_data, ]
      
      # After adding the covariate
      if (all(is.na(pw$matched_covs))) {
        showNotification("The selected covariate contains only missing values. Please select another covariate.", type = "error")
        return(NULL)
      }
      
      # For categorical covariates, convert to factor
      if (input$covariate_type == "categorical") {
        pw$matched_covs <- as.factor(pw$matched_covs)
        
        # Update category selection dropdown
        levels <- sort(unique(as.character(pw$matched_covs)))
        if (length(levels) > 0) {
          updateSelectInput(session, "category_value", choices = levels, selected = levels[1])
          updateSelectInput(session, "ranking_category_value", choices = levels, selected = levels[1])
        } else {
          showNotification("No valid categories found for the selected covariate.", type = "error")
          return(NULL)
        }
      }
      
      # Handle treatment-by-covariate interaction
      treatments <- nma_result()$trts
      reference <- if (!is.null(input$interaction_reference)) input$interaction_reference else treatments[1]
      
      progress$set(message = "Running meta-regression model...", value = 0.5)
      
      # Run meta-regression with metafor
      if (input$include_interaction) {
        # Create treatment indicators
        pw$reference <- as.numeric(pw$treat1 == reference | pw$treat2 == reference)
        
        # Create interaction term based on covariate type
        if (input$covariate_type == "categorical") {
          # For categorical covariates
          mr_model <- rma(yi = TE, sei = seTE, mods = ~ matched_covs * reference, 
                          data = pw, method = "REML")
        } else {
          # For numeric covariates
          pw$interaction <- pw$matched_covs * pw$reference
          
          # Model with interaction
          mr_model <- rma(yi = TE, sei = seTE, mods = ~ matched_covs + reference + interaction, 
                          data = pw, method = "REML")
        }
      } else {
        # Simple model without interaction
        mr_model <- rma(yi = TE, sei = seTE, mods = ~ matched_covs, 
                        data = pw, method = "REML")
      }
      
      # Generate predictions at different covariate values
      if (input$covariate_type == "continuous") {
        cov_range <- range(pw$matched_covs, na.rm = TRUE)
        predictions <- calculate_metareg_predictions(mr_model, cov_range, treatments)
      } else {
        # For categorical, predictions are handled differently
        predictions <- NULL
      }
      
      # Store the results
      metareg_result(list(
        covariate = covariate_name,
        type = input$covariate_type,
        data = pw,
        model = mr_model,
        method = "metafor",
        reference_treatment = reference,
        include_interaction = input$include_interaction,
        predictions = predictions
      ))
      
      # Show success notification
      showNotification("Meta-regression completed", type = "message")
      
    }, error = function(e) {
      showNotification(paste("Error in meta-regression:", e$message), type = "error")
    })
    
    progress$set(message = "Meta-regression complete!", value = 1)
  })
  
  # Covariate regression plot
  output$covariate_regression_plot <- renderPlot({
    req(metareg_result())
    
    mr <- metareg_result()
    
    if (!is.null(mr$model)) {
      if (mr$type == "continuous") {
        # For continuous covariates - create a scatter plot with regression line
        # Calculate weights based on precision (1/variance)
        weights <- 1 / (mr$data$seTE^2)
        point_sizes <- weights / max(weights) * 5 + 1  # Scale for point sizes
        
        # Create a scatter plot with varying point sizes based on precision
        plot(mr$data$matched_covs, mr$data$TE, 
             xlab = mr$covariate, 
             ylab = paste("Effect (", input$effect_measure, ")", sep = ""),
             main = paste("Meta-regression with", mr$covariate, 
                          ifelse(mr$include_interaction, 
                                 paste("\nWith interaction using", mr$reference_treatment, "as reference"),
                                 "")),
             pch = 19, col = "blue", cex = point_sizes)
        
        # Add legend for point sizes
        legend("topright", 
               legend = c("Larger points = higher precision", 
                          "Solid line = regression", 
                          "Dashed lines = 95% CI"),
               bty = "n", cex = 0.8)
        
        # Add regression line
        if (!mr$include_interaction) {
          # Simple regression line
          abline(mr$model$beta[1], mr$model$beta[2], col = "red", lwd = 2)
        } else {
          # For interaction models, show the average effect
          x_range <- range(mr$data$matched_covs, na.rm = TRUE)
          x_seq <- seq(x_range[1], x_range[2], length.out = 100)
          y_seq <- mr$model$beta[1] + mr$model$beta[2] * x_seq
          lines(x_seq, y_seq, col = "red", lwd = 2)
        }
        
        # Add horizontal line at y = 0
        abline(h = 0, lty = 2)
        
        # Add confidence bands
        x_range <- range(mr$data$matched_covs, na.rm = TRUE)
        newdata <- data.frame(matched_covs = seq(x_range[1], x_range[2], length.out = 100))
        
        tryCatch({
          if (!mr$include_interaction) {
            pred <- predict(mr$model, newmods = cbind(newdata$matched_covs), level = 0.95)
          } else {
            # For interaction model, predict at average reference value
            newdata$reference <- mean(mr$data$reference)
            newdata$interaction <- newdata$matched_covs * newdata$reference
            pred <- predict(mr$model, newmods = cbind(newdata$matched_covs, newdata$reference, newdata$interaction), level = 0.95)
          }
          
          # Add the prediction interval
          lines(newdata$matched_covs, pred$pred, col = "red", lwd = 2)
          lines(newdata$matched_covs, pred$ci.lb, col = "red", lty = 2)
          lines(newdata$matched_covs, pred$ci.ub, col = "red", lty = 2)
        }, error = function(e) {
          # Continue without confidence bands if prediction fails
        })
      } else {
        # For categorical covariates - create a box plot
        # Convert matched_covs to factor if it's not already
        factor_data <- factor(mr$data$matched_covs)
        
        # Create boxplot
        boxplot(mr$data$TE ~ factor_data,
                xlab = mr$covariate,
                ylab = paste("Effect (", input$effect_measure, ")", sep = ""),
                main = paste("Meta-regression with", mr$covariate, "(Categorical)"),
                col = "lightblue")
        
        # Add horizontal line at y = 0
        abline(h = 0, lty = 2, col = "red")
        
        # Add sample sizes
        n_per_category <- table(factor_data)
        mtext(paste("n =", n_per_category), at = 1:length(n_per_category), 
              side = 3, line = 0, cex = 0.8)
      }
    } else {
      # Create placeholder plot
      plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(-2, 2),
           xlab = input$covariate_var, ylab = "Effect",
           main = paste("Meta-Regression with", input$covariate_var, "\n(Run analysis first)"))
      text(0.5, 0, "Run meta-regression analysis", cex = 1.5)
    }
  })
  
  # Covariate forest plot
  output$covariate_forest_plot <- renderPlot({
    req(metareg_result(), nma_result())
    
    mr <- metareg_result()
    nma <- nma_result()
    
    if (mr$type == "continuous") {
      # For continuous covariates, use the slider value
      covariate_value <- input$covariate_value
    } else {
      # For categorical covariates, use the selected category
      covariate_value <- input$category_value
    }
    
    # Use the helper function to create a covariate-adjusted forest plot
    create_metareg_forest_plot(
      nma = nma,
      metareg_model = mr$model,
      covariate_value = covariate_value,
      reference_treatment = mr$reference_treatment,
      sm = input$effect_measure,
      cov_name = mr$covariate
    )
  })
  
  # Covariate comparison table
  output$covariate_comparison_table <- renderDT({
    req(metareg_result())
    
    mr <- metareg_result()
    
    # Extract model results
    results <- summary(mr$model)
    
    # Create data frame of results
    result_df <- data.frame(
      Parameter = rownames(results$beta),
      Estimate = round(results$beta, 4),
      SE = round(results$se, 4),
      Z_value = round(results$zval, 3),
      P_value = format.pval(results$pval, digits = 3),
      CI_lower = round(results$ci.lb, 4),
      CI_upper = round(results$ci.ub, 4),
      stringsAsFactors = FALSE
    )
    
    # Add row names that are more descriptive
    if (mr$include_interaction) {
      if (nrow(result_df) >= 4) {
        result_df$Parameter[1] <- "Intercept"
        result_df$Parameter[2] <- mr$covariate
        result_df$Parameter[3] <- paste("Reference treatment (", mr$reference_treatment, ")")
        result_df$Parameter[4] <- paste(mr$covariate, "× Reference interaction")
      }
    } else {
      if (nrow(result_df) >= 2) {
        result_df$Parameter[1] <- "Intercept"
        result_df$Parameter[2] <- mr$covariate
      }
    }
    
    datatable(result_df, 
              caption = paste("Meta-Regression Results:", mr$covariate),
              options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # Covariate ranking plot
  output$covariate_ranking_plot <- renderPlot({
    req(nma_result())
    
    if (!is.null(metareg_result())) {
      mr <- metareg_result()
      
      # Get covariate value based on type
      if (mr$type == "continuous") {
        covariate_value <- input$ranking_covariate_value
      } else {
        covariate_value <- input$ranking_category_value
      }
      
      # Generate adjusted rankings using the helper function
      adjusted_rankings <- generate_adjusted_rankings(
        nma = nma_result(),
        metareg_model = mr$model,
        covariate_value = covariate_value,
        cov_name = mr$covariate
      )
      
      # Create barplot with adjusted P-scores
      barplot(adjusted_rankings$p_scores * 100, 
              main = paste("Adjusted P-Score Rankings at", mr$covariate, "=", round(covariate_value, 2),
                           "\n(", ifelse(input$values_direction == "good", 
                                         "smaller values are better", 
                                         "larger values are better"), ")"),
              xlab = "Treatment", ylab = "Adjusted P-Score (%)",
              col = "lightgreen",
              ylim = c(0, 100),
              names.arg = adjusted_rankings$treatments,
              las = 2)
      
      # Add a note about the adjustment
      if (mr$include_interaction) {
        mtext(paste("Rankings adjusted for covariate-by-treatment interaction",
                    "using", mr$reference_treatment, "as reference"), 
              side = 3, line = 0.5, cex = 0.7)
      } else {
        mtext("Rankings adjusted for covariate effect", 
              side = 3, line = 0.5, cex = 0.7)
      }
    } else {
      # Show placeholder asking to run meta-regression
      plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
           xlab = "", ylab = "",
           main = "Covariate-Adjusted Rankings")
      text(0.5, 0.5, "Run meta-regression first", cex = 1.5)
    }
  })
  
  # Bubble plot for meta-regression
  output$metareg_bubble_plot <- renderPlot({
    req(metareg_result(), input$bubble_treatment)
    
    mr <- metareg_result()
    
    # Filter data for the selected treatment
    bubble_data <- mr$data %>%
      filter(treat1 == input$bubble_treatment | treat2 == input$bubble_treatment)
    
    # Create placeholder if not enough data
    if (nrow(bubble_data) < 2) {
      plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
           xlab = "", ylab = "",
           main = paste("Bubble Plot for", input$bubble_treatment))
      text(0.5, 0.5, "Not enough data points for this treatment", cex = 1.2)
      return()
    }
    
    # Create bubble plot using the helper function
    print(create_metareg_bubble_plot(
      data = bubble_data,
      xvar = "matched_covs",
      yvar = "TE",
      sizevar = "seTE",
      title = paste("Meta-Regression Bubble Plot for", input$bubble_treatment),
      xlab = mr$covariate,
      ylab = paste("Effect (", input$effect_measure, ")")
    ))
  })
  
  # Meta-regression model summary
  output$metareg_model_summary <- renderPrint({
    req(metareg_result())
    
    mr <- metareg_result()
    
    cat("Meta-Regression Model Summary:\n")
    cat("==============================\n\n")
    
    print(summary(mr$model))
    
    cat("\nInterpretation:\n")
    cat("================\n")
    
    if (mr$type == "continuous") {
      cat("- The intercept represents the effect when", mr$covariate, "= 0\n")
      cat("- The coefficient for", mr$covariate, "indicates how the effect changes\n")
      cat("  with each unit increase in", mr$covariate, "\n")
      
      if (mr$model$beta[2] > 0) {
        cat("- Positive coefficient: treatment effect increases with higher", mr$covariate, "\n")
      } else if (mr$model$beta[2] < 0) {
        cat("- Negative coefficient: treatment effect decreases with higher", mr$covariate, "\n")
      }
    } else {
      cat("- The intercept represents the effect for the reference category\n")
      cat("- Other coefficients represent the difference from the reference category\n")
    }
    
    if (mr$include_interaction) {
      cat("\nNote: This model includes treatment-by-covariate interaction using\n")
      cat(mr$reference_treatment, "as the reference treatment.\n")
      cat("The interaction term indicates whether the covariate effect differs\n")
      cat("between treatments involving the reference and other comparisons.\n")
    }
  })
  
  # Predicted effects table
  output$predicted_effects_table <- renderDT({
    req(metareg_result())
    
    mr <- metareg_result()
    
    if (mr$type == "continuous" && !is.null(mr$predictions)) {
      # Create a table of predicted effects at different covariate values
      pred_data <- do.call(rbind, mr$predictions)
      
      # Format the data frame
      pred_data <- pred_data %>%
        mutate(
          Covariate_Value = round(covariate_value, 2),
          Effect = round(effect, 3),
          Lower_CI = round(lower_ci, 3),
          Upper_CI = round(upper_ci, 3),
          CI = paste0("[", Lower_CI, ", ", Upper_CI, "]")
        ) %>%
        select(Covariate_Value, Effect, CI)
      
      datatable(pred_data,
                caption = paste("Predicted Effects at Different Values of", mr$covariate),
                options = list(pageLength = 10, scrollX = TRUE))
    } else if (mr$type == "categorical") {
      # For categorical variables, show the effect for each category
      coef_summary <- summary(mr$model)
      
      # Extract category effects
      categories <- levels(factor(mr$data$matched_covs))
      effects_df <- data.frame(
        Category = categories,
        Effect = NA,
        CI = NA,
        stringsAsFactors = FALSE
      )
      
      # Reference category
      effects_df$Effect[1] <- round(coef_summary$beta[1], 3)
      effects_df$CI[1] <- paste0("[", round(coef_summary$ci.lb[1], 3), ", ", 
                                 round(coef_summary$ci.ub[1], 3), "]")
      
      # Other categories
      if (length(categories) > 1) {
        for (i in 2:length(categories)) {
          if (i <= nrow(coef_summary$beta)) {
            effects_df$Effect[i] <- round(coef_summary$beta[1] + coef_summary$beta[i], 3)
            # Approximate CI (this is simplified)
            effects_df$CI[i] <- paste0("[", 
                                       round(coef_summary$beta[1] + coef_summary$beta[i] - 1.96 * coef_summary$se[i], 3), 
                                       ", ", 
                                       round(coef_summary$beta[1] + coef_summary$beta[i] + 1.96 * coef_summary$se[i], 3), 
                                       "]")
          }
        }
      }
      
      datatable(effects_df,
                caption = paste("Predicted Effects by", mr$covariate, "Category"),
                options = list(pageLength = 10, scrollX = TRUE))
    } else {
      # Return empty table if no predictions
      datatable(data.frame(Message = "No predictions available"),
                caption = "Predicted Effects",
                options = list(dom = 't'))
    }
  })
  
  # Model diagnostics - Goodness of fit
  output$goodness_of_fit <- renderPrint({
    req(metareg_result())
    
    mr <- metareg_result()
    
    cat("Goodness of Fit for Meta-Regression Model:\n")
    cat("==========================================\n\n")
    
    # Extract model fit statistics
    cat("Model Fit Statistics:\n")
    cat("---------------------\n")
    cat("AIC:", round(AIC(mr$model), 2), "\n")
    cat("BIC:", round(BIC(mr$model), 2), "\n")
    cat("Log-likelihood:", round(logLik(mr$model), 2), "\n\n")
    
    # Test of moderators
    cat("Test of Moderators (Omnibus Test):\n")
    cat("-----------------------------------\n")
    cat("QM =", round(mr$model$QM, 2), "on", mr$model$m, "df, p =", 
        format.pval(mr$model$QMp, digits = 3), "\n")
    
    if (mr$model$QMp < 0.05) {
      cat("The moderator(s) significantly explain heterogeneity in the effect sizes.\n\n")
    } else {
      cat("The moderator(s) do not significantly explain heterogeneity in the effect sizes.\n\n")
    }
    
    # Variance components 
    cat("Variance Components:\n")
    cat("--------------------\n")
    cat("Tau² (estimated amount of residual heterogeneity):", round(mr$model$tau2, 4), "\n")
    cat("I² (residual heterogeneity / unaccounted variability):", 
        round(mr$model$I2, 1), "%\n")
    cat("H² (unaccounted variability / sampling variability):", round(mr$model$H2, 2), "\n")
    
    # R² statistic
    if (!is.null(mr$model$R2)) {
      cat("\nProportion of Heterogeneity Explained:\n")
      cat("--------------------------------------\n")
      cat("R² (amount of heterogeneity accounted for):", round(mr$model$R2, 1), "%\n")
    }
  })
  
  # Model diagnostics - Residual plot
  output$residual_plot <- renderPlot({
    req(metareg_result())
    
    mr <- metareg_result()
    
    # Get residuals and fitted values
    fitted_vals <- fitted(mr$model)
    resid_vals <- residuals(mr$model)
    
    # Create a residual plot
    plot(fitted_vals, resid_vals, 
         xlab = "Fitted Values", ylab = "Residuals",
         main = "Residuals vs Fitted Values",
         pch = 19, col = "blue")
    
    # Add a horizontal line at y = 0
    abline(h = 0, lty = 2, col = "red", lwd = 2)
    
    # Add a smoothed line to show any patterns
    if (length(fitted_vals) > 10) {
      lines(lowess(fitted_vals, resid_vals), col = "green", lwd = 2)
      legend("topright", 
             legend = c("Residuals", "Zero line", "Lowess smooth"),
             col = c("blue", "red", "green"),
             lty = c(NA, 2, 1),
             pch = c(19, NA, NA),
             lwd = c(NA, 2, 2),
             bty = "n")
    }
    
    # Add text about patterns
    if (length(fitted_vals) > 10) {
      # Simple check for heteroscedasticity
      cor_test <- cor.test(abs(resid_vals), fitted_vals, method = "spearman")
      if (cor_test$p.value < 0.05) {
        mtext("Warning: Potential heteroscedasticity detected", 
              side = 3, line = 0.5, col = "red", cex = 0.9)
      }
    }
  })
  
  # Model diagnostics - Q-Q plot
  output$qq_plot <- renderPlot({
    req(metareg_result())
    
    mr <- metareg_result()
    
    # Get standardized residuals
    res <- rstudent(mr$model)
    
    # Create a Q-Q plot
    qqnorm(res$z, main = "Normal Q-Q Plot of Standardized Residuals",
           pch = 19, col = "blue")
    qqline(res$z, col = "red", lwd = 2)
    
    # Add confidence bands
    n <- length(res$z)
    if (n > 10) {
      # Calculate approximate confidence bands
      sorted_res <- sort(res$z)
      theoretical_quantiles <- qnorm(ppoints(n))
      
      # Standard error for order statistics
      se <- sqrt(ppoints(n) * (1 - ppoints(n)) / n) / dnorm(theoretical_quantiles)
      
      # 95% confidence bands
      upper_band <- theoretical_quantiles + 1.96 * se
      lower_band <- theoretical_quantiles - 1.96 * se
      
      lines(theoretical_quantiles, upper_band, lty = 2, col = "gray")
      lines(theoretical_quantiles, lower_band, lty = 2, col = "gray")
      
      legend("topleft", 
             legend = c("Data", "Expected", "95% CI"),
             col = c("blue", "red", "gray"),
             lty = c(NA, 1, 2),
             pch = c(19, NA, NA),
             lwd = c(NA, 2, 1),
             bty = "n")
    }
    
    # Test for normality
    if (n > 7) {
      shapiro_test <- shapiro.test(res$z)
      mtext(paste("Shapiro-Wilk test: p =", format.pval(shapiro_test$p.value, digits = 3)), 
            side = 3, line = 0.5, cex = 0.9)
    }
  })
  
  # Download handlers for all meta-regression outputs
  output$download_covariate_comparison <- downloadHandler(
    filename = function() {
      paste0("metareg_comparison_", gsub(" ", "_", input$covariate_var), ".csv")
    },
    content = function(file) {
      req(metareg_result())
      
      mr <- metareg_result()
      
      # Extract model results
      results <- summary(mr$model)
      
      # Create data frame of results
      result_df <- data.frame(
        Parameter = rownames(results$beta),
        Estimate = results$beta,
        SE = results$se,
        Z_value = results$zval,
        P_value = results$pval,
        CI_lower = results$ci.lb,
        CI_upper = results$ci.ub,
        stringsAsFactors = FALSE
      )
      
      write.csv(result_df, file, row.names = FALSE)
    }
  )
  
  output$download_covariate_ranking <- downloadHandler(
    filename = function() {
      paste0("metareg_rankings_", gsub(" ", "_", input$covariate_var), ".csv")
    },
    content = function(file) {
      req(nma_result(), metareg_result())
      
      mr <- metareg_result()
      
      # Get covariate value based on type
      if (mr$type == "continuous") {
        covariate_value <- input$ranking_covariate_value
      } else {
        covariate_value <- input$ranking_category_value
      }
      
      # Generate adjusted rankings using the helper function
      adjusted_rankings <- generate_adjusted_rankings(
        nma = nma_result(),
        metareg_model = mr$model,
        covariate_value = covariate_value,
        cov_name = mr$covariate
      )
      
      # Create ranking dataframe
      ranking_df <- data.frame(
        Treatment = adjusted_rankings$treatments,
        Adjusted_P_Score = round(adjusted_rankings$p_scores * 100, 1),
        Covariate = mr$covariate,
        Covariate_Value = covariate_value,
        stringsAsFactors = FALSE
      )
      
      write.csv(ranking_df, file, row.names = FALSE)
    }
  )
  
  output$download_predicted_effects <- downloadHandler(
    filename = function() {
      paste0("metareg_predictions_", gsub(" ", "_", input$covariate_var), ".csv")
    },
    content = function(file) {
      req(metareg_result())
      
      mr <- metareg_result()
      
      if (mr$type == "continuous" && !is.null(mr$predictions)) {
        # Create a table of predicted effects at different covariate values
        pred_data <- do.call(rbind, mr$predictions)
        write.csv(pred_data, file, row.names = FALSE)
      } else {
        # Return empty table if no predictions
        write.csv(data.frame(Message = "No predictions available"), file, row.names = FALSE)
      }
    }
  )
  
  output$download_baseline_comparison <- downloadHandler(
    filename = function() {
      "baseline_risk_comparison.csv"
    },
    content = function(file) {
      req(baseline_result(), nma_result())
      
      # Get all treatments
      treatments <- nma_result()$trts
      n_treats <- length(treatments)
      
      # Create comparison matrix
      comparison_matrix <- matrix(NA, nrow = n_treats, ncol = n_treats,
                                  dimnames = list(treatments, treatments))
      
      # Get baseline risk value
      baseline_value <- input$baseline_value
      
      # Calculate adjusted effects
      br <- baseline_result()
      adjustment <- 0
      
      if (!is.null(br$model)) {
        # Calculate adjustment based on baseline risk
        mean_baseline <- mean(br$data$baseline_risk, na.rm = TRUE)
        adjustment <- br$model$beta[2] * (baseline_value - mean_baseline)
      }
      
      # Fill comparison matrix with adjusted NMA estimates
      for (i in 1:n_treats) {
        for (j in 1:n_treats) {
          if (i != j) {
            # Find the comparison in NMA results
            idx <- which((nma_result()$treat1 == treatments[i] & nma_result()$treat2 == treatments[j]) |
                         (nma_result()$treat1 == treatments[j] & nma_result()$treat2 == treatments[i]))
            
            if (length(idx) > 0) {
              if (nma_result()$treat1[idx] == treatments[i]) {
                comparison_matrix[i, j] <- round(nma_result()$TE.nma.random[idx] + adjustment, 3)
              } else {
                comparison_matrix[i, j] <- round(-nma_result()$TE.nma.random[idx] + adjustment, 3)
              }
            }
          }
        }
      }
      
      # Convert to data frame
      comparison_df <- as.data.frame(comparison_matrix)
      comparison_df$Treatment <- rownames(comparison_df)
      comparison_df <- comparison_df[, c("Treatment", treatments)]
      
      write.csv(comparison_df, file, row.names = FALSE)
    }
  )
  
  output$download_baseline_ranking <- downloadHandler(
    filename = function() {
      paste0("baseline_risk_rankings_", input$baseline_ranking_value, ".csv")
    },
    content = function(file) {
      req(baseline_result(), nma_result(), input$baseline_ranking_value)
      
      # Generate adjusted rankings
      adjusted_rankings <- generate_adjusted_rankings(
        nma = nma_result(),
        metareg_model = baseline_result()$model,
        covariate_value = input$baseline_ranking_value,
        cov_name = "Baseline Risk"
      )
      
      # Create ranking dataframe
      ranking_df <- data.frame(
        Treatment = adjusted_rankings$treatments,
        P_Score = round(adjusted_rankings$p_scores * 100, 1),
        Baseline_Risk = input$baseline_ranking_value,
        stringsAsFactors = FALSE
      )
      
      write.csv(ranking_df, file, row.names = FALSE)
    }
  )
  
  # Download baseline forest plot
  output$download_baseline_forest <- downloadHandler(
    filename = function() {
      "baseline_forest_plot.png"
    },
    content = function(file) {
      req(baseline_result(), nma_result(), input$baseline_value)
      
      # Save the plot
      png(file, width = 800, height = 600, res = 100)
      
      # Create adjusted forest plot
      create_metareg_forest_plot(
        nma = nma_result(),
        metareg_model = baseline_result()$model,
        covariate_value = input$baseline_value,
        reference_treatment = baseline_result()$reference,
        sm = input$effect_measure,
        cov_name = "Baseline Risk"
      )
      
      dev.off()
    }
  )
  
  # Download covariate forest plot
  output$download_covariate_forest <- downloadHandler(
    filename = function() {
      paste0("covariate_forest_plot_", gsub(" ", "_", input$covariate_var), ".png")
    },
    content = function(file) {
      req(metareg_result(), nma_result())
      
      # Save the plot
      png(file, width = 800, height = 600, res = 100)
      
      mr <- metareg_result()
      nma <- nma_result()
      
      if (mr$type == "continuous") {
        # For continuous covariates, use the slider value
        covariate_value <- input$covariate_value
      } else {
        # For categorical covariates, use the selected category
        covariate_value <- input$category_value
      }
      
      # Use the helper function to create a covariate-adjusted forest plot
      create_metareg_forest_plot(
        nma = nma,
        metareg_model = mr$model,
        covariate_value = covariate_value,
        reference_treatment = mr$reference_treatment,
        sm = input$effect_measure,
        cov_name = mr$covariate
      )
      
      dev.off()
    }
  )
}

# Run the application
shinyApp(ui = ui, server = server)
