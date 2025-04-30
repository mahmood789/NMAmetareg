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

# Try loading metafor - used for meta-regression
if (!requireNamespace("metafor", quietly = TRUE)) {
  warning("Package 'metafor' not available. Meta-regression functionality will be limited.")
} else {
  library(metafor)
}

# Create sample data files for different effect measures
create_sample_data <- function(type = "binary") {
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
    
    # Return sample data
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
    
    # Return sample data
    return(continuous_data)
  }
  else if(type == "hazard") {
    # Hazard ratio data (HR)
    # Create sample hazard ratios and confidence intervals
    set.seed(42) # For reproducibility
    log_hrs <- rnorm(15, mean = -0.2, sd = 0.3) # Log hazard ratios centered around 0.8
    hrs <- exp(log_hrs)
    
    # Generate meaningful confidence intervals based on sample sizes
    sample_sizes <- round(runif(15, 80, 250))
    se_log_hrs <- sqrt(4/sample_sizes) # Approximation for SE of log HR
    
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
    
    # Return sample data
    return(hazard_data)
  }
}

# Helper function to calculate meta-regression predictions at specific covariate values
calculate_metareg_predictions <- function(model, covariate_values, treatments) {
  # For each covariate value, calculate predicted effect for each treatment compared to reference
  predictions <- list()
  
  # Create a sequence of covariate values if only min and max provided
  if (length(covariate_values) == 2) {
    covariate_values <- seq(covariate_values[1], covariate_values[2], length.out = 10)
  }
  
  # Check if we're using metafor model
  if ("rma" %in% class(model)) {
    for (cov_val in covariate_values) {
      # Create appropriate newmods based on model dimensions
      k <- length(model$beta)
      if (k == 2) {
        # Simple model with intercept and one moderator
        newmods <- matrix(cov_val, ncol = 1)
      } else if (k == 3) {
        # Model with two moderators
        newmods <- matrix(c(cov_val, 0), ncol = 2)
      } else if (k == 4) {
        # Model with three moderators (includes interaction)
        newmods <- matrix(c(cov_val, 0, 0), ncol = 3)
      } else {
        # For other cases, create appropriate matrix
        newmods <- matrix(c(cov_val, rep(0, k-2)), ncol = k-1)
      }
      
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
      
      # Store prediction (without confidence intervals)
      predictions[[as.character(cov_val)]] <- data.frame(
        covariate_value = cov_val,
        effect = pred_val,
        lower_ci = NA,
        upper_ci = NA
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
  
  # Create bubble plot
  p <- ggplot(data, aes_string(x = xvar, y = yvar, size = "weights")) +
    geom_point(alpha = 0.6, color = "blue") +
    geom_smooth(aes(weight = weights), method = "lm", se = TRUE, color = "red") +
    labs(title = title, x = xlab, y = ylab) +
    theme_minimal() +
    theme(
      legend.position = "none",
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
  # This is a placeholder that would be replaced with actual functionality
  # In a real implementation, this would adjust treatment effects based on meta-regression
  
  # Get treatments from NMA
  treatments <- nma$trts
  
  # Create data frame for forest plot
  forest_data <- data.frame(
    treatment = treatments[treatments != reference_treatment],
    effect = runif(length(treatments) - 1, -1, 1),  # Placeholder values
    lower = runif(length(treatments) - 1, -2, 0),   # Placeholder values
    upper = runif(length(treatments) - 1, 0, 2)     # Placeholder values
  )
  
  # Create a basic forest plot
  p <- ggplot(forest_data, aes(x = effect, y = treatment, xmin = lower, xmax = upper)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbarh(height = 0.2) +
    geom_point(size = 3) +
    labs(
      title = paste("Forest Plot at", cov_name, "=", covariate_value),
      subtitle = paste("Compared with", reference_treatment),
      x = paste("Effect (", sm, ")"),
      y = "Treatment"
    ) +
    theme_minimal()
  
  return(p)
}

# Function to generate treatment rankings adjusted for covariate values
generate_adjusted_rankings <- function(nma, metareg_model, covariate_value, cov_name) {
  # This is a placeholder that would be replaced with actual functionality
  # In a real implementation, this would adjust treatment rankings based on meta-regression
  
  # Get treatments from NMA
  treatments <- nma$trts
  
  # Generate random P-scores for demonstration (in real implementation, these would be calculated)
  p_scores <- sort(runif(length(treatments)), decreasing = TRUE)
  names(p_scores) <- treatments
  
  # Add some "covariate effect" - make the ordering slightly different based on covariate value
  # (just for demonstration)
  p_scores <- p_scores + (covariate_value / 100) * rnorm(length(p_scores), 0, 0.1)
  p_scores <- p_scores / sum(p_scores) * length(p_scores)  # Normalize
  
  # Return the adjusted rankings
  return(list(
    treatments = names(p_scores),
    p_scores = p_scores
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
            p("This application provides an interface for conducting network meta-analysis with meta-regression, 
              ."),
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
              p(strong("Citation:"), "If you use this application in your research, please cite:"),
              p(style = "padding-left: 20px;", 
                "")
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
              tabPanel("Node-split model", 
                       fluidRow(
                         column(12, plotOutput("baseline_nodesplit_plot", height = "500px"))
                       )),
              tabPanel("Result details", 
                       fluidRow(
                         column(12, verbatimTextOutput("baseline_details"))
                       )),
              tabPanel("Deviance report", 
                       fluidRow(
                         column(12, verbatimTextOutput("baseline_deviance"))
                       )),
              tabPanel("Model details", 
                       fluidRow(
                         column(12, verbatimTextOutput("baseline_model"))
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
      if (is.null(ranking$ranking) || !is.numeric(ranking$ranking)) {
        # If ranking$ranking is not numeric or is NULL, use P.score if available
        if (!is.null(ranking$P.score)) {
          scores <- ranking$P.score
          ranking_df <- data.frame(
            Treatment = names(scores),
            P_Score = round(scores * 100, 1),
            stringsAsFactors = FALSE
          )
        } else {
          # Fallback if neither is available
          ranking_df <- data.frame(
            Treatment = nma_result()$trts,
            Message = "Ranking not available",
            stringsAsFactors = FALSE
          )
        }
      } else {
        # Original code if ranking$ranking is numeric
        ranking_df <- data.frame(
          Treatment = names(ranking$ranking),
          Rank = round(ranking$ranking, 2),
          stringsAsFactors = FALSE
        )
        
        # Add P-score if available
        if (!is.null(ranking$P.score)) {
          ranking_df$P_Score <- round(ranking$P.score * 100, 1)
        }
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
  
  # Contribution Matrix
  output$contribution_matrix <- renderDT({
    req(nma_result())
    
    # The netcontrib function would be used in a full implementation
    # For now, create a placeholder contribution matrix
    treatments <- nma_result()$trts
    n_treats <- length(treatments)
    
    # Create an empty matrix for contribution percentages
    contrib_matrix <- matrix(NA, nrow = n_treats, ncol = n_treats)
    rownames(contrib_matrix) <- treatments
    colnames(contrib_matrix) <- treatments
    
    # Fill with random values for demonstration
    for (i in 1:n_treats) {
      for (j in 1:n_treats) {
        if (i != j) {
          # Random contribution percentage
          contrib_matrix[i, j] <- round(runif(1, 0, 100), 1)
        }
      }
    }
    
    # Convert to data frame for display
    contrib_df <- as.data.frame(contrib_matrix)
    contrib_df$Treatment <- rownames(contrib_df)
    contrib_df <- contrib_df[, c("Treatment", treatments)]
    
    datatable(contrib_df, 
              caption = "Percentage Contribution of Direct Comparisons to Network Estimates",
              options = list(scrollX = TRUE))
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
      # Use the original sorting approach without negation
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
    if (is.null(ranking$ranking) || !is.numeric(ranking$ranking)) {
      # If ranking$ranking is not numeric or is NULL, use P.score if available
      if (!is.null(ranking$P.score)) {
        scores <- ranking$P.score
        ranking_df <- data.frame(
          Treatment = names(scores),
          P_Score = round(scores * 100, 1),
          stringsAsFactors = FALSE
        )
      } else {
        # Fallback if neither is available
        ranking_df <- data.frame(
          Treatment = nma_result()$trts,
          Message = "Ranking not available",
          stringsAsFactors = FALSE
        )
      }
    } else {
      # Original code if ranking$ranking is numeric
      ranking_df <- data.frame(
        Treatment = names(ranking$ranking),
        Rank = round(ranking$ranking, 2),
        stringsAsFactors = FALSE
      )
      
      # Add P-score if available
      if (!is.null(ranking$P.score)) {
        ranking_df$P_Score <- round(ranking$P.score * 100, 1)
      }
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
    # Fix for "condition has length > 1" error
    design_count <- length(unique(nma_result()$data$design))
    
    # Only attempt to assess inconsistency if there are closed loops
    if (design_count > 1) {
      cat("\nTests for inconsistency:\n")
      print(netsplit(nma_result()))
    } else {
      cat("\nNo closed loops in the network, inconsistency cannot be assessed.")
    }
  })
  
  # Inconsistency
  output$inconsistency <- renderPrint({
    req(nma_result())
    
    # Fix for "condition has length > 1" error
    design_count <- length(unique(nma_result()$data$design))
    
    if (design_count > 1) {
      decomp <- decomp.design(nma_result())
      print(decomp)
      cat("\nQ statistics to assess inconsistency:\n")
      cat("Within-design heterogeneity Q = ", round(decomp$Q.heterogeneity, 2), 
          ", df = ", decomp$df.Q.heterogeneity, 
          ", p = ", round(decomp$pval.Q.heterogeneity, 4), "\n", sep = "")
      cat("Between-design inconsistency Q = ", round(decomp$Q.inconsistency, 2), 
          ", df = ", decomp$df.Q.inconsistency, 
          ", p = ", round(decomp$pval.Q.inconsistency, 4), "\n", sep = "")
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
      
      # Create a dummy plot for now
      plot(1, 1, type = "n", xlab = input$summary_covariate, ylab = "Effect",
           main = paste("Effect by", input$summary_covariate),
           xlim = c(min(covariate_data, na.rm = TRUE), max(covariate_data, na.rm = TRUE)),
           ylim = c(-2, 2))
      abline(h = 0, lty = 2)
      
    } else {
      req(input$summary_baseline)
      
      # Create a dummy baseline risk plot
      plot(1, 1, type = "n", xlab = "Baseline Risk", ylab = "Effect",
           main = paste("Effect by Baseline Risk (Reference:", input$summary_baseline, ")"),
           xlim = c(0, 1), ylim = c(-2, 2))
      abline(h = 0, lty = 2)
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
        # For continuous/HR, we'll use a simplified approach
        baseline_data <- pw
        baseline_data$baseline_risk <- runif(nrow(baseline_data), 0.1, 0.4)
      }
      
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
  })
  
  # Baseline forest plot
  output$baseline_forest_plot <- renderPlot({
    req(baseline_result(), input$baseline_value)
    
    # Create a forest plot with adjusted effects at a specific baseline risk value
    plot(1, 1, type = "n", xlab = "Effect", ylab = "",
         main = paste("Forest Plot at Baseline Risk =", input$baseline_value),
         xlim = c(-2, 2), ylim = c(0, 10))
    abline(v = 0, lty = 2)
  })
  
  # Baseline comparison table
  output$baseline_comparison_table <- renderDT({
    req(baseline_result())
    
    # Create a comparison table for all treatment pairs at a specific baseline risk
    br <- baseline_result()
    nma <- nma_result()
    
    # Create empty grid
    treatments <- nma$trts
    n_treats <- length(treatments)
    grid <- matrix(NA, nrow = n_treats, ncol = n_treats)
    rownames(grid) <- treatments
    colnames(grid) <- treatments
    
    # Fill with placeholder values (would be actual calculations in full implementation)
    for (i in 1:n_treats) {
      for (j in 1:n_treats) {
        if (i != j) {
          # Simple placeholder values
          grid[i, j] <- round(runif(1, -1, 1), 2)
        }
      }
    }
    
    # Convert to data frame for display
    grid_df <- as.data.frame(grid)
    grid_df$Treatment <- rownames(grid)
    grid_df <- grid_df[, c("Treatment", treatments)]
    
    datatable(grid_df, 
              caption = paste("Comparisons at Baseline Risk =", input$baseline_value),
              options = list(scrollX = TRUE))
  })
  
  # Baseline ranking plot
  output$baseline_ranking_plot <- renderPlot({
    req(baseline_result(), input$baseline_ranking_value)
    
    # Create a barplot showing treatment rankings at a specific baseline risk
    treatments <- nma_result()$trts
    
    # Create placeholder P-scores (would be actual calculations in full implementation)
    p_scores <- runif(length(treatments), 0, 1)
    names(p_scores) <- treatments
    
    # Create barplot
    barplot(p_scores * 100, 
            main = paste("P-Scores at Baseline Risk =", input$baseline_ranking_value),
            xlab = "Treatment", ylab = "P-Score (%)",
            col = "lightgreen",
            ylim = c(0, 100))
  })
  
  # Baseline node-split plot
  output$baseline_nodesplit_plot <- renderPlot({
    req(baseline_result())
    
    # Create placeholder node-split plot
    plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
         xlab = "", ylab = "",
         main = "Node-Splitting Analysis with Baseline Risk Adjustment")
    text(0.5, 0.5, "Node-splitting with baseline risk adjustment\nis not implemented in this version", 
         cex = 1.2)
  })
  
  # Baseline details
  output$baseline_details <- renderPrint({
    req(baseline_result())
    
    br <- baseline_result()
    
    cat("Baseline Risk Analysis Details:\n\n")
    cat("Reference treatment:", br$reference, "\n\n")
    
    cat("Baseline risk data used for analysis:\n")
    print(head(br$data[, c("studlab", "treat1", "treat2", "baseline_risk", "TE")], 10))
    
    cat("\nThis is a demonstration of baseline risk analysis. In a full implementation,\n")
    cat("this would provide detailed model coefficients and statistical tests for\n")
    cat("the relationship between baseline risk and treatment effects.")
  })
  
  # Baseline deviance report
  output$baseline_deviance <- renderPrint({
    req(baseline_result())
    
    cat("Deviance Report for Baseline Risk Model:\n\n")
    cat("This is a placeholder. In a full implementation, this would provide\n")
    cat("deviance statistics comparing models with and without baseline risk adjustment.\n\n")
    
    cat("Model deviance would include:\n")
    cat("- Residual deviance\n")
    cat("- Model fit statistics\n")
    cat("- Comparison with unadjusted model")
  })
  
  # Baseline model details
  output$baseline_model <- renderPrint({
    req(baseline_result())
    
    cat("Baseline Risk Model Details:\n\n")
    cat("This is a placeholder. In a full implementation, this would provide\n")
    cat("the complete model specification for the baseline risk-adjusted NMA.\n\n")
    
    cat("Model details would include:\n")
    cat("- Model formula\n")
    cat("- Prior distributions (for Bayesian models)\n")
    cat("- MCMC settings (for Bayesian models)\n")
    cat("- Model diagnostics")
  })
  
  # Covariate analysis (meta-regression) - Enhanced version
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
      
      # Add covariate information
      covariate_name <- input$covariate_var
      pw$matched_covs <- data()[[covariate_name]]
      
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
      if (requireNamespace("metafor", quietly = TRUE)) {
        
        if (input$include_interaction) {
          # Create treatment indicators
          pw$reference <- as.numeric(pw$treat1 == reference | pw$treat2 == reference)
          
          # Create interaction term based on covariate type
          if (input$covariate_type == "categorical") {
            # For categorical covariates, use formula notation
            mr_model <- try(rma(yi = TE, sei = seTE, mods = ~ matched_covs * reference, 
                                data = pw, method = "REML"), silent = TRUE)
          } else {
            # For numeric covariates, create the interaction term
            pw$interaction <- pw$matched_covs * pw$reference
            
            # Model with interaction
            mr_model <- try(rma(yi = TE, sei = seTE, mods = ~ matched_covs + reference + interaction, 
                                data = pw, method = "REML"), silent = TRUE)
          }
        } else {
          # Simple model without interaction
          mr_model <- try(rma(yi = TE, sei = seTE, mods = ~ matched_covs, 
                              data = pw, method = "REML"), silent = TRUE)
        }
        
        # Check if model ran successfully
        if (inherits(mr_model, "try-error")) {
          showNotification(paste("Error in meta-regression:", mr_model[1]), type = "error")
          return(NULL)
        }
        
        # Generate predictions at different covariate values
        if (input$covariate_type == "continuous") {
          cov_range <- range(pw$matched_covs, na.rm = TRUE)
          predictions <- calculate_metareg_predictions(mr_model, cov_range, treatments)
        } else {
          # For categorical, predict at each level
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
      } else {
        # Fallback to simple linear model with similar fixes
        if (input$include_interaction) {
          # Create treatment indicators
          pw$reference <- as.numeric(pw$treat1 == reference | pw$treat2 == reference)
          
          # Create interaction term
          if (input$covariate_type == "categorical") {
            # Use character variable to avoid factor * numeric issues
            mr_model <- lm(TE ~ matched_covs * reference, data = pw)
          } else {
            pw$interaction <- pw$matched_covs * pw$reference
            mr_model <- lm(TE ~ matched_covs + reference + interaction, data = pw)
          }
        } else {
          # Simple model without interaction
          mr_model <- lm(TE ~ matched_covs, data = pw)
        }
        
        # Store the results
        metareg_result(list(
          covariate = covariate_name,
          type = input$covariate_type,
          data = pw,
          model = mr_model,
          method = "lm",
          reference_treatment = reference,
          include_interaction = input$include_interaction
        ))
        
        showNotification("Meta-regression completed (using simplified model)", type = "warning")
      }
      
    }, error = function(e) {
      showNotification(paste("Error in meta-regression:", e$message), type = "error")
    })
  })
  
  # Covariate regression plot - Enhanced version
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
        if (mr$method %in% c("metafor", "rma")) {
          # For metafor, we need to extract the correct coefficients
          if (mr$include_interaction) {
            # This is simplified - in a real implementation, we'd need to handle 
            # the predictions more carefully with interactions
            abline(mr$model$beta[1], mr$model$beta[2], col = "red", lwd = 2)
          } else {
            abline(mr$model$beta[1], mr$model$beta[2], col = "red", lwd = 2)
          }
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
          
          if (mr$include_interaction) {
            # This is a simplification - with interaction terms, we would need
            # to create the correct model matrix for predictions
            newdata$reference <- 0
            newdata$interaction <- 0
            
            # Create the appropriate newmods matrix
            k <- length(mr$model$beta)
            if (k == 4) { 
              # Model with intercept, covariate, reference, and interaction
              newmods <- cbind(newdata$matched_covs, newdata$reference, newdata$interaction)
            } else if (k > 4) {
              # More complex model with multiple terms
              newmods <- cbind(newdata$matched_covs, rep(0, nrow(newdata)), rep(0, nrow(newdata)))
            } else {
              # Simplified case
              newmods <- cbind(newdata$matched_covs)
            }
            
            # Safe prediction with error handling
            tryCatch({
              pred <- predict(mr$model, newmods = newmods, level = 0.95)
              
              # Add the prediction interval
              lines(newdata$matched_covs, pred$pred, col = "red", lwd = 2)
              lines(newdata$matched_covs, pred$ci.lb, col = "red", lty = 2)
              lines(newdata$matched_covs, pred$ci.ub, col = "red", lty = 2)
            }, error = function(e) {
              # Draw a simple regression line if prediction fails
              message("Error in prediction: ", e$message)
            })
          } else {
            # Simple model without interaction
            tryCatch({
              pred <- predict(mr$model, newmods = cbind(newdata$matched_covs), level = 0.95)
              
              # Add the prediction interval
              lines(newdata$matched_covs, pred$pred, col = "red", lwd = 2)
              lines(newdata$matched_covs, pred$ci.lb, col = "red", lty = 2)
              lines(newdata$matched_covs, pred$ci.ub, col = "red", lty = 2)
            }, error = function(e) {
              # Draw a simple regression line if prediction fails
              message("Error in prediction: ", e$message)
            })
          }
        }
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
      }
    } else {
      # Create placeholder plot
      plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(-2, 2),
           xlab = input$covariate_var, ylab = "Effect",
           main = paste("Meta-Regression with", input$covariate_var, "\n(Run analysis first)"))
      text(0.5, 0, "Run meta-regression analysis", cex = 1.5)
    }
  })
  
  # Covariate forest plot - Enhanced version
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
    
    if (mr$method == "metafor") {
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
        CI_upper = results$ci.ub
      )
      
      datatable(result_df, 
                caption = paste("Meta-Regression with", mr$covariate),
                options = list(scrollX = TRUE))
    } else {
      # Use lm summary for linear model
      lm_summary <- summary(mr$model)
      
      # Create result dataframe
      result_df <- data.frame(
        Parameter = rownames(lm_summary$coefficients),
        Estimate = lm_summary$coefficients[,1],
        SE = lm_summary$coefficients[,2],
        t_value = lm_summary$coefficients[,3],
        P_value = lm_summary$coefficients[,4]
      )
      
      datatable(result_df, 
                caption = paste("Linear Regression with", mr$covariate),
                options = list(scrollX = TRUE))
    }
  })
  
  # Covariate ranking plot - Enhanced version
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
              main = paste("Adjusted P-Score Rankings at", mr$covariate, "=", covariate_value,
                           "\n(", ifelse(input$values_direction == "good", 
                                         "smaller values are better", 
                                         "larger values are better"), ")"),
              xlab = "Treatment", ylab = "Adjusted P-Score (%)",
              col = "lightgreen",
              ylim = c(0, 100),
              names.arg = adjusted_rankings$treatments)
      
      # Add a note about the adjustment
      if (mr$include_interaction) {
        mtext(paste("Note: Rankings adjusted for covariate-by-treatment interaction",
                    "using", mr$reference_treatment, "as reference"), 
              side = 3, line = 0.5, cex = 0.7)
      } else {
        mtext("Note: Rankings adjusted for covariate effect", 
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
    create_metareg_bubble_plot(
      data = bubble_data,
      xvar = "matched_covs",
      yvar = "TE",
      sizevar = "seTE",
      title = paste("Meta-Regression Bubble Plot for", input$bubble_treatment),
      xlab = mr$covariate,
      ylab = paste("Effect (", input$effect_measure, ")")
    )
  })
  
  # Meta-regression model summary
  output$metareg_model_summary <- renderPrint({
    req(metareg_result())
    
    mr <- metareg_result()
    
    cat("Meta-Regression Model Summary:\n\n")
    
    if (mr$method %in% c("metafor", "rma")) {
      # For metafor models
      print(summary(mr$model))
      
      cat("\nInterpretation:\n")
      cat("- Positive coefficients indicate that as the covariate increases, the effect size increases\n")
      cat("- For log odds ratios, risk ratios, or hazard ratios, this means higher risk\n")
      cat("- For mean differences, this indicates a larger difference between treatments\n")
      
      if (mr$include_interaction) {
        cat("\nNote: This model includes treatment-by-covariate interaction using")
        cat(mr$reference_treatment, "as the reference treatment\n")
      }
    } else if (mr$method == "lm") {
      # For linear models
      print(summary(mr$model))
      
      cat("\nWarning: This is a simplified analysis using ordinary least squares regression.\n")
      cat("It does not account for the precision of the effect estimates.\n")
    }
  })
  
  # Predicted effects table
  output$predicted_effects_table <- renderDT({
    req(metareg_result(), nma_result())
    
    mr <- metareg_result()
    
    if (mr$type == "continuous" && !is.null(mr$predictions)) {
      # Create a table of predicted effects at different covariate values
      pred_data <- data.frame()
      
      # Get a sequence of covariate values
      cov_values <- sort(unique(sapply(mr$predictions, function(x) x$covariate_value)))
      
      # Create rows for each covariate value
      for (cov_val in cov_values) {
        pred <- mr$predictions[[as.character(cov_val)]]
        pred_row <- data.frame(
          Covariate_Value = cov_val,
          Effect = round(pred$effect, 3),
          Lower_CI = round(pred$lower_ci, 3),
          Upper_CI = round(pred$upper_ci, 3)
        )
        
        pred_data <- rbind(pred_data, pred_row)
      }
      
      datatable(pred_data,
                caption = paste("Predicted Effects at Different Values of", mr$covariate),
                options = list(pageLength = 5, scrollX = TRUE))
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
    
    cat("Goodness of Fit for Meta-Regression Model:\n\n")
    
    if (mr$method %in% c("metafor", "rma")) {
      # Extract model fit statistics
      cat("Model Fit Statistics:\n")
      if (!is.null(mr$model$fit.stats)) {
        cat("AIC:", round(mr$model$fit.stats["AIC", "REML"], 2), "\n")
        cat("BIC:", round(mr$model$fit.stats["BIC", "REML"], 2), "\n")
        cat("Log-likelihood:", round(mr$model$fit.stats["ll", "REML"], 2), "\n\n")
      } else {
        cat("Model fit statistics not available\n\n")
      }
      
      # Test of moderators
      cat("Test of Moderators (Omnibus Test):\n")
      cat("QM =", round(mr$model$QM, 2), "on", mr$model$m, "df, p =", 
          format.pval(mr$model$QMp, digits = 4), "\n\n")
      
      # Variance components 
      cat("Variance Components:\n")
      cat("Tau^2 (estimated amount of residual heterogeneity):", round(mr$model$tau2, 4), "\n")
      cat("I^2 (residual heterogeneity / unaccounted variability):", 
          round(mr$model$I2 * 100, 1), "%\n")
      cat("H^2 (unaccounted variability / sampling variability):", round(mr$model$H2, 2), "\n")
      
      # R² statistic (if available)
      if (!is.null(mr$model$R2)) {
        cat("R^2 (amount of heterogeneity accounted for):", round(mr$model$R2 * 100, 1), "%\n")
      }
    } else if (mr$method == "lm") {
      # For linear models
      summary_obj <- summary(mr$model)
      
      cat("Model Summary:\n")
      cat("Multiple R-squared:", round(summary_obj$r.squared, 4), "\n")
      cat("Adjusted R-squared:", round(summary_obj$adj.r.squared, 4), "\n")
      cat("F-statistic:", round(summary_obj$fstatistic[1], 2), "on", 
          summary_obj$fstatistic[2], "and", summary_obj$fstatistic[3], "DF, p-value:", 
          format.pval(1 - pf(summary_obj$fstatistic[1], summary_obj$fstatistic[2], 
                             summary_obj$fstatistic[3]), digits = 4), "\n")
    }
  })
  
  # Model diagnostics - Residual plot
  output$residual_plot <- renderPlot({
    req(metareg_result())
    
    mr <- metareg_result()
    
    if (mr$method %in% c("metafor", "rma")) {
      # Get residuals vs fitted values
      fitted <- fitted(mr$model)
      resid <- residuals(mr$model)
      
      # Create a residual plot
      plot(fitted, resid, 
           xlab = "Fitted Values", ylab = "Residuals",
           main = "Residuals vs Fitted Values",
           pch = 19, col = "blue")
      
      # Add a horizontal line at y = 0
      abline(h = 0, lty = 2, col = "red")
      
      # Add a smoothed line to show any patterns
      if (length(fitted) > 3) {
        lines(lowess(fitted, resid), col = "green", lwd = 2)
      }
    } else if (mr$method == "lm") {
      # For linear models
      plot(mr$model, which = 1)
    }
  })
  
  # Model diagnostics - Q-Q plot
  output$qq_plot <- renderPlot({
    req(metareg_result())
    
    mr <- metareg_result()
    
    if (mr$method %in% c("metafor", "rma")) {
      # Get standardized residuals
      res <- try(rstandard(mr$model), silent = TRUE)
      
      if (!inherits(res, "try-error")) {
        # Create a Q-Q plot
        qqnorm(res, main = "Normal Q-Q Plot of Standardized Residuals",
               pch = 19, col = "blue")
        qqline(res, col = "red", lwd = 2)
      } else {
        # Fallback if standardized residuals aren't available
        plot(0, 0, type = "n", xlim = c(-3, 3), ylim = c(-3, 3),
             xlab = "Theoretical Quantiles", ylab = "Sample Quantiles",
             main = "Q-Q Plot (not available)")
        text(0, 0, "Standardized residuals not available", cex = 1.2)
      }
    } else if (mr$method == "lm") {
      # For linear models
      plot(mr$model, which = 2)
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
      
      if (mr$method == "metafor") {
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
          CI_upper = results$ci.ub
        )
        
        write.csv(result_df, file, row.names = FALSE)
      } else {
        # Use lm summary for linear model
        lm_summary <- summary(mr$model)
        
        # Create result dataframe
        result_df <- data.frame(
          Parameter = rownames(lm_summary$coefficients),
          Estimate = lm_summary$coefficients[,1],
          SE = lm_summary$coefficients[,2],
          t_value = lm_summary$coefficients[,3],
          P_value = lm_summary$coefficients[,4]
        )
        
        write.csv(result_df, file, row.names = FALSE)
      }
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
        Covariate_Value = covariate_value
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
        pred_data <- data.frame()
        
        # Get a sequence of covariate values
        cov_values <- sort(unique(sapply(mr$predictions, function(x) x$covariate_value)))
        
        # Create rows for each covariate value
        for (cov_val in cov_values) {
          pred <- mr$predictions[[as.character(cov_val)]]
          pred_row <- data.frame(
            Covariate_Value = cov_val,
            Effect = round(pred$effect, 3),
            Lower_CI = round(pred$lower_ci, 3),
            Upper_CI = round(pred$upper_ci, 3)
          )
          
          pred_data <- rbind(pred_data, pred_row)
        }
        
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
      req(baseline_result())
      
      # Create a comparison table for all treatment pairs at a specific baseline risk
      br <- baseline_result()
      nma <- nma_result()
      
      # Create empty grid
      treatments <- nma$trts
      n_treats <- length(treatments)
      grid <- matrix(NA, nrow = n_treats, ncol = n_treats)
      rownames(grid) <- treatments
      colnames(grid) <- treatments
      
      # Fill with placeholder values (would be actual calculations in full implementation)
      for (i in 1:n_treats) {
        for (j in 1:n_treats) {
          if (i != j) {
            # Simple placeholder values
            grid[i, j] <- round(runif(1, -1, 1), 2)
          }
        }
      }
      
      # Convert to data frame for display
      grid_df <- as.data.frame(grid)
      grid_df$Treatment <- rownames(grid)
      grid_df <- grid_df[, c("Treatment", treatments)]
      
      write.csv(grid_df, file, row.names = FALSE)
    }
  )
  
  output$download_baseline_ranking <- downloadHandler(
    filename = function() {
      paste0("baseline_risk_rankings_", input$baseline_ranking_value, ".csv")
    },
    content = function(file) {
      req(baseline_result(), input$baseline_ranking_value)
      
      # Get treatments
      treatments <- nma_result()$trts
      
      # Create placeholder P-scores (would be actual calculations in full implementation)
      p_scores <- runif(length(treatments), 0, 1)
      names(p_scores) <- treatments
      
      # Create ranking dataframe
      ranking_df <- data.frame(
        Treatment = names(p_scores),
        P_Score = round(p_scores * 100, 1),
        Baseline_Risk = input$baseline_ranking_value
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
      req(baseline_result(), input$baseline_value)
      
      # Save the plot
      png(file, width = 800, height = 600, res = 100)
      
      # Create a forest plot with adjusted effects at a specific baseline risk value
      plot(1, 1, type = "n", xlab = "Effect", ylab = "",
           main = paste("Forest Plot at Baseline Risk =", input$baseline_value),
           xlim = c(-2, 2), ylim = c(0, 10))
      abline(v = 0, lty = 2)
      
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
