# ---------------------------
# 1. Load Packages & Setup
# ---------------------------
# Ensure these packages are installed: install.packages(c("shiny", "bs4Dash", "shinycssloaders", "promises", "future", "netmeta", "visNetwork", "DT", "readr", "igraph", "dplyr"))
library(shiny)
library(bs4Dash)
library(shinycssloaders)
library(promises)
library(future)
library(netmeta)
library(visNetwork)
library(DT)
library(readr)
library(igraph)
library(dplyr) # Added for data manipulation like rename

plan(multisession) # Set up asynchronous processing

# ---------------------------
# 2. Data Upload Module (Revised & Corrected)
# ---------------------------
dataUploadUI <- function(id) {
  ns <- NS(id)
  tagList(
    fileInput(ns("file"), "Upload CSV Data", accept = ".csv"),
    radioButtons(ns("format"), "Data Format:",
                 choices = c("Arm-level (Binary or Survival/Ratio)" = "arm",
                             "Contrast-level (Log-Ratio)" = "contrast")),
    # Moved outcome type selection here
    selectInput(ns("outcomeType"), "Outcome Type:",
                choices = c("Binary (Events/N)" = "binary",
                            "Survival/Ratio (Log-HR/Log-Ratio)" = "ratio_cont"), # Clarified expected input
                selected = "binary"),
    # Conditional UI for summary measure (RR/OR) for binary data
    uiOutput(ns("summaryMeasureUI_upload")),
    checkboxInput(ns("logTrans"), "Input TE is on Ratio Scale (HR/OR/RR)? (Check to log-transform)", FALSE), # Clarified log transform meaning
    helpText("If using contrast-level data, TE should be log(HR), log(OR), or log(RR) unless the box above is checked."),
    helpText("For Survival/Ratio arm-level data, input log(HR) or log(Ratio) as 'mean' and SE(log-HR/log-Ratio) as 'sd'."),
    DTOutput(ns("dataPreview"))
  )
}

dataUploadServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    dataset <- reactiveVal(NULL)
    summaryMeasure <- reactiveVal("RR") # Default for binary
    
    # Dynamic UI for RR/OR selection
    output$summaryMeasureUI_upload <- renderUI({
      ns <- session$ns
      req(input$outcomeType)
      if (input$outcomeType == "binary") {
        radioButtons(ns("summaryMeasure"), "Binary Summary Measure:",
                     choices = c("Risk Ratio (RR)" = "RR",
                                 "Odds Ratio (OR)" = "OR"),
                     selected = summaryMeasure()) # Use reactiveVal default
      } else {
        NULL # No selection needed for survival/ratio
      }
    })
    
    # Update reactiveVal when radio button changes
    observeEvent(input$summaryMeasure, {
      # Only update if the input exists (i.e., outcomeType is binary)
      if(!is.null(input$summaryMeasure)){
        summaryMeasure(input$summaryMeasure)
      }
    })
    
    
    observeEvent(input$file, {
      req(input$file)
      raw <- tryCatch({
        read_csv(input$file$datapath, show_col_types = FALSE) %>%
          rename_all(tolower) # Convert all columns to lower case for consistency
      }, error = function(e) {
        showNotification(paste("Error reading CSV:", e$message), type = "error", duration=10)
        return(NULL)
      })
      
      if (is.null(raw)) return()
      
      processed <- NULL
      error_msg <- NULL
      sm_to_store <- NULL # Variable to store the determined summary measure
      
      if (input$format == "arm") {
        if (input$outcomeType == "binary") {
          req_cols <- c("study", "treatment", "event", "n")
          if (!all(req_cols %in% names(raw))) {
            error_msg <- paste("Binary arm-level data must contain:", paste(req_cols, collapse = ", "))
          } else {
            # Ensure columns are numeric
            raw$event <- as.numeric(raw$event)
            raw$n <- as.numeric(raw$n)
            if(any(is.na(raw$event)) || any(is.na(raw$n))) {
              error_msg <- "Non-numeric values found in 'event' or 'n' columns for binary arm-level data."
            } else if (any(raw$event > raw$n, na.rm=TRUE)) {
              error_msg <- "Error: 'event' count cannot be greater than 'n' (sample size)."
            } else {
              # Use the selected summary measure (RR or OR)
              sm_selected <- isolate(summaryMeasure()) # Get current selection
              sm_to_store <- sm_selected # Store for later use
              pw <- tryCatch({
                pairwise(treat = treatment, event = event, n = n,
                         studlab = study, data = raw, sm = sm_selected, allstudies = TRUE)
              }, error = function(e) {
                error_msg <<- paste("Error in pairwise calculation for binary data:", e$message)
                NULL
              })
              if (!is.null(pw)) {
                contrast <- as.data.frame(pw)
                # Ensure standard column names for netmeta
                contrast <- contrast %>% rename(te = TE, sete = seTE)
                processed <- list(raw = raw, contrast = contrast, format = "arm", sm = sm_to_store)
              }
            }
          }
        } else { # Survival/Ratio arm-level
          req_cols <- c("study", "treatment", "mean", "sd", "n")
          if (!all(req_cols %in% names(raw))) {
            error_msg <- paste("Survival/Ratio arm-level data must contain:", paste(req_cols, collapse = ", "))
          } else {
            # Ensure columns are numeric
            raw$mean <- as.numeric(raw$mean)
            raw$sd <- as.numeric(raw$sd)
            raw$n <- as.numeric(raw$n)
            if(any(is.na(raw$mean)) || any(is.na(raw$sd)) || any(is.na(raw$n))) {
              error_msg <- "Non-numeric values found in 'mean', 'sd', or 'n' columns for survival/ratio arm-level data."
            } else {
              # Use ROM, assuming user input log-HR as mean and SE(log-HR) as sd
              sm_to_store <- "HR" # Assume HR if using this pathway
              pw <- tryCatch({
                pairwise(treat = treatment, mean = mean, sd = sd, n = n,
                         studlab = study, data = raw, sm = "ROM", allstudies = TRUE)
              }, error = function(e) {
                error_msg <<- paste("Error in pairwise calculation for survival/ratio data:", e$message)
                NULL
              })
              if (!is.null(pw)) {
                contrast <- as.data.frame(pw)
                # Ensure standard column names for netmeta
                contrast <- contrast %>% rename(te = TE, sete = seTE)
                processed <- list(raw = raw, contrast = contrast, format = "arm", sm = sm_to_store)
              }
            }
          }
        }
      } else { # Contrast-level
        req_cols <- c("studlab", "treat1", "treat2", "te", "sete") # Use lower case
        if (!all(req_cols %in% names(raw))) {
          error_msg <- paste("Contrast-level data must contain:", paste(req_cols, collapse = ", "))
        } else {
          # Ensure columns are numeric
          raw$te <- as.numeric(raw$te)
          raw$sete <- as.numeric(raw$sete)
          if(any(is.na(raw$te)) || any(is.na(raw$sete))) {
            error_msg <- "Non-numeric values found in 'te' or 'sete' columns for contrast-level data."
          } else {
            # Determine sm based on outcome type selection
            sm_to_store <- if(isolate(input$outcomeType) == "binary") {
              isolate(summaryMeasure())
            } else {
              "HR" # Assume HR for ratio_cont
            }
            processed <- list(raw = raw, contrast = raw, format = "contrast", sm = sm_to_store)
          }
        }
      }
      
      if (!is.null(error_msg)) {
        showNotification(error_msg, type = "error", duration = 10)
        dataset(NULL)
        output$dataPreview <- renderDT(NULL)
        return()
      }
      
      # Log transformation for contrast data if requested and TE is positive
      if (input$format == "contrast" && input$logTrans && "te" %in% names(processed$contrast)) {
        if(all(processed$contrast$te > 0, na.rm = TRUE)) {
          processed$contrast$te <- log(processed$contrast$te)
          # Warning about SE transformation
          if(any(processed$contrast$sete / abs(processed$contrast$te) > 0.5, na.rm = TRUE)) {
            showNotification("Warning: Log-transforming TE. Standard errors (seTE) were assumed to be for the log-transformed scale. If they were for the original ratio scale, results may be inaccurate.", type="warning", duration=15)
          }
        } else {
          showNotification("Cannot log-transform TE because some values are zero or negative.", type="error", duration=10)
          dataset(NULL) # Prevent further processing
          output$dataPreview <- renderDT(NULL)
          return()
        }
      }
      
      
      dataset(processed)
      output$dataPreview <- renderDT({
        # Ensure column names are consistent before display
        disp_data <- processed$contrast
        # Use standard names for display if they exist
        if("te" %in% names(disp_data)) disp_data <- rename(disp_data, TE = te)
        if("sete" %in% names(disp_data)) disp_data <- rename(disp_data, seTE = sete)
        datatable(head(disp_data, 10), options = list(scrollX = TRUE))
      })
    })
    
    # Return the reactive list containing raw data, contrast data, format, and sm
    return(reactive({ dataset() }))
  })
}


# ---------------------------
# 3. Network Plot Module (No changes needed from previous version)
# ---------------------------
networkUI <- function(id) {
  ns <- NS(id)
  tagList(
    selectInput(ns("focus"), "Focus on Treatment:", choices = NULL),
    visNetworkOutput(ns("networkPlot")) %>% withSpinner(),
    br(),
    h4("Customize Static Network Plot (igraph)"),
    fluidRow(
      column(4,
             selectInput(ns("layoutOption"), "Layout:",
                         choices = c("Fruchterman-Reingold" = "layout_with_fr",
                                     "Kamada-Kawai" = "layout_with_kk",
                                     "Reingold-Tilford" = "layout_as_tree",
                                     "Random" = "layout_randomly"),
                         selected = "layout_with_fr")
      ),
      column(4,
             numericInput(ns("vertexSize"), "Vertex Size:", value = 30, min = 10, max = 100)
      ),
      column(4,
             selectInput(ns("vertexColor"), "Vertex Color:",
                         choices = c("skyblue", "red", "green", "blue", "orange", "purple", "gray", "black", "yellow", "pink"),
                         selected = "skyblue")
      )
    ),
    fluidRow(
      column(4,
             numericInput(ns("vertexLabelCex"), "Vertex Label Size:", value = 1.2, min = 0.5, max = 3)
      ),
      column(4,
             selectInput(ns("edgeColor"), "Edge Color:",
                         choices = c("gray", "black", "red", "blue", "green", "orange", "purple"),
                         selected = "gray")
      ),
      column(4,
             numericInput(ns("edgeWidth"), "Edge Width:", value = 1, min = 0.5, max = 5)
      )
    ),
    br(),
    h4("Static Network Plot Preview"),
    plotOutput(ns("staticNetworkPlot")) %>% withSpinner(),
    br(),
    downloadButton(ns("dlNetworkIgraphPNG"), "Download Static Network Plot")
  )
}

networkServer <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    observe({ # Use observe to react whenever data() changes
      d <- data() # Get the reactive value
      req(d) # Ensure data is not NULL
      # Use consistent lowercase column names
      if (!is.null(d$contrast) && all(c("treat1", "treat2") %in% names(d$contrast))) {
        updateSelectInput(session, "focus", choices = unique(c(d$contrast$treat1, d$contrast$treat2)))
      }
    })
    
    networkData <- reactive({
      d <- data() # Get the reactive value
      req(d) # Ensure data is not NULL
      # Use consistent lowercase column names
      if (!is.null(d$contrast) && all(c("treat1", "treat2") %in% names(d$contrast))) {
        nodes <- data.frame(id = unique(c(d$contrast$treat1, d$contrast$treat2)),
                            label = unique(c(d$contrast$treat1, d$contrast$treat2)))
        edges <- data.frame(from = d$contrast$treat1, to = d$contrast$treat2)
        if (!is.null(input$focus) && input$focus != "") {
          edges <- edges[edges$from == input$focus | edges$to == input$focus, ]
        }
        list(nodes = nodes, edges = edges)
      } else {
        NULL # Return NULL if data is not ready
      }
    })
    
    output$networkPlot <- renderVisNetwork({
      net <- networkData()
      validate(need(!is.null(net), "No network data available or data format incorrect."))
      visNetwork(net$nodes, net$edges, height = "500px", width = "100%") %>%
        visNodes(shape = "dot", size = 20) %>%
        visEdges(smooth = FALSE)
    })
    
    output$staticNetworkPlot <- renderPlot({
      d <- data() # Get the reactive value
      req(d) # Ensure data is not NULL
      # Use consistent lowercase column names
      if (!is.null(d$contrast) && all(c("treat1", "treat2") %in% names(d$contrast))) {
        edges <- d$contrast[, c("treat1", "treat2")]
        g <- igraph::graph_from_data_frame(edges, directed = FALSE)
        layoutFunc <- switch(input$layoutOption,
                             layout_with_fr = igraph::layout_with_fr,
                             layout_with_kk = igraph::layout_with_kk,
                             layout_as_tree = igraph::layout_as_tree,
                             layout_randomly = igraph::layout_randomly)
        lay <- layoutFunc(g)
        plot(g, layout = lay,
             vertex.size = input$vertexSize,
             vertex.color = input$vertexColor,
             vertex.label.cex = input$vertexLabelCex,
             edge.color = input$edgeColor,
             edge.width = input$edgeWidth,
             main = "Static Network Plot (igraph)")
      } else {
        plot.new(); text(0.5, 0.5, "Contrast data with 'treat1' and 'treat2' columns needed.")
      }
    })
    
    output$dlNetworkIgraphPNG <- downloadHandler(
      filename = function() { "NetworkPlot_Static.png" },
      content = function(file) {
        png(file, width = 800, height = 600, units = "px", res = 96)
        d <- data() # Get the reactive value
        if (!is.null(d) && !is.null(d$contrast) && all(c("treat1", "treat2") %in% names(d$contrast))) {
          edges <- d$contrast[, c("treat1", "treat2")]
          g <- igraph::graph_from_data_frame(edges, directed = FALSE)
          layoutFunc <- switch(input$layoutOption,
                               layout_with_fr = igraph::layout_with_fr,
                               layout_with_kk = igraph::layout_with_kk,
                               layout_as_tree = igraph::layout_as_tree,
                               layout_randomly = igraph::layout_randomly)
          lay <- layoutFunc(g)
          plot(g, layout = lay,
               vertex.size = input$vertexSize,
               vertex.color = input$vertexColor,
               vertex.label.cex = input$vertexLabelCex,
               edge.color = input$edgeColor,
               edge.width = input$edgeWidth,
               main = "Static Network Plot (igraph)")
          dev.off()
        } else {
          # Create a blank plot with error message if download is attempted without data
          plot.new()
          text(0.5, 0.5, "No network data available for download.")
          dev.off()
          showNotification("No network data available for download.", type = "error")
        }
      }
    )
  })
}


# ---------------------------
# 4. Analysis Module (Corrected - accepts selected_sm)
# ---------------------------
analysisUI <- function(id) {
  ns <- NS(id)
  tagList(
    # Removed method selection as only frequentist is used
    radioButtons(ns("effectModel"), "Effect Model:",
                 choices = c("Fixed-effect" = "fixed", "Random-effects" = "random"),
                 selected = "random"),
    sliderInput(ns("confLevel"), "Confidence Level:",
                min = 0.80, max = 0.99, value = 0.95, step = 0.01),
    verbatimTextOutput(ns("status"))
  )
}

analysisServer <- function(id, data, selected_sm) { # Added selected_sm argument
  moduleServer(id, function(input, output, session) {
    resultContainer <- reactiveVal(NULL)
    
    observe({ # Use observe instead of observeEvent to react to data OR model changes
      d_all <- data() # Get the reactive list
      req(d_all, d_all$contrast) # Ensure data and contrast data are available
      # If sm comes from the data list (arm-level processing), use it.
      # Otherwise, get it from the reactive selected_sm (for contrast-level).
      sm_from_data <- d_all$sm
      sm_from_input <- selected_sm()
      
      localSM <- if (!is.null(sm_from_data)) {
        sm_from_data
      } else if (!is.null(sm_from_input)) {
        sm_from_input
      } else {
        "HR" # Fallback default if neither is available
      }
      
      req(localSM) # Ensure we have a summary measure
      
      # Isolate reactive inputs before the future call
      localEffectModel <- isolate(input$effectModel)
      localConfLevel   <- isolate(input$confLevel)
      
      output$status <- renderText("Running frequentist analysis...")
      
      future({
        # Ensure correct column names (lowercase)
        contrast_data <- d_all$contrast %>%
          rename_with(tolower)
        
        # Make sure TE and seTE columns exist after renaming
        required_contrast_cols <- c("studlab", "treat1", "treat2", "te", "sete")
        if (!all(required_contrast_cols %in% names(contrast_data))) {
          stop(paste("Required contrast columns missing. Need:", paste(required_contrast_cols, collapse=", ")))
        }
        
        netfit <- tryCatch({
          netmeta(
            TE = te, seTE = sete, treat1 = treat1, treat2 = treat2, studlab = studlab,
            data = contrast_data,
            common = (localEffectModel == "fixed"),
            random = (localEffectModel == "random"),
            level = localConfLevel,
            sm = localSM, # Pass the determined summary measure
            reference.group = "", # Let netmeta determine reference
            details.chkmultiarm = TRUE,
            tol.chkmultiarm = 1e-4 # Added tolerance
          )
        }, error = function(e){
          print(paste("Error in netmeta:", e$message)) # Print error to console
          # Try to provide more info if possible
          if(grepl("Treatment names must be identical", e$message)){
            print("Check consistency of treatment names across treat1 and treat2 columns.")
          }
          if(grepl("Need at least three different treatments", e$message)){
            print("Network meta-analysis requires at least 3 distinct treatments in the data.")
          }
          NULL # Return NULL on error
        })
        
        if(is.null(netfit)){
          return(list(method = "freq", model = NULL, nma_text = "Error during netmeta execution. Check console/debug output.", sm = localSM))
        }
        
        # Capture output safely
        nma_text <- capture.output(print(netfit))
        list(method = "freq", model = netfit, nma_text = paste(nma_text, collapse = "\n"), sm = localSM) # Store sm in result
      }) %...>% resultContainer
    }) # End observe
    
    return(reactive({ resultContainer() }))
  }) # End moduleServer
}


# ---------------------------
# 5. Diagnostics Module (No changes needed from previous version)
# ---------------------------
diagnosticsUI <- function(id) {
  ns <- NS(id)
  tagList(
    actionButton(ns("runDiagnostics"), "Run Diagnostics"),
    verbatimTextOutput(ns("diagText")),
    plotOutput(ns("funnelPlot")) %>% withSpinner(),
    br(),
    downloadButton(ns("dlFunnelPNG"), "Download Funnel Plot as PNG")
  )
}

diagnosticsServer <- function(id, analysisResult) {
  moduleServer(id, function(input, output, session) {
    observeEvent(input$runDiagnostics, {
      req(analysisResult())
      res <- analysisResult()
      if (is.null(res$model) || res$method != "freq") {
        output$diagText <- renderText("Diagnostics require a successful frequentist analysis.")
      } else {
        diag_res <- tryCatch(decomp.design(res$model), error = function(e) e)
        if(inherits(diag_res, "error")){
          output$diagText <- renderPrint({ paste("Error running decomp.design:", diag_res$message) })
        } else {
          output$diagText <- renderPrint({ diag_res })
        }
      }
    })
    output$funnelPlot <- renderPlot({
      req(analysisResult())
      res <- analysisResult()
      if (!is.null(res$model) && res$method == "freq") {
        # Check if treatments exist and there are at least 2 before plotting
        if(!is.null(res$model$trts) && length(res$model$trts) >= 2){
          ref <- res$model$trts[1] # Use first treatment as reference
          tryCatch(funnel(res$model, order = ref), error = function(e){
            plot.new(); text(0.5, 0.5, paste("Funnel plot error:", e$message))
          })
        } else {
          plot.new(); text(0.5, 0.5, "Funnel plot requires at least 2 treatments.")
        }
      } else {
        plot.new()
        text(0.5, 0.5, "Funnel plot not available.", cex = 1.3)
      }
    })
    output$dlFunnelPNG <- downloadHandler(
      filename = function() { "FunnelPlot.png" },
      content = function(file) {
        png(file, width = 800, height = 600, units = "px", res = 96)
        res <- analysisResult()
        if (!is.null(res$model) && res$method == "freq") {
          if(!is.null(res$model$trts) && length(res$model$trts) >= 2){
            ref <- res$model$trts[1]
            tryCatch(funnel(res$model, order = ref), error = function(e){
              plot.new(); text(0.5, 0.5, paste("Funnel plot error:", e$message))
            })
          } else {
            plot.new(); text(0.5, 0.5, "Funnel plot requires at least 2 treatments.")
          }
        } else {
          plot.new()
          text(0.5, 0.5, "Funnel plot not available.", cex = 1.3)
        }
        dev.off()
      }
    )
  })
}


# ---------------------------
# 6. Results Module (Revised for dynamic labels and checks)
# ---------------------------
resultsUI <- function(id) {
  ns <- NS(id)
  tabsetPanel(
    tabPanel("Forest Plot",
             fluidRow(
               column(3, textInput(ns("forestTitle"), "Plot Title:", value = "Forest Plot")),
               column(3, uiOutput(ns("forestXlabUI"))), # Dynamic X-label
               column(3, sliderInput(ns("forestXlim"), "X-axis Range:",
                                     min = 0.1, max = 10, value = c(0.5, 3), step = 0.1)),
               column(3, numericInput(ns("forestCex"), "Label Size (cex):", value = 1, min = 0.5, max = 3, step = 0.1))
             ),
             fluidRow(
               column(3, numericInput(ns("forestDigits"), "Digits (Ratio):", value = 2, min = 0, max = 5, step = 1)),
               column(3, numericInput(ns("forestDigitsSE"), "Digits (SE):", value = 2, min = 0, max = 5, step = 1)),
               column(3, numericInput(ns("forestDigitsPval"), "Digits (p-val):", value = 3, min = 0, max = 5, step = 1))
             ),
             fluidRow(
               column(3, selectInput(ns("forestSquareColor"), "Square Color:",
                                     choices = c("black", "blue", "red", "green", "purple"), selected = "black")),
               column(3, selectInput(ns("forestDiamondColor"), "Diamond Color:",
                                     choices = c("black", "blue", "red", "green", "purple"), selected = "blue"))
             ),
             br(),
             plotOutput(ns("forestPlot")) %>% withSpinner(),
             br(),
             downloadButton(ns("dlForestPNG"), "Download Forest Plot as PNG")
    ),
    tabPanel("Relative Effects",
             verbatimTextOutput(ns("relativeEffectsText")),
             br(),
             downloadButton(ns("dlRelativeEffects"), "Download Relative Effects Summary")
    ),
    tabPanel("NMA Results",
             verbatimTextOutput(ns("nmaResults"))
    ),
    tabPanel("Net Heat Plot",
             plotOutput(ns("netHeatPlot")) %>% withSpinner(),
             br(),
             downloadButton(ns("dlNetHeatPNG"), "Download Net Heat Plot as PNG")
    ),
    tabPanel("Net Splitting",
             plotOutput(ns("netSplitPlot")) %>% withSpinner(),
             br(),
             downloadButton(ns("dlNetSplitPNG"), "Download Net Split Plot as PNG")
    ),
    tabPanel("Comparison-Adjusted Funnel Plot",
             plotOutput(ns("funnelPlotAdj")) %>% withSpinner(),
             br(),
             downloadButton(ns("dlFunnelAdjPNG"), "Download Funnel Plot as PNG")
    ),
    tabPanel("Leave-One-Out Analysis",
             DTOutput(ns("leaveOneOut")),
             br(),
             downloadButton(ns("dlLeaveOneOut"), "Download Leave-One-Out Table as CSV")
    ),
    tabPanel("Detailed Model Summary",
             verbatimTextOutput(ns("detailedModelSummary")),
             br(),
             downloadButton(ns("dlDetailedModelSummary"), "Download Model Summary as TXT")
    ),
    tabPanel("Network Graph",
             fluidRow(
               column(6, textInput(ns("netGraphTitle"), "Graph Title:", value = "Network Graph")),
               column(6, textInput(ns("longLabels"), "Full Treatment Labels (comma separated):", value = "DrugA,DrugB,DrugC,DrugD"))
             ),
             plotOutput(ns("netGraphPlot")) %>% withSpinner(),
             br(),
             downloadButton(ns("dlNetGraphPNG"), "Download Network Graph as PNG")
    ),
    tabPanel("Effect Estimate Table",
             DTOutput(ns("effectTable")),
             br(),
             downloadButton(ns("dlEffectTable"), "Download Effect Table as CSV")
    ),
    tabPanel("Treatment Ranking",
             DTOutput(ns("treatmentRanking")),
             br(),
             downloadButton(ns("dlTreatmentRanking"), "Download Treatment Ranking as CSV")
    )
  )
}

resultsServer <- function(id, analysisResult) {
  moduleServer(id, function(input, output, session) {
    
    # Dynamic UI for Forest Plot X-axis Label
    output$forestXlabUI <- renderUI({
      ns <- session$ns
      res <- analysisResult()
      default_lab <- "Effect Measure (Log Scale)" # Fallback
      if(!is.null(res) && !is.null(res$sm)){ # Check analysis result and sm
        if(res$sm == "RR") default_lab <- "Risk Ratio (RR)"
        if(res$sm == "OR") default_lab <- "Odds Ratio (OR)"
        if(res$sm %in% c("HR", "ROM")) default_lab <- "Hazard Ratio / Ratio of Means"
      }
      textInput(ns("forestXlab"), "X-axis Label:", value = default_lab)
    })
    
    # Forest Plot
    output$forestPlot <- renderPlot({
      res <- analysisResult()
      validate(need(!is.null(res) && !is.null(res$model), "NMA model not available. Please run analysis first."))
      if (res$method == "freq") {
        # Use tryCatch for plotting robustness
        tryCatch({
          # Ensure xlim is valid
          xlim_val <- input$forestXlim
          if(is.null(xlim_val) || length(xlim_val) != 2 || xlim_val[1] >= xlim_val[2]){
            xlim_val <- NULL # Use default if invalid
          }
          
          netmeta:::forest.netmeta(
            res$model,
            main = input$forestTitle,
            xlab = input$forestXlab,
            xlim = xlim_val,
            cex = input$forestCex,
            digits = input$forestDigits,
            digits.se = input$forestDigitsSE,
            digits.pval = input$forestDigitsPval,
            col.square = input$forestSquareColor,
            col.diamond = input$forestDiamondColor,
            transf = exp # Exponentiate log-ratios
          )
        }, error = function(e){
          plot.new(); text(0.5, 0.5, paste("Forest plot error:", e$message))
        })
      } else {
        plot.new()
        text(0.5, 0.5, "Forest plot not available for this method.", cex = 1.3)
      }
    })
    
    output$dlForestPNG <- downloadHandler(
      filename = function() { "NMA_ForestPlot.png" },
      content = function(file) {
        png(file, width = 800, height = 600, units = "px", res = 96)
        res <- analysisResult()
        if (!is.null(res$model) && res$method == "freq") {
          tryCatch({
            xlim_val <- input$forestXlim
            if(is.null(xlim_val) || length(xlim_val) != 2 || xlim_val[1] >= xlim_val[2]){
              xlim_val <- NULL
            }
            netmeta:::forest.netmeta(
              res$model,
              main = input$forestTitle,
              xlab = input$forestXlab,
              xlim = xlim_val,
              cex = input$forestCex,
              digits = input$forestDigits,
              digits.se = input$forestDigitsSE,
              digits.pval = input$forestDigitsPval,
              col.square = input$forestSquareColor,
              col.diamond = input$forestDiamondColor,
              transf = exp
            )
          }, error = function(e){
            plot.new(); text(0.5, 0.5, paste("Forest plot error:", e$message))
          })
        } else {
          plot.new()
          text(0.5, 0.5, "Forest plot not available.", cex = 1.3)
        }
        dev.off()
      }
    )
    
    # Relative Effects Table
    output$relativeEffectsText <- renderPrint({
      res <- analysisResult()
      validate(need(!is.null(res) && !is.null(res$model), "NMA model not available."))
      if (res$method == "freq" && !is.null(res$model$trts) &&
          length(unique(res$model$trts)) >= 2) {
        league <- tryCatch(netleague(res$model), error=function(e) NULL)
        if(is.null(league)) return("Could not generate league table.")
        print(league)
      } else {
        cat("Relative effects summary not available (insufficient treatments or non-frequentist model).")
      }
    })
    
    output$dlRelativeEffects <- downloadHandler(
      filename = function() { "RelativeEffectsSummary.txt" },
      content = function(file) {
        res <- analysisResult()
        if (!is.null(res$model) && res$method == "freq" && !is.null(res$model$trts) &&
            length(unique(res$model$trts)) >= 2) {
          league <- tryCatch(netleague(res$model), error=function(e) NULL)
          if(is.null(league)){
            writeLines("Could not generate league table.", con = file)
          } else {
            txt <- capture.output(print(league))
            writeLines(txt, con = file)
          }
        } else {
          writeLines("Relative effects summary not available.", con = file)
        }
      }
    )
    
    # NMA Results Summary
    output$nmaResults <- renderPrint({
      res <- analysisResult()
      validate(need(!is.null(res) && !is.null(res$model), "NMA model not available."))
      cat("=== NMA Analysis Summary ===\n")
      # Print captured summary text if available, otherwise print model object
      if(!is.null(res$nma_text) && nzchar(res$nma_text)){
        cat(res$nma_text, "\n")
      } else {
        print(res$model) # Fallback
      }
      # Add key stats explicitly if available in model object
      if (!is.null(res$model$k)) cat("Number of studies:", res$model$k, "\n")
      if (!is.null(res$model$n)) cat("Number of observations:", res$model$n, "\n") # If available
      if (!is.null(res$model$Q)) cat("Q statistic: ", round(res$model$Q,2), "\n")
      if (!is.null(res$model$I2)) cat("I-squared: ", scales::percent(res$model$I2, accuracy=0.1), "\n")
      if (!is.null(res$model$tau)) cat("Tau: ", round(res$model$tau,3), "\n")
      if (!is.null(res$model$tau2)) cat("Tau^2: ", round(res$model$tau2,3), "\n") # Use tau2 if available
      if (!is.null(res$model$pval.Q)) cat("p-value (heterogeneity): ", format.pval(res$model$pval.Q, digits=3), "\n")
      if (!is.null(res$model$studies.excluded) && length(res$model$studies.excluded) > 0) {
        cat("Excluded studies in this analysis:", paste(res$model$studies.excluded, collapse = ", "), "\n")
      }
    })
    
    # --- Net Heat Plot ---
    output$netHeatPlot <- renderPlot({
      res <- analysisResult()
      validate(need(!is.null(res) && !is.null(res$model), "NMA model not available."))
      tryCatch(netheat(res$model, random = (res$model$random)), error = function(e){
        plot.new(); text(0.5, 0.5, paste("Net Heat plot error:", e$message))
      })
    })
    output$dlNetHeatPNG <- downloadHandler(
      filename = function() { "NetHeatPlot.png" },
      content = function(file) {
        png(file, width = 800, height = 600, units = "px", res = 96)
        res <- analysisResult()
        if(!is.null(res$model)){
          tryCatch(netheat(res$model, random = (res$model$random)), error = function(e){
            plot.new(); text(0.5, 0.5, paste("Net Heat plot error:", e$message))
          })
        } else { plot.new(); text(0.5, 0.5, "Model not available.") }
        dev.off()
      }
    )
    
    # --- Net Splitting ---
    output$netSplitPlot <- renderPlot({
      res <- analysisResult()
      validate(need(!is.null(res) && !is.null(res$model), "NMA model not available."))
      nsplit <- tryCatch(netsplit(res$model), error = function(e) NULL)
      validate(need(!is.null(nsplit), "Net splitting analysis failed or not applicable (e.g., insufficient paths)."))
      tryCatch(plot(nsplit, main = "Net Splitting Analysis"), error = function(e){
        plot.new(); text(0.5, 0.5, paste("Net Split plot error:", e$message))
      })
    })
    output$dlNetSplitPNG <- downloadHandler(
      filename = function() { "NetSplitPlot.png" },
      content = function(file) {
        png(file, width = 800, height = 600, units = "px", res = 96)
        res <- analysisResult()
        if(!is.null(res$model)){
          nsplit <- tryCatch(netsplit(res$model), error = function(e) NULL)
          if(!is.null(nsplit)){
            tryCatch(plot(nsplit, main = "Net Splitting Analysis"), error = function(e){
              plot.new(); text(0.5, 0.5, paste("Net Split plot error:", e$message))
            })
          } else { plot.new(); text(0.5, 0.5, "Net splitting failed or not applicable.") }
        } else { plot.new(); text(0.5, 0.5, "Model not available.") }
        dev.off()
      }
    )
    
    # --- Comparison-Adjusted Funnel Plot ---
    output$funnelPlotAdj <- renderPlot({
      res <- analysisResult()
      validate(need(!is.null(res) && !is.null(res$model), "NMA model not available."))
      if (!is.null(res$model$trts) && length(res$model$trts) >= 2) {
        tryCatch(funnel(res$model, order = res$model$trts, pch = 19, method.bias = "Egger"), error = function(e){
          plot.new(); text(0.5, 0.5, paste("Funnel plot error:", e$message))
        })
      } else {
        plot.new()
        text(0.5, 0.5, "Comparison-adjusted funnel plot requires at least 2 treatments.", cex = 1.3)
      }
    })
    output$dlFunnelAdjPNG <- downloadHandler(
      filename = function() { "ComparisonAdjustedFunnelPlot.png" },
      content = function(file) {
        png(file, width = 800, height = 600, units = "px", res = 96)
        res <- analysisResult()
        if (!is.null(res$model) && !is.null(res$model$trts) && length(res$model$trts) >= 2) {
          tryCatch(funnel(res$model, order = res$model$trts, pch = 19, method.bias = "Egger"), error = function(e){
            plot.new(); text(0.5, 0.5, paste("Funnel plot error:", e$message))
          })
        } else {
          plot.new()
          text(0.5, 0.5, "Funnel plot requires at least 2 treatments.", cex = 1.3)
        }
        dev.off()
      }
    )
    
    # --- Leave-One-Out Analysis ---
    output$leaveOneOut <- renderDT({
      res <- analysisResult()
      validate(need(!is.null(res) && !is.null(res$model), "NMA model not available."))
      loo_res <- tryCatch({
        # Use the internal function if available, otherwise return message
        if(exists("leave1out", envir = asNamespace("netmeta"))) {
          netmeta:::leave1out(res$model)
        } else {
          data.frame(Message = "leave1out function not currently available in loaded netmeta version.")
        }
      }, error = function(e) {
        data.frame(Error = paste("Leave-one-out analysis failed:", e$message))
      })
      datatable(loo_res, options = list(pageLength = 10, scrollX = TRUE))
    })
    
    output$dlLeaveOneOut <- downloadHandler(
      filename = function() { "LeaveOneOutAnalysis.csv" },
      content = function(file) {
        res <- analysisResult()
        validate(need(!is.null(res) && !is.null(res$model), "NMA model not available."))
        loo_res <- tryCatch({
          if(exists("leave1out", envir = asNamespace("netmeta"))) {
            netmeta:::leave1out(res$model)
          } else {
            data.frame(Message = "leave1out function not currently available in loaded netmeta version.")
          }
        }, error = function(e) {
          data.frame(Error = paste("Leave-one-out analysis failed:", e$message))
        })
        write.csv(loo_res, file, row.names = FALSE)
      }
    )
    
    # --- Detailed Model Summary ---
    output$detailedModelSummary <- renderPrint({
      res <- analysisResult()
      validate(need(!is.null(res) && !is.null(res$model), "NMA model not available."))
      summary(res$model)
    })
    output$dlDetailedModelSummary <- downloadHandler(
      filename = function() { "DetailedModelSummary.txt" },
      content = function(file) {
        res <- analysisResult()
        validate(need(!is.null(res) && !is.null(res$model), "NMA model not available."))
        txt <- capture.output(summary(res$model))
        writeLines(txt, con = file)
      }
    )
    
    # --- Network Graph ---
    output$netGraphPlot <- renderPlot({
      res <- analysisResult()
      validate(need(!is.null(res) && !is.null(res$model), "NMA model not available."))
      long.labels <- trimws(unlist(strsplit(input$longLabels, split = ",")))
      # Use model treatments if labels don't match or are empty
      if (length(long.labels) == 0 || length(long.labels) != length(res$model$trts)) {
        long.labels <- res$model$trts
        updateTextInput(session, "longLabels", value = paste(long.labels, collapse=",")) # Update UI
      }
      tryCatch(netgraph(res$model, labels = long.labels, main = input$netGraphTitle), error = function(e){
        plot.new(); text(0.5, 0.5, paste("Network graph error:", e$message))
      })
    })
    output$dlNetGraphPNG <- downloadHandler(
      filename = function() { "NetworkGraph.png" },
      content = function(file) {
        png(file, width = 800, height = 600, units = "px", res = 96)
        res <- analysisResult()
        if(!is.null(res$model)){
          long.labels <- trimws(unlist(strsplit(input$longLabels, split = ",")))
          if (length(long.labels) == 0 || length(long.labels) != length(res$model$trts)) long.labels <- res$model$trts
          tryCatch(netgraph(res$model, labels = long.labels, main = input$netGraphTitle), error = function(e){
            plot.new(); text(0.5, 0.5, paste("Network graph error:", e$message))
          })
        } else { plot.new(); text(0.5, 0.5, "Model not available.") }
        dev.off()
      }
    )
    
    # --- Effect Estimate Table (League Table) ---
    output$effectTable <- renderDT({
      res <- analysisResult()
      validate(need(!is.null(res) && !is.null(res$model), "NMA model not available."))
      league <- tryCatch(netleague(res$model), error=function(e) NULL)
      validate(need(!is.null(league), "Could not generate league table."))
      # Decide whether to show fixed or random based on model run
      table_to_show <- if(!is.null(league$random)) league$random else league$fixed
      # Ensure table is a data frame for datatable
      if(!is.data.frame(table_to_show)) table_to_show <- as.data.frame(table_to_show)
      datatable(table_to_show, options = list(pageLength = 10, scrollX = TRUE))
    })
    output$dlEffectTable <- downloadHandler(
      filename = function() { "EffectEstimateTable.csv" },
      content = function(file) {
        res <- analysisResult()
        validate(need(!is.null(res) && !is.null(res$model), "NMA model not available."))
        league <- tryCatch(netleague(res$model), error=function(e) NULL)
        validate(need(!is.null(league), "Could not generate league table."))
        table_to_show <- if(!is.null(league$random)) league$random else league$fixed
        write.csv(table_to_show, file, row.names = TRUE) # Keep row names for league table
      }
    )
    
    # --- Treatment Ranking ---
    output$treatmentRanking <- renderDT({
      res <- analysisResult()
      validate(need(!is.null(res) && !is.null(res$model), "NMA model not available."))
      ranking <- tryCatch(netrank(res$model, small.values = "good"), error=function(e) NULL)
      validate(need(!is.null(ranking), "Could not calculate treatment rankings."))
      # Use common or random depending on model
      rank_scores <- if(!is.null(ranking$ranking.random)) ranking$ranking.random else ranking$ranking.common
      p_scores <- if(!is.null(ranking$Pscore.random)) ranking$Pscore.random else ranking$Pscore.common
      
      # Check if scores are NULL or empty before creating data frame
      validate(need(!is.null(rank_scores) && !is.null(p_scores) && length(rank_scores) > 0, "Ranking scores are unavailable."))
      
      ranking_df <- data.frame(
        Treatment = names(rank_scores),
        PScore = as.numeric(p_scores)
      ) %>% arrange(desc(PScore)) # Arrange by P-score
      
      datatable(ranking_df, options = list(pageLength = 10, order = list(list(1, 'desc')))) # Order by PScore column (index 1)
    })
    output$dlTreatmentRanking <- downloadHandler(
      filename = function() { "TreatmentRanking.csv" },
      content = function(file) {
        res <- analysisResult()
        validate(need(!is.null(res) && !is.null(res$model), "NMA model not available."))
        ranking <- tryCatch(netrank(res$model, small.values = "good"), error=function(e) NULL)
        validate(need(!is.null(ranking), "Could not calculate treatment rankings."))
        rank_scores <- if(!is.null(ranking$ranking.random)) ranking$ranking.random else ranking$ranking.common
        p_scores <- if(!is.null(ranking$Pscore.random)) ranking$Pscore.random else ranking$Pscore.common
        validate(need(!is.null(rank_scores) && !is.null(p_scores) && length(rank_scores) > 0, "Ranking scores are unavailable."))
        ranking_df <- data.frame(
          Treatment = names(rank_scores),
          PScore = as.numeric(p_scores)
        ) %>% arrange(desc(PScore))
        write.csv(ranking_df, file, row.names = FALSE)
      }
    )
    # --- Add other results outputs as needed ---
    
  }) # End moduleServer
}



# ---------------------------
# 7. Inconsistency Module (No changes needed from previous version)
# ---------------------------
inconsistencyUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Design-based Inconsistency (Decomposition of Q)"),
    verbatimTextOutput(ns("decompText")),
    br(),
    downloadButton(ns("dlDecompText"), "Download Inconsistency Summary as TXT")
  )
}

inconsistencyServer <- function(id, analysisResult) {
  moduleServer(id, function(input, output, session) {
    decomp_res <- reactive({
      res <- analysisResult()
      req(res, res$model) # Ensure result and model are available
      if (res$method == "freq") {
        tryCatch(decomp.design(res$model), error=function(e) e)
      } else {
        NULL
      }
    })
    
    output$decompText <- renderPrint({
      res_decomp <- decomp_res()
      if(is.null(res_decomp)) {
        cat("Inconsistency analysis not applicable or model not ready.")
      } else if(inherits(res_decomp, "error")) {
        cat("Error running inconsistency analysis:", res_decomp$message)
      } else {
        print(res_decomp)
      }
    })
    
    output$dlDecompText <- downloadHandler(
      filename = function() { "InconsistencySummary.txt" },
      content = function(file) {
        res_decomp <- decomp_res()
        if(is.null(res_decomp)) {
          writeLines("Inconsistency analysis not applicable or model not ready.", con = file)
        } else if(inherits(res_decomp, "error")) {
          writeLines(paste("Error running inconsistency analysis:", res_decomp$message), con = file)
        } else {
          txt <- capture.output(print(res_decomp))
          writeLines(txt, con = file)
        }
      }
    )
  })
}


# ---------------------------
# 8. Example CSV Files Module (Revised for Binary/Ratio)
# ---------------------------
exampleCSVUI <- function(id) {
  ns <- NS(id)
  tagList(
    h3("Example CSV Files"),
    h4("Arm-level Binary Data Example (Events/N)"),
    verbatimTextOutput(ns("armCSVBin")),
    downloadButton(ns("dlArmCSVBin"), "Download Binary Arm-level CSV"),
    br(), br(),
    h4("Contrast-level Log-Ratio Data Example (Log-RR/Log-OR/Log-HR)"),
    verbatimTextOutput(ns("contrastCSV")),
    downloadButton(ns("dlContrastCSV"), "Download Contrast-level CSV")
  )
}

exampleCSVServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # Sample for Arm-level Binary (event/n)
    armCSVBin <- "study,treatment,event,n
Study1,DrugA,10,50
Study1,Placebo,15,50
Study2,DrugB,8,60
Study2,Placebo,12,60
Study3,DrugC,5,55
Study3,Placebo,10,55
Study4,DrugA,7,40
Study4,DrugB,10,40
Study5,DrugC,6,45
Study5,DrugD,9,45"
    
    # Sample for Contrast-level (Log-Ratio TE/seTE)
    contrastCSV <- "studlab,treat1,treat2,TE,seTE
Study1,DrugA,Placebo,-0.51,0.37
Study2,DrugB,Placebo,-0.41,0.41
Study3,DrugC,Placebo,-0.74,0.51
Study4,DrugA,DrugB,-0.36,0.52
Study5,DrugC,DrugD,-0.45,0.60" # Example log-ratios
    
    output$armCSVBin <- renderText({ armCSVBin })
    output$contrastCSV <- renderText({ contrastCSV })
    
    output$dlArmCSVBin <- downloadHandler(
      filename = function() { "example_binary_arm_data.csv" },
      content = function(file) {
        writeLines(armCSVBin, con = file)
      }
    )
    
    output$dlContrastCSV <- downloadHandler(
      filename = function() { "example_contrast_data.csv" },
      content = function(file) {
        writeLines(contrastCSV, con = file)
      }
    )
  })
}


# ---------------------------
# 9. Main App UI & Server (bs4Dash)
# ---------------------------
ui <- bs4DashPage(
  title = "786-MIII NMA Shiny App (Frequentist - Ratio Outcomes)", # Updated title
  header = bs4DashNavbar(title = "786-MIII NMA Shiny App"),
  sidebar = bs4DashSidebar(
    bs4SidebarMenu(
      id = "sidebarMenu", # Added id for observing tab changes
      bs4SidebarMenuItem("Data Upload", tabName = "dataUpload", icon = icon("upload")),
      bs4SidebarMenuItem("Network Plot", tabName = "network", icon = icon("project-diagram")),
      bs4SidebarMenuItem("Analysis Setup", tabName = "analysis", icon = icon("sliders-h")), # Renamed
      bs4SidebarMenuItem("Diagnostics", tabName = "diagnostics", icon = icon("stethoscope")),
      bs4SidebarMenuItem("Results", tabName = "results", icon = icon("table")),
      bs4SidebarMenuItem("Inconsistency", tabName = "inconsistency", icon = icon("exclamation-triangle")),
      bs4SidebarMenuItem("Example CSV Files", tabName = "exampleCSV", icon = icon("file-alt"))
    )
  ),
  body = bs4DashBody(
    bs4TabItems(
      bs4TabItem(tabName = "dataUpload", dataUploadUI("dataUpload")),
      bs4TabItem(tabName = "network", networkUI("network")),
      bs4TabItem(tabName = "analysis", analysisUI("analysis")),
      bs4TabItem(tabName = "diagnostics", diagnosticsUI("diagnostics")),
      bs4TabItem(tabName = "results", resultsUI("results")),
      bs4TabItem(tabName = "inconsistency", inconsistencyUI("inconsistency")),
      bs4TabItem(tabName = "exampleCSV", exampleCSVUI("exampleCSV"))
    )
  )
)

server <- function(input, output, session) {
  # Instantiate modules and pass reactive summary measure from data upload to analysis
  data_output <- dataUploadServer("dataUpload")
  
  # Create a reactive expression specifically for the selected summary measure
  selected_sm_reactive <- reactive({
    # Need to access the specific input from the module UI
    # Use %||% to provide a default if the input doesn't exist yet (e.g., before UI renders)
    input$`dataUpload-summaryMeasure` %||% "HR" # Default to HR if binary UI not shown or not selected
  })
  
  analysisRes <- analysisServer("analysis", data = data_output, selected_sm = selected_sm_reactive) # Pass the reactive expression
  
  networkServer("network", data = data_output)
  resultsServer("results", analysisResult = analysisRes)
  diagnosticsServer("diagnostics", analysisResult = analysisRes)
  inconsistencyServer("inconsistency", analysisResult = analysisRes)
  exampleCSVServer("exampleCSV")
}


shinyApp(ui = ui, server = server)
