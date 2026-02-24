### R/CorrelationUI.R
CorrelationUI <- function(id) {
  ns <- NS(id)
  
  tabPanel("Colocalization analysis",
           value="tab8",
           sidebarLayout(
             sidebarPanel(
               radioButtons(ns("corr_source"), "Type of data", choices = list(
                 "Read imzML / .rds"="from_file",
                 "From Stats"="from_stats"
               ), selected = "from_stats"),
               uiOutput(ns("source_ui")),
               
               uiOutput(ns("corr_vars"), label="select_variable_for_correlation"),
               uiOutput(ns("corr_mz"), label="m/z selection for correlation"),
               radioButtons(ns("corr_sort"), "Sorted by", choices = list(
                 "Pearson's correlation"="cor",
                 "Manders overlap coefficient"="MOC",
                 "Manders’ colocalization coefficients (M1)"="M1",
                 "Manders’ colocalization coefficients (M2)"="M2",
                 "Dice similarity coefficient"="Dice"
               ), selected = "cor"),
               radioButtons(ns("corr_table_order"), "Result table order", choices = list(
                 "Descending" = "desc",
                 "Ascending" = "asc"
               ), selected = "desc", inline = TRUE),
               checkboxInput(ns("corr_full_table"), "Compute full correlation table (slower, needed to inspect negative Pearson values)", value = FALSE),
               numericInput(ns("corr_plot_n"), "Number of top correlations to plot", min=1, value=3),
               checkboxInput(ns("corr_plot_include_source"), "Include source ion in plotted panels", value = FALSE),
               checkboxInput(ns("corr_enhance_hist"), "Histogram contrast enhancement", value = TRUE),
               checkboxInput(ns("corr_smooth_gaussian"), "Gaussian smoothing", value = FALSE),
               tags$div(style = "font-size: 0.9em; color: #666; margin-bottom: 8px;",
                        "Correlation results only update when 'Generate colocalization' is clicked."),
               actionButton(ns("action_corr"), "Generate colocalization")
               
               
             ),
             mainPanel(
               tags$div(
                 style = "overflow-x:auto; max-width:100%; border:1px solid #eee; padding:6px; background:#fafafa;",
                 imageOutput(ns("plot15"), width="100%", height="600px")
               ),
               tags$hr(),
               fluidRow(
                 column(4, actionButton(ns("corr_use_selected"), "Use selected row m/z")),
                 column(4, actionButton(ns("corr_rerun_selected"), "Re-run from selected row"))
               ),
               tags$br(),
               DT::DTOutput(ns("corr_table"))
               
             )
             
           )
           
  )
  
  
}
