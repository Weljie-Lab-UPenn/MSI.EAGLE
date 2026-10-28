MassCalibrationUI <- function(id) {
  ns <- NS(id)
  tabPanel("Mass Calibration", value=ns("tab"), sidebarLayout(sidebarPanel(
    selectInput(ns("dataset_file"),"Saved dataset in working directory",choices=character(0)),
    actionButton(ns("load"),"Load dataset"),
    selectInput(ns("polarity"),"Acquisition polarity",c("Select acquisition polarity"="unknown","Positive"="positive","Negative"="negative")),
    uiOutput(ns("adduct_controls")),
    fileInput(ns("references"),"Custom reference TSV (optional)",accept=".tsv"),
    helpText("Custom columns: name, adduct, formula, theo. Competing ion compositions are checked automatically."),
    checkboxInput(ns("extended"),"Extended low-mass search and models",FALSE),
    actionButton(ns("find"),"Find anchors"),
    checkboxInput(ns("show_all_matches"),"Show weaker and ambiguous matches",FALSE),
    helpText("Click candidate rows in the anchor table to select or deselect them, then click Compare fits. Only selected candidates that pass the quality checks are used."),
    selectInput(ns("model"),"Model selection",c("Automatic"="auto","Unchanged"="unchanged","Constant"="constant","Linear"="linear","Quadratic"="quadratic","Cubic"="cubic","Inverse + linear"="inverse_linear","Inverse + quadratic"="inverse_quadratic","Log quadratic"="log_quadratic","Log cubic"="log_cubic")),
    actionButton(ns("fit"),"Compare fits"),
    actionButton(ns("apply"),"Apply to new dataset"),
    actionButton(ns("undo"),"Restore original masses"),
    actionButton(ns("open"),"Open result in Peak Picking"),
    shinyFiles::shinySaveButton(ns("save_imzml"),"Save calibrated imzML bundle","Save",filetype=list("")), width=3),
    mainPanel(verbatimTextOutput(ns("status")),DT::dataTableOutput(ns("anchors")),
      DT::dataTableOutput(ns("models")),plotOutput(ns("errors")),
      helpText("LOO is an internal model-selection diagnostic after anchor filtering. Missing isotope evidence is shown explicitly. Correction is held constant outside the anchor range. This workflow applies one shared feature axis to all runs."),
      plot_card_UI(ns("preview")),width=9)))
}
