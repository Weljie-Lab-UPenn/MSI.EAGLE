MzIdentificationUI <- function(id) {
  ns <- NS(id)
  tabPanel("m/z Identification", value = ns("tab"),
    sidebarLayout(
      sidebarPanel(
        h4("Dataset"),
        selectInput(ns("dataset_file"), "Saved dataset in working directory", choices = character(0)),
        actionButton(ns("load_dataset"), "Load saved dataset"),
        actionButton(ns("use_calibration_result"), "Use current Mass Calibration result"),
        helpText("Load a calibrated .imzML bundle or .rds dataset. The saved dataset becomes the active source for identification and Statistics."),
        h4("Local libraries"),
        uiOutput(ns("library_status")),
        actionButton(ns("install_hmdb"), "Install / update HMDB"),
        actionButton(ns("open_hmdb"), "Open HMDB download page"),
        actionButton(ns("install_hmdb_downloaded"), "Install HMDB from Downloads"),
        actionButton(ns("open_lmsd"), "Open LIPID MAPS download page"),
        actionButton(ns("install_lmsd_downloaded"), "Install LIPID MAPS from Downloads"),
        helpText("Installers build compact local indexes from official releases. If a provider requires browser access, complete its official download, then install the downloaded archive; no upload is required."),
        uiOutput(ns("library_select")),
        selectInput(ns("polarity"), "Polarity", choices = c("Positive" = "positive", "Negative" = "negative")),
        uiOutput(ns("adducts")),
        numericInput(ns("ppm"), "Mass tolerance (ppm)", value = 5, min = 0.1, max = 100, step = 0.1),
        actionButton(ns("match"), "Match corrected m/z values"),
        actionButton(ns("accept"), "Accept selected candidate"),
        shinyFiles::shinySaveButton(ns("save_annotated"), "Save annotated imzML bundle", "Save", filetype = list("")),
        helpText("All candidates in every selected database are retained. LIPID MAPS candidates are shown at sum-composition level because MS1 MSI cannot resolve individual acyl chains. An MS1 match is a tentative annotation, not a confirmed identification."),
        width = 3),
      mainPanel(
        verbatimTextOutput(ns("status")),
        DT::dataTableOutput(ns("features")),
        h4("Candidates for selected feature"),
        DT::dataTableOutput(ns("candidates")),
        width = 9)))
}
