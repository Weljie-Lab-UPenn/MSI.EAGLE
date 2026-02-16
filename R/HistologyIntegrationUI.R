HistologyIntegrationUI <- function(id) {
  ns <- NS(id)
  polygon_color_input <- if (requireNamespace("colourpicker", quietly = TRUE)) {
    colourpicker::colourInput(
      ns("polygon_outline_color"),
      "Polygon outline color",
      value = "#73FFFF",
      allowTransparent = FALSE
    )
  } else {
    textInput(ns("polygon_outline_color"), "Polygon outline color (hex)", value = "#73FFFF")
  }

  tabPanel(
    "Histology Overlay",
    value = ns("tab9"),
    tagList(
      tags$style(HTML("
        .compact-controls .form-group { margin-bottom: 6px; }
        .compact-controls .control-label { margin-bottom: 2px; font-size: 12px; }
        .compact-controls .btn { margin-right: 4px; margin-bottom: 4px; }
        .nudge-pad .btn { min-width: 42px; padding: 4px 8px; }
      ")),
      sidebarLayout(
        sidebarPanel(
          width = 3,
          div(
            class = "compact-controls",
            tags$h5("Inputs"),
            radioButtons(
              ns("msi_source"),
              "MSI Source",
              choices = c(
                "Use current MSI data" = "current",
                "Upload MSI .rds file" = "upload"
              ),
              selected = "current"
            ),
            uiOutput(ns("msi_upload_ui")),
            fileInput(
              ns("histology_upload"),
              "Histology image",
              accept = c(".png", ".jpg", ".jpeg", ".tif", ".tiff")
            ),
            fileInput(
              ns("cluster_overlay_upload"),
              "Cluster-color image (optional)",
              accept = c(".png", ".jpg", ".jpeg", ".tif", ".tiff")
            ),
            fileInput(
              ns("polygon_file"),
              "Polygon file (.geojson/.json)",
              accept = c(".geojson", ".json")
            ),
            selectizeInput(
              ns("mz_select"),
              "Select m/z",
              choices = character(0),
              selected = NULL,
              options = list(placeholder = "Load MSI data to populate m/z list")
            ),
            selectInput(
              ns("intensity_transform"),
              "Intensity transform",
              choices = c("none", "sqrt", "log1p", "asinh"),
              selected = "asinh"
            ),
            selectInput(
              ns("msi_palette"),
              "MSI palette",
              choices = c("Inferno", "Magma", "Plasma", "Viridis", "Cividis", "YlOrRd", "BluYl", "Temps", "Terrain 2"),
              selected = "Inferno"
            ),
            selectInput(
              ns("cluster_palette"),
              "Discrete palette",
              choices = c("Set 2", "Dark 3", "Dynamic", "Warm", "Cold", "Harmonic"),
              selected = "Set 2"
            ),
            checkboxInput(ns("enhance_contrast"), "Enhance contrast", value = TRUE),
            checkboxInput(ns("gaussian_smooth"), "Gaussian smoothing", value = TRUE),
            numericInput(ns("gaussian_sigma"), "Gaussian sigma", value = 1, min = 0.1, step = 0.1),
            checkboxInput(ns("flip_histology_y"), "Flip histology Y", value = FALSE),
            tags$hr(),
            tags$h5("Mapping & Save"),
            radioButtons(
              ns("mapping_source"),
              "Mapping source",
              choices = c("Cluster image" = "cluster", "Polygon GeoJSON" = "polygon"),
              selected = "polygon"
            ),
            uiOutput(ns("mapping_source_ui")),
            actionButton(ns("map_to_pdata"), "Map selected source to pData"),
            shinyFiles::shinySaveButton(ns("save_mapped_imzml"), "Save mapped imzML", "Save", filetype = list("imzML" = "imzML")),
            tags$hr(),
            uiOutput(ns("pdata_field_ui"))
          )
        ),
        mainPanel(
          width = 9,
          checkboxInput(ns("show_overlay_info"), "Show overlay details", value = FALSE),
          uiOutput(ns("overlay_info_ui")),
          radioButtons(
            ns("overlay_layer"),
            "Overlay shown",
            choices = c("Polygon outlines" = "polygon", "Histology image" = "histology", "Cluster-color image" = "cluster"),
            selected = "polygon",
            inline = TRUE
          ),
          plotOutput(ns("overlay_plot"), height = "680px"),
          tags$hr(),
          div(
            class = "compact-controls",
            tags$h5("Registration"),
            fluidRow(
              column(
                3,
                selectInput(
                  ns("overlay_scale_mode"),
                  "Overlay scale mode",
                  choices = c(
                    "Absolute ratio (microscopy px to MSI px)" = "absolute",
                    "Relative to fit-to-canvas" = "fit"
                  ),
                  selected = "absolute"
                ),
                numericInput(ns("histology_um_per_px"), "Histology um/pixel", value = 0.1138, min = 1e-06, step = 1e-04),
                numericInput(ns("msi_um_per_px"), "MSI um/pixel", value = 2, min = 1e-06, step = 0.1),
                numericInput(ns("histology_resample_factor"), "Histology export factor", value = 1, min = 0.01, step = 0.01),
                actionButton(ns("set_scale_from_resolution"), "Set scale = histology/MSI")
              ),
              column(
                3,
                sliderInput(ns("scale_x"), "Scale X", min = 0.001, max = 50, value = 1, step = 0.0005),
                sliderInput(ns("scale_y"), "Scale Y", min = 0.001, max = 50, value = 1, step = 0.0005),
                sliderInput(ns("rotate_deg"), "Rotation (degrees)", min = -180, max = 180, value = 0, step = 0.5),
                sliderInput(ns("translate_x"), "Translate X", min = -1000, max = 1000, value = 0, step = 1),
                sliderInput(ns("translate_y"), "Translate Y", min = -1000, max = 1000, value = 0, step = 1)
              ),
              column(
                3,
                checkboxInput(ns("fine_move"), "Fine move (1; otherwise 5)", value = FALSE),
                div(
                  class = "nudge-pad",
                  tags$table(
                    tags$tr(tags$td(""), tags$td(actionButton(ns("move_up"), "\u2191")), tags$td("")),
                    tags$tr(
                      tags$td(actionButton(ns("move_left"), "\u2190")),
                      tags$td(tags$span("Nudge")),
                      tags$td(actionButton(ns("move_right"), "\u2192"))
                    ),
                    tags$tr(tags$td(""), tags$td(actionButton(ns("move_down"), "\u2193")), tags$td(""))
                  )
                ),
                numericInput(ns("translate_x_num"), "Translate X (num)", value = 0, min = -1000, max = 1000, step = 0.1),
                numericInput(ns("translate_y_num"), "Translate Y (num)", value = 0, min = -1000, max = 1000, step = 0.1)
              ),
              column(
                3,
                selectInput(
                  ns("polygon_axis_mode"),
                  "Polygon axis mode",
                  choices = c("Auto (infer 90-degree swap)" = "auto", "Standard (X,Y)" = "xy", "Swap axes (Y,X)" = "yx"),
                  selected = "yx"
                ),
                checkboxInput(ns("polygon_color_by_label"), "Polygon color by label", value = FALSE),
                polygon_color_input,
                sliderInput(ns("histology_alpha"), "Histology alpha", min = 0, max = 1, value = 0.5, step = 0.01),
                sliderInput(ns("cluster_alpha"), "Cluster alpha", min = 0, max = 1, value = 0.7, step = 0.01)
              )
            ),
            fluidRow(
              column(3, numericInput(ns("scale_x_num"), "Scale X (num)", value = 1, min = 0.001, max = 50, step = 0.0005)),
              column(3, numericInput(ns("scale_y_num"), "Scale Y (num)", value = 1, min = 0.001, max = 50, step = 0.0005)),
              column(3, numericInput(ns("rotate_deg_num"), "Rotation (num)", value = 0, min = -180, max = 180, step = 0.1)),
              column(3, actionButton(ns("reset_transform"), "Reset transform"))
            ),
            fluidRow(
              column(3, numericInput(ns("histology_alpha_num"), "Histology alpha (num)", value = 0.5, min = 0, max = 1, step = 0.01)),
              column(3, numericInput(ns("cluster_alpha_num"), "Cluster alpha (num)", value = 0.7, min = 0, max = 1, step = 0.01)),
              column(3, downloadButton(ns("download_registration_params"), "Save params (.txt)")),
              column(3, fileInput(ns("registration_params_upload"), "Load params (.txt/.csv)", accept = c(".txt", ".csv")))
            )
          ),
          tags$hr(),
          plotOutput(ns("pdata_plot"), height = "460px"),
          tags$h5("Mapping Lookup"),
          DT::dataTableOutput(ns("cluster_lookup_table")),
          tags$h5("pData Preview (first 200 rows)"),
          DT::dataTableOutput(ns("pdata_preview"))
        )
      )
    )
  )
}
