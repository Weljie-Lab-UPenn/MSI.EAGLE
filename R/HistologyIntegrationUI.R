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
            fileInput(
              ns("roi_polygon_file"),
              "ROI rectangle file (optional)",
              accept = c(".geojson", ".json")
            ),
            fileInput(
              ns("nucleus_polygon_file"),
              "Nucleus polygon file (optional)",
              accept = c(".geojson", ".json")
            ),
            fileInput(
              ns("histology_metadata_file"),
              "Histology metadata sidecar (optional .json)",
              accept = c(".json")
            ),
            tags$small("If provided, polygon mapping will also annotate nucleus-versus-cytoplasm pixels in pData."),
            selectInput(
              ns("msi_plot_mode"),
              "MSI display mode",
              choices = c(
                "Single m/z ion" = "mz",
                "RGB ions (2-3)" = "rgb",
                "pData field" = "pdata"
              ),
              selected = "mz"
            ),
            selectizeInput(
              ns("mz_select"),
              "Select m/z",
              choices = character(0),
              selected = NULL,
              options = list(placeholder = "Load MSI data to populate m/z list")
            ),
            uiOutput(ns("msi_mode_controls")),
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
              choices = c("Alphabet", "Set 2", "Dark 3", "Dynamic", "Warm", "Cold", "Harmonic"),
              selected = "Alphabet"
            ),
            checkboxInput(ns("enhance_contrast"), "Enhance contrast", value = TRUE),
            checkboxInput(ns("gaussian_smooth"), "Gaussian smoothing", value = FALSE),
            numericInput(ns("gaussian_sigma"), "Gaussian sigma", value = 1, min = 0.1, step = 0.1),
            uiOutput(ns("registration_frame_status_ui")),
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
            tags$div(style = "margin-top: 8px;", actionButton(ns("reset_transform"), "Reset transform")),
            tags$div(style = "margin-top: 8px;", downloadButton(ns("download_registration_params"), "Save params (.txt)")),
            tags$div(
              style = "margin-top: 8px;",
              fileInput(ns("registration_params_upload"), "Load params (.txt/.csv)", accept = c(".txt", ".csv"))
            ),
            shinyFiles::shinySaveButton(ns("save_mapped_imzml"), "Save mapped imzML", "Save", filetype = list("")),
            tags$hr(),
            uiOutput(ns("pdata_field_ui")),
            checkboxInput(ns("hide_pdata_legend"), "Hide pData legend", value = TRUE)
          )
        ),
        mainPanel(
          width = 9,
          checkboxInput(ns("show_overlay_info"), "Show overlay details", value = FALSE),
          uiOutput(ns("overlay_info_ui")),
          radioButtons(
            ns("overlay_layer"),
            "Overlay shown",
            choices = c("Combined: histology + polygons" = "combined", "Polygon outlines" = "polygon", "Histology image" = "histology", "Cluster-color image" = "cluster"),
            selected = "combined",
            inline = TRUE
          ),
          plotOutput(ns("overlay_plot"), height = "680px"),
          downloadButton(ns("download_overlay_pdf"), "Download overlay (PDF)"),
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
                numericInput(ns("scale_x_num"), "Scale X (num)", value = 1, min = 0.001, max = 50, step = 0.0005),
                sliderInput(ns("scale_y"), "Scale Y", min = 0.001, max = 50, value = 1, step = 0.0005),
                numericInput(ns("scale_y_num"), "Scale Y (num)", value = 1, min = 0.001, max = 50, step = 0.0005),
                sliderInput(ns("rotate_deg"), "Rotation (degrees)", min = -180, max = 180, value = 0, step = 0.5),
                numericInput(ns("rotate_deg_num"), "Rotation (num)", value = 0, min = -180, max = 180, step = 0.1),
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
                checkboxInput(ns("show_advanced_registration"), "Show advanced registration controls", value = FALSE),
                conditionalPanel(
                  condition = sprintf("input['%s']", ns("show_advanced_registration")),
                  selectInput(
                    ns("polygon_axis_mode"),
                    "Polygon axis mode",
                    choices = c("Auto (infer 90-degree swap)" = "auto", "Standard (X,Y)" = "xy", "Swap axes (Y,X)" = "yx"),
                    selected = "auto"
                  ),
                  checkboxInput(ns("flip_histology_y"), "Flip histology Y", value = FALSE)
                ),
                checkboxInput(ns("polygon_color_by_label"), "Polygon color by label", value = FALSE),
                polygon_color_input,
                sliderInput(ns("polygon_linewidth"), "Polygon line width", min = 0.2, max = 6, value = 1, step = 0.1),
                sliderInput(ns("histology_alpha"), "Histology alpha", min = 0, max = 1, value = 0.5, step = 0.01),
                numericInput(ns("histology_alpha_num"), "Histology alpha (num)", value = 0.5, min = 0, max = 1, step = 0.01),
                sliderInput(ns("cluster_alpha"), "Cluster alpha", min = 0, max = 1, value = 0.7, step = 0.01),
                numericInput(ns("cluster_alpha_num"), "Cluster alpha (num)", value = 0.7, min = 0, max = 1, step = 0.01)
              )
            ),
            checkboxInput(ns("show_fit_info"), "Show fit diagnostics/info", value = TRUE),
            tabsetPanel(
              id = ns("registration_fit_tabs"),
              type = "tabs",
              selected = "histology_fit",
              tabPanel(
                "Histology Fit",
                value = "histology_fit",
                tags$small("Local XY refinement using the transformed histology image against the current MSI spatial signal. Current scale and rotation are held fixed."),
                fluidRow(
                  column(2, numericInput(ns("histology_fit_range"), "Range", value = 30, min = 1, max = 500, step = 1)),
                  column(2, numericInput(ns("histology_fit_step"), "Step", value = 2, min = 1, max = 50, step = 1)),
                  column(
                    3,
                    selectInput(
                      ns("histology_feature_mode"),
                      "Histology feature",
                      choices = c(
                        "Hematoxylin-like + darkness" = "hematoxylin",
                        "Darkness / tissue density" = "darkness",
                        "Purple chromatin emphasis" = "purple"
                      ),
                      selected = "purple"
                    )
                  ),
                  column(
                    3,
                    selectInput(
                      ns("histology_fit_signal_source"),
                      "MSI target signal",
                      choices = c(
                        "Current MSI display" = "current",
                        "Fused MSI (top spatial ions, default)" = "multi",
                        "PCA component (self-contained)" = "pca",
                        "pData field (numeric or categorical)" = "pdata"
                      ),
                      selected = "multi"
                    )
                  ),
                  column(
                    2,
                    tags$div(
                      style = "margin-top: 24px;",
                      actionButton(ns("run_histology_fit"), "Run Histology Fit")
                    )
                  )
                ),
                fluidRow(
                  column(4, uiOutput(ns("histology_fit_pdata_field_ui"))),
                  column(
                    2,
                    selectInput(
                      ns("histology_fit_intensity_relation"),
                      "Intensity relation",
                      choices = c(
                        "Either / unsigned (recommended)" = "either",
                        "Inverse" = "inverse",
                        "Direct" = "direct"
                      ),
                      selected = "inverse"
                    )
                  ),
                  column(4, uiOutput(ns("histology_fit_preview_ui"))),
                  column(
                    2,
                    tags$div(
                      style = "margin-top: 24px;",
                      actionButton(ns("apply_histology_fit_choice"), "Apply selected")
                    )
                  )
                ),
                conditionalPanel(
                  condition = sprintf("input['%s']", ns("show_fit_info")),
                  fluidRow(
                    column(6, uiOutput(ns("histology_fit_heatmap_ui"))),
                    column(6, plotOutput(ns("histology_fit_contour"), height = "240px"))
                  ),
                  fluidRow(
                    column(8, plotOutput(ns("histology_fit_target_plot"), height = "240px")),
                    column(
                      4,
                      uiOutput(ns("histology_fit_target_info_ui")),
                      tags$div(style = "margin-top: 8px;", downloadButton(ns("download_histology_fit_target_png"), "Download target (.png)"))
                    )
                  ),
                  verbatimTextOutput(ns("histology_fit_summary"))
                )
              ),
              tabPanel(
                "Stat Fit",
                value = "stat_fit",
                tags$small("Color-map independent fit using inside-vs-outside or polygon-cluster-group ion statistics."),
                fluidRow(
                  column(2, numericInput(ns("stat_fit_range"), "Range", value = 20, min = 1, max = 300, step = 1)),
                  column(2, numericInput(ns("stat_fit_step"), "Step", value = 2, min = 1, max = 50, step = 1)),
                  column(2, numericInput(ns("stat_fit_buffer_px"), "Outside buffer (px)", value = 3, min = 1, max = 25, step = 1)),
                  column(2, numericInput(ns("stat_fit_min_pixels"), "Min pixels/group", value = 10, min = 2, max = 1000, step = 1)),
                  column(2, numericInput(ns("stat_fit_max_ions"), "Max ions", value = 200, min = 20, max = 2000, step = 10)),
                  column(2, numericInput(ns("stat_fit_alpha"), "Alpha", value = 0.1, min = 0.001, max = 0.5, step = 0.01))
                ),
                fluidRow(
                  column(
                    3,
                    selectInput(
                      ns("stat_fit_metric_mode"),
                      "Stat-fit metric",
                      choices = c(
                        "Inside vs outside polygon (t-test)" = "inside_outside",
                        "Polygon cluster groups (ANOVA; exclude outside/unassigned)" = "polygon_cluster_groups"
                      ),
                      selected = "inside_outside"
                    )
                  ),
                  column(
                    4,
                    selectInput(
                      ns("stat_fit_outside_mode"),
                      "Outside group",
                      choices = c(
                        "BBox complement (recommended)" = "bbox",
                        "Local ring around polygon" = "local",
                        "Global complement (all non-polygon pixels)" = "global"
                      ),
                      selected = "bbox"
                    )
                  ),
                  column(3, selectInput(
                    ns("stat_fit_objective"),
                    "Stat objective",
                    choices = c("Minimize score" = "min", "Maximize score" = "max"),
                      selected = "min"
                  )),
                  column(2, numericInput(ns("stat_fit_bbox_pad"), "BBox pad (px)", value = 25, min = 0, max = 500, step = 1))
                ),
                fluidRow(
                  column(6, uiOutput(ns("stat_fit_group_field_ui")))
                ),
                fluidRow(
                  column(
                    2,
                    tags$div(
                      style = "margin-top: 24px;",
                      actionButton(ns("run_stat_fit"), "Run Stat Fit")
                    )
                  ),
                  column(4, uiOutput(ns("stat_fit_preview_ui"))),
                  column(
                    2,
                    tags$div(
                      style = "margin-top: 24px;",
                      actionButton(ns("apply_stat_fit_choice"), "Apply selected")
                    )
                  ),
                  column(
                    2,
                    tags$div(
                      style = "margin-top: 24px;",
                      checkboxInput(ns("stat_fit_use_adjusted"), "Use BH-adjusted p", value = TRUE)
                    )
                  ),
                  column(
                    2,
                    tags$div(
                      style = "margin-top: 24px;",
                      checkboxInput(ns("stat_fit_use_abs_lfc"), "Use |log2FC|", value = TRUE)
                    )
                  )
                ),
                fluidRow(
                  column(2, numericInput(ns("stat_fit_top_n"), "Top features", value = 10, min = 1, max = 100, step = 1)),
                  column(
                    3,
                    tags$div(
                      style = "margin-top: 24px;",
                      actionButton(ns("stat_fit_max_info"), "Max info to console")
                    )
                  )
                ),
                conditionalPanel(
                  condition = sprintf("input['%s']", ns("show_fit_info")),
                  fluidRow(
                    column(6, uiOutput(ns("stat_fit_heatmap_ui"))),
                    column(6, plotOutput(ns("stat_fit_contour"), height = "240px"))
                  ),
                  verbatimTextOutput(ns("stat_fit_summary"))
                )
              ),
              tabPanel(
                "Edge Fit",
                value = "edge_fit",
                fluidRow(
                  column(2, numericInput(ns("optimize_xy_range"), "Range", value = 40, min = 1, max = 500, step = 1)),
                  column(2, numericInput(ns("optimize_xy_step"), "Step", value = 2, min = 1, max = 50, step = 1)),
                  column(2, numericInput(ns("optimize_edge_band"), "Boundary band (px)", value = 4, min = 1, max = 20, step = 1)),
                  column(
                    2,
                    tags$div(
                      style = "margin-top: 24px;",
                      actionButton(ns("optimize_xy"), "Auto-fit XY")
                    )
                  ),
                  column(2, uiOutput(ns("optimize_xy_preview_ui"))),
                  column(
                    2,
                    tags$div(
                      style = "margin-top: 24px;",
                      actionButton(ns("apply_optimize_xy_choice"), "Apply selected")
                    )
                  )
                ),
                fluidRow(
                  column(
                    4,
                    selectInput(
                      ns("edge_fit_signal_source"),
                      "Edge-fit optimization signal",
                      choices = c(
                        "Current MSI display (default)" = "current",
                        "pData field (numeric or categorical)" = "pdata"
                      ),
                      selected = "current"
                    )
                  ),
                  column(8, uiOutput(ns("edge_fit_pdata_field_ui")))
                )
              ),
              tabPanel(
                "Polygon clustering",
                value = "polygon_clustering",
                tags$small("Cluster polygons using QuPath measurements (cells/nuclei) or geometry-only features for manual polygons."),
                fluidRow(
                  column(
                    4,
                    fileInput(
                      ns("polygon_cluster_table"),
                      "External measurements table (optional override)",
                      accept = c(".csv", ".tsv", ".txt")
                    ),
                    tags$small("By default, clustering uses attributes already stored in the loaded polygon GeoJSON."),
                    selectInput(
                      ns("polygon_cluster_feature_mode"),
                      "Feature source",
                      choices = c(
                        "Auto (GeoJSON attributes if available, else geometry)" = "auto",
                        "Attributes only" = "measurements",
                        "Attributes + geometry" = "measurements_plus_geometry",
                        "Geometry only (manual polygons)" = "geometry"
                      ),
                      selected = "auto"
                    ),
                    uiOutput(ns("polygon_cluster_id_field_ui")),
                    uiOutput(ns("polygon_cluster_features_ui"))
                  ),
                  column(
                    4,
                    uiOutput(ns("polygon_cluster_join_field_ui")),
                    selectInput(
                      ns("polygon_cluster_join_mode"),
                      "Attribute-to-polygon matching",
                      choices = c(
                        "Exact match on selected IDs" = "exact",
                        "Row order (same polygon order in file/table)" = "row_order"
                      ),
                      selected = "exact"
                    ),
                    numericInput(ns("polygon_cluster_k"), "Number of clusters (k)", value = 4, min = 2, max = 50, step = 1),
                    numericInput(ns("polygon_cluster_keep_prop"), "Keep proportion per cluster (distinctive polygons)", value = 1, min = 0.01, max = 1, step = 0.05),
                    numericInput(ns("polygon_cluster_pca_dims"), "PCA dimensions", value = 4, min = 1, max = 20, step = 1),
                    numericInput(ns("polygon_cluster_seed"), "Random seed", value = 123, min = 1, step = 1),
                    numericInput(ns("polygon_cluster_nstart"), "k-means nstart", value = 25, min = 1, max = 200, step = 1),
                    checkboxInput(ns("polygon_cluster_scale"), "Scale features", value = TRUE),
                    checkboxInput(ns("polygon_cluster_drop_unmatched"), "Exclude unmatched polygons from clustering", value = FALSE)
                  ),
                  column(
                    4,
                    textInput(ns("polygon_cluster_pdata_col"), "Clustered pData column name", value = "polygon_cluster_class"),
                    tags$div(style = "margin-top: 24px;", actionButton(ns("run_polygon_clustering"), "Run polygon clustering")),
                    fluidRow(
                      column(6, numericInput(ns("polygon_cluster_outlier_sd"), "Outlier SD threshold", value = 2, min = 0.5, max = 10, step = 0.25)),
                      column(6, tags$div(style = "margin-top: 25px;", actionButton(ns("remove_polygon_cluster_outliers"), "Remove PCA outliers")))
                    ),
                    fluidRow(
                      column(6, numericInput(ns("polygon_cluster_profile_top_n"), "Top defining features/cluster", value = 10, min = 1, max = 100, step = 1)),
                      column(6, tags$div(style = "margin-top: 25px;", actionButton(ns("summarize_polygon_clusters"), "Summarize cluster features")))
                    ),
                    tags$div(style = "margin-top: 6px;", downloadButton(ns("download_polygon_cluster_profile"), "Download summary (.csv)")),
                    tags$div(style = "margin-top: 8px;", actionButton(ns("use_polygon_clusters_for_labels"), "Use clusters for overlay/mapping")),
                    tags$small("Then set overlay to polygon + color-by-label, or use the existing polygon mapping button to write clustered labels into pData."),
                    tags$br(), tags$br(),
                    verbatimTextOutput(ns("polygon_cluster_summary"))
                  )
                ),
                fluidRow(
                  column(6, plotOutput(ns("polygon_cluster_plot"), height = "260px")),
                  column(6, DT::dataTableOutput(ns("polygon_cluster_counts_table")))
                ),
                DT::dataTableOutput(ns("polygon_cluster_preview_table")),
                tags$h5("Cluster-defining feature summary"),
                DT::dataTableOutput(ns("polygon_cluster_profile_table"))
              )
            ),
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
