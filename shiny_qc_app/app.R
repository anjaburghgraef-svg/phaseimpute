library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(plotly)
library(jsonlite)
library(dplyr)
library(tidyr)

# Helper function to find files flexibly in the output directory
find_file <- function(base_dir, patterns, subdirs = c(".", "*/")) {
  for (subdir in subdirs) {
    for (pattern in patterns) {
      search_path <- file.path(base_dir, subdir, pattern)
      files <- Sys.glob(search_path)
      if (length(files) > 0) {
        # Return the most recent file if multiple found
        return(files[which.max(file.info(files)$mtime)])
      }
    }
  }
  return(NULL)
}

# Load MultiQC data (main data source)
load_multiqc_data <- function(base_dir) {
  multiqc_file <- find_file(base_dir, "multiqc_data/multiqc_data.json",
                            c("multiqc/", "*/multiqc/", "."))
  if (is.null(multiqc_file) || !file.exists(multiqc_file)) return(NULL)

  tryCatch({
    data <- fromJSON(multiqc_file)
    return(data)
  }, error = function(e) {
    return(NULL)
  })
}

# Load parameters
load_params <- function(base_dir) {
  params_file <- find_file(base_dir, "params*.json",
                           c("pipeline_info/", "*/pipeline_info/", "."))
  if (is.null(params_file)) return(NULL)

  tryCatch({
    params <- fromJSON(params_file)

    # Auto-detect steps from available directories
    step_dirs <- c("prep_panel", "simulate", "imputation", "validation")
    detected_steps <- c()
    for (step in step_dirs) {
      step_path <- file.path(base_dir, step)
      if (dir.exists(step_path)) {
        detected_steps <- c(detected_steps, step)
      }
    }
    if (length(detected_steps) > 0) {
      params$steps <- paste(detected_steps, collapse = ",")
    }

    # Auto-detect tools from imputation directory
    imputation_path <- file.path(base_dir, "imputation")
    if (dir.exists(imputation_path)) {
      tool_dirs <- list.dirs(imputation_path, full.names = FALSE, recursive = FALSE)
      # Filter out common non-tool directories
      tool_dirs <- tool_dirs[!tool_dirs %in% c("csv", "bcftools", "stats", "logs")]
      if (length(tool_dirs) > 0) {
        params$tools <- paste(tool_dirs, collapse = ",")
      }
    }

    # Remove deprecated parameters
    params$normalize <- NULL
    params$phase <- NULL
    params$compute_freq <- NULL

    return(params)
  }, error = function(e) {
    return(NULL)
  })
}

# Load MAF spectrum
load_maf_spectrum <- function(base_dir) {
  maf_file <- find_file(base_dir, "maf_spectrum*.json", c("maf/", "*/maf/", "."))
  if (is.null(maf_file)) return(NULL)

  tryCatch({
    json_data <- fromJSON(maf_file)
    maf_data <- json_data$data

    if (is.list(maf_data)) {
      tool_names <- names(maf_data)
      df_list <- lapply(tool_names, function(tool) {
        tool_data <- maf_data[[tool]]
        if (is.matrix(tool_data) || is.data.frame(tool_data)) {
          df <- as.data.frame(tool_data)
          names(df) <- c("MAF", "Dosage_r2")
          df$Tool <- ifelse(tool == "unknown", "Imputation", tool)
          return(df)
        }
        return(NULL)
      })
      df <- do.call(rbind, df_list[!sapply(df_list, is.null)])
      return(df)
    }
    return(NULL)
  }, error = function(e) {
    return(NULL)
  })
}

# Load tool accuracy
load_tool_accuracy <- function(base_dir) {
  accuracy_file <- find_file(base_dir, "tool_accuracy*.txt", c("tool/", "*/tool/", "."))
  if (is.null(accuracy_file)) return(NULL)

  tryCatch({
    data <- read.table(accuracy_file, header = TRUE, sep = "\t", skip = 4, comment.char = "")
    return(data)
  }, error = function(e) {
    return(NULL)
  })
}

# Load GLIMPSE accuracy from MultiQC
load_glimpse_accuracy <- function(base_dir) {
  # Try to find multiqc_glimpse_err_spl.txt
  glimpse_file <- find_file(base_dir, "multiqc_glimpse_err_spl.txt",
                            c("multiqc/multiqc_data/", "*/multiqc/multiqc_data/", "."))
  if (is.null(glimpse_file)) return(NULL)

  tryCatch({
    # Read the file - it has Python dict-like data in GCsS, GCsI, GCsV columns
    data <- read.table(glimpse_file, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, comment.char = "", quote = "")

    # Parse the GCsV column which contains the validation stats
    if ("GCsV" %in% names(data)) {
      # Extract r-squared and concordance metrics from Python dict strings
      parsed_data <- lapply(data$GCsV, function(dict_str) {
        tryCatch({
          # Convert Python dict to proper JSON - just replace single quotes with double quotes
          json_str <- gsub("'", '"', dict_str)

          parsed <- fromJSON(json_str)

          # Calculate total variants and matches from validation genotype counts
          n_variants <- as.numeric(parsed$val_gt_RR) +
                       as.numeric(parsed$val_gt_RA) +
                       as.numeric(parsed$val_gt_AA)

          # Calculate general concordance (overall match rate)
          total_matches <- as.numeric(parsed$RR_hom_matches) +
                          as.numeric(parsed$RA_het_matches) +
                          as.numeric(parsed$AA_hom_matches)
          general_concordance_pct <- ifelse(n_variants > 0,
                                           100 * total_matches / n_variants,
                                           NA)

          # Convert non-ref discordance to concordance for easier comparison
          non_ref_concordance_pct <- 100 - as.numeric(parsed$non_reference_discordance_rate_percent)

          return(list(
            best_gt_rsquared = as.numeric(parsed$best_gt_rsquared),
            imputed_ds_rsquared = as.numeric(parsed$imputed_ds_rsquared),
            non_ref_discordance = as.numeric(parsed$non_reference_discordance_rate_percent),
            non_ref_concordance = non_ref_concordance_pct,
            general_concordance = general_concordance_pct,
            n_variants = n_variants
          ))
        }, error = function(e) {
          return(list(best_gt_rsquared = NA, imputed_ds_rsquared = NA,
                     non_ref_discordance = NA, non_ref_concordance = NA,
                     general_concordance = NA, n_variants = NA))
        })
      })

      # Convert to dataframe - keep full sample name with tool suffix
      result <- data.frame(
        Sample = data$Sample,  # Keep full sample name (e.g., "DUR_10004940.beagle5")
        N_Variants = sapply(parsed_data, function(x) x$n_variants),
        Best_GT_r2 = sapply(parsed_data, function(x) x$best_gt_rsquared),
        Dosage_r2 = sapply(parsed_data, function(x) x$imputed_ds_rsquared),
        General_Concordance_Pct = sapply(parsed_data, function(x) x$general_concordance),
        NonRef_Concordance_Pct = sapply(parsed_data, function(x) x$non_ref_concordance),
        NonRef_Discordance_Pct = sapply(parsed_data, function(x) x$non_ref_discordance),
        stringsAsFactors = FALSE
      )

      # Remove rows with all NA
      result <- result[!is.na(result$Dosage_r2), ]

      return(result)
    }
    return(NULL)
  }, error = function(e) {
    return(NULL)
  })
}

# Load MAF-stratified accuracy from GLIMPSE
load_maf_accuracy <- function(base_dir) {
  maf_file <- find_file(base_dir, "glimpse-err-grp-plot_Imputed_dosage_r-squared_SNPs_.txt",
                        c("multiqc/multiqc_data/", "*/multiqc/multiqc_data/", "."))
  if (is.null(maf_file)) return(NULL)

  tryCatch({
    data <- read.table(maf_file, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, comment.char = "")

    # Parse the data - columns are MAF bins, rows are samples
    # Column names are like "X0", "X1", "X2", etc. representing bins
    # Values are tuples like "(MAF, r2)"
    # Each row represents a sample-tool combination (e.g., DUR_10004940_PREF_PANEL_Tbeagle5)

    # Get bin columns - could be "X0", "X1" or just "0", "1" depending on R version
    bin_cols <- grep("^X?[0-9]+$", names(data), value = TRUE)

    # If no columns found, return NULL
    if (length(bin_cols) == 0) return(NULL)

    # Extract tool from sample names
    result_list <- list()

    for (i in 1:nrow(data)) {
      sample_name <- data$Sample[i]
      # Extract tool from sample name (e.g., "beagle5" or "glimpse2")
      tool <- if (grepl("Tbeagle5$", sample_name)) {
        "beagle5"
      } else if (grepl("Tglimpse2$", sample_name)) {
        "glimpse2"
      } else {
        "unknown"
      }

      # Parse each bin for this sample
      for (col in bin_cols) {
        val <- data[[col]][i]
        if (!is.na(val) && val != "") {
          # Remove parentheses and split
          val <- gsub("[()]", "", val)
          parts <- as.numeric(strsplit(val, ",")[[1]])
          if (length(parts) == 2 && !is.na(parts[1]) && !is.na(parts[2])) {
            result_list[[length(result_list) + 1]] <- data.frame(
              MAF = parts[1],
              Dosage_r2 = parts[2],
              Tool = tool,
              Sample = sample_name
            )
          }
        }
      }
    }

    if (length(result_list) == 0) return(NULL)

    # Combine all results
    all_data <- do.call(rbind, result_list)

    # Average by MAF bin and tool using base R
    all_data$MAF <- as.numeric(all_data$MAF)
    all_data$Dosage_r2 <- as.numeric(all_data$Dosage_r2)

    # Round MAF to consistent precision (3 decimal places) to group properly
    all_data$MAF_rounded <- round(all_data$MAF, 3)

    # Get unique combinations of MAF and Tool
    unique_tools <- unique(all_data$Tool)
    unique_mafs <- unique(all_data$MAF_rounded)

    agg_list <- list()
    for (tool in unique_tools) {
      for (maf in unique_mafs) {
        subset_data <- all_data[all_data$Tool == tool & all_data$MAF_rounded == maf, ]
        if (nrow(subset_data) > 0) {
          agg_list[[length(agg_list) + 1]] <- data.frame(
            MAF = maf,
            Tool = tool,
            Dosage_r2 = mean(subset_data$Dosage_r2, na.rm = TRUE),
            stringsAsFactors = FALSE
          )
        }
      }
    }

    agg_data <- do.call(rbind, agg_list)

    return(agg_data)
  }, error = function(e) {
    return(NULL)
  })
}

# Load project summary
load_project_summary <- function(base_dir) {
  summary_file <- find_file(base_dir, "project_summary*.txt", c("project/", "*/project/", "."))
  if (is.null(summary_file)) return(NULL)

  tryCatch({
    data <- read.table(summary_file, header = TRUE, sep = "\t", skip = 9, comment.char = "")
    return(data)
  }, error = function(e) {
    return(NULL)
  })
}

# Load summary stats (panel samples, input SNPs, etc.)
load_summary_stats <- function(base_dir) {
  stats_file <- find_file(base_dir, "summary_stats.txt", c("qc_stats/", "*/qc_stats/", "."))
  if (is.null(stats_file)) return(NULL)

  tryCatch({
    data <- read.table(stats_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                      colClasses = c("character", "numeric"))

    # Try to get actual panel sample count from REF_PANEL bcftools stats
    panel_stats_file <- find_file(base_dir, "REF_PANEL.panel.bcftools_stats.txt",
                                   c("prep_panel/stats/", "*/prep_panel/stats/", "."))
    if (!is.null(panel_stats_file) && file.exists(panel_stats_file)) {
      # Read the bcftools stats to get actual number of samples in panel VCF
      stats_lines <- readLines(panel_stats_file, n = 50)
      # Look for "number of samples:" line
      sample_line <- grep("number of samples:", stats_lines, value = TRUE)
      if (length(sample_line) > 0) {
        # Extract number from line like "SN	0	number of samples:	15"
        parts <- strsplit(sample_line[1], "\t")[[1]]
        if (length(parts) >= 4) {
          n_samples <- as.numeric(parts[4])
          # Update Panel_Samples in data
          if (!is.na(n_samples) && n_samples > 0) {
            data$Value[data$Metric == "Panel_Samples"] <- n_samples
          }
        }
      }
    }

    return(data)
  }, error = function(e) {
    return(NULL)
  })
}

# Load validation data
load_validation_data <- function(base_dir) {
  val_file <- find_file(base_dir, "AllSamples.txt", c("validation/stats/", "*/validation/stats/", "."))
  if (is.null(val_file)) return(NULL)

  tryCatch({
    data <- read.table(val_file, header = TRUE, sep = "\t", comment.char = "")
    return(data)
  }, error = function(e) {
    return(NULL)
  })
}

# Chromosome accuracy removed - not feasible with current pipeline

# Load INFO score distribution from pipeline-generated file
load_info_scores <- function(base_dir) {
  # Look for info_scores.txt generated by pipeline
  info_file <- find_file(base_dir, "info_scores.txt",
                         c("qc_stats/", "*/qc_stats/", ".", "multiqc/"))

  if (is.null(info_file)) return(NULL)

  tryCatch({
    data <- read.table(info_file, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, comment.char = "")

    # Check if we have the expected columns
    if (!all(c("INFO_Range", "Variant_Count", "Percentage") %in% names(data))) return(NULL)

    return(data)
  }, error = function(e) {
    return(NULL)
  })
}

# Extract panel sample count from input CSV or output
get_panel_sample_count <- function(panel_path, base_dir) {
  # First try the original input CSV
  if (!is.null(panel_path) && file.exists(panel_path)) {
    tryCatch({
      panel_csv <- read.csv(panel_path, stringsAsFactors = FALSE, header = TRUE)
      if (ncol(panel_csv) > 0) {
        unique_samples <- length(unique(panel_csv[, 1]))
        return(unique_samples)
      }
    }, error = function(e) {})
  }

  # Try the output CSV from prep_panel
  output_panel_csv <- find_file(base_dir, "panel.csv",
                                 c("prep_panel/csv/", "*/prep_panel/csv/", "."))
  if (!is.null(output_panel_csv)) {
    tryCatch({
      panel_csv <- read.csv(output_panel_csv, stringsAsFactors = FALSE, header = TRUE)
      # This just shows REF_PANEL per chromosome, not individual samples
      # The original input CSV would have the sample count
      return(NULL)
    }, error = function(e) {})
  }

  return(NULL)
}

# Extract input/target SNP count from stats or CSV
get_input_snp_count <- function(input_path, base_dir) {
  # First try the original input CSV
  if (!is.null(input_path) && file.exists(input_path)) {
    tryCatch({
      input_csv <- read.csv(input_path, stringsAsFactors = FALSE, header = TRUE)
      # Count unique samples
      if (ncol(input_csv) > 0) {
        n_samples <- length(unique(input_csv[, 1]))
        # We have sample count but not SNP count from CSV alone
      }
    }, error = function(e) {})
  }

  # Try to find bcftools stats for input/target files
  # These aren't typically in the output, so return NULL
  return(NULL)
}

# UI
ui <- fluidPage(
  title = "QC Report for nf-core/phaseimpute",

  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
    tags$script(HTML("
      $(document).ready(function() {
        // Track current tab
        window.currentTab = 'overview';

        // Tab switching function
        window.switchTab = function(tabName) {
          $('.tab-content').hide();
          $('#' + tabName + '-content').show();
          $('.navbar-tab').removeClass('active');
          $('button[data-tab=\"' + tabName + '\"]').addClass('active');
          window.currentTab = tabName;

          // Trigger update of sidebar navigation
          Shiny.setInputValue('active_tab', tabName, {priority: 'event'});
        };

        // Smooth scroll to section
        window.scrollToSection = function(sectionId) {
          $('html, body').animate({
            scrollTop: $('#' + sectionId).offset().top - 90
          }, 500);
        };

        // Initialize - show overview by default
        switchTab('overview');
      });
    "))
  ),

  # Fixed top navbar
  tags$div(
    class = "top-navbar",
    tags$img(src = "logo_hendrix.png", class = "navbar-logo"),
    tags$span("phaseimpute QC Report", class = "navbar-title"),
    tags$div(
      class = "navbar-tabs",
      tags$button("Overview", class = "navbar-tab active",
                 `data-tab` = "overview",
                 onclick = "switchTab('overview')"),
      tags$button("Quality & Accuracy", class = "navbar-tab",
                 `data-tab` = "quality",
                 onclick = "switchTab('quality')")
  ),
    tags$div(
      class = "navbar-links",
      tags$a(href = "https://nf-co.re/phaseimpute/latest", target = "_blank",
            class = "navbar-link", "Documentation"),
      tags$span(style = "color: white; margin: 0 5px;", "|"),
      tags$a(href = "https://github.com/nf-core/phaseimpute/tree/1.1.0", target = "_blank",
            class = "navbar-link", "Original Pipeline"),
      tags$span(style = "color: white; margin: 0 5px;", "|"),
      tags$a(href = "https://github.com/anjaburghgraef-svg/phaseimpute", target = "_blank",
            class = "navbar-link", "Fork Repository")
    )
  ),

  # Fixed left sidebar
  tags$div(
    class = "left-sidebar",

    # Data loading section
    tags$div(
      class = "sidebar-section",
      tags$h4("Data Loading"),
      textInput("data_dir", NULL,
               value = "your_data/output",
               placeholder = "Results directory path"),
      actionButton("load_data", "Load Data", icon = icon("sync"),
                  class = "btn-primary", style = "width: 100%;"),
      tags$small(style = "color: #666; display: block; margin-top: 10px;",
                "Point to the top-level results directory from the pipeline.")
    ),

    # Navigation links - dynamically updated based on active tab
    tags$div(
      class = "sidebar-section",
      tags$h4("Navigation"),
      uiOutput("sidebar_nav")
    )
  ),

  # Main content area
  tags$div(
    class = "main-content",

    # Overview tab content
    tags$div(
      id = "overview-content",
      class = "tab-content",

      tags$div(
        class = "section",
        id = "overview-metrics",
        tags$h3("Overview Metrics"),
        uiOutput("overview_value_boxes")
      ),

      tags$div(
        class = "section",
        id = "pipeline-params",
        tags$h3("Pipeline Parameters"),
        tags$div(
          class = "content-box",
          DTOutput("params_table")
        )
      ),

      tags$div(
        class = "section",
        id = "sample-overview",
        tags$h3("Sample Overview"),
        tags$div(
          class = "content-box",
          htmlOutput("sample_overview_text"),
          br(),
          actionButton("show_full_sample_list", "Show Full Sample List",
                      icon = icon("list"))
        )
      )
    ),

    # Quality & Accuracy tab content
    tags$div(
      id = "quality-content",
      class = "tab-content",
      style = "display: none;",

      # Top row: MAF plot (left) + Accuracy metrics (right)
      tags$div(
        style = "display: flex; gap: 10px; margin-bottom: 10px;",

        # Left column - MAF plot
        tags$div(
          style = "flex: 1;",
          tags$div(
            class = "section",
            id = "maf-plot",
            tags$h3("Imputation Quality by MAF"),
            tags$div(
              class = "content-box",
              plotlyOutput("maf_spectrum_plot", height = "550px")
            )
          )
        ),

        # Right column - Accuracy metrics stacked
        tags$div(
          style = "flex: 1; display: flex; flex-direction: column; gap: 0px;",

          tags$div(
            class = "section",
            id = "overall-accuracy",
            style = "margin-bottom: 0;",
            tags$h3("Overall Accuracy Metrics"),
            tags$div(
              class = "content-box",
              uiOutput("tool_accuracy_display")
            )
          ),

          tags$div(
            class = "section",
            id = "quality-summary",
            style = "margin-bottom: 0;",
            tags$h3("Accuracy Plots"),
            tags$div(
              class = "content-box",
              plotlyOutput("quality_summary_plot", height = "280px")
            )
          )
        )
      ),

      # Bottom row: Per-sample accuracy (full width)
      tags$div(
        class = "section",
        id = "per-sample-accuracy",
        tags$h3("Per-Sample Accuracy"),
        tags$div(
          class = "content-box",
          uiOutput("sample_accuracy_display")
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {

  # Reactive values to store loaded data
  data_store <- reactiveValues(
    multiqc_data = NULL,
    params = NULL,
    project_summary = NULL,
    tool_accuracy = NULL,
    maf_spectrum = NULL,
    validation_data = NULL,
    tool_used = NULL,
    base_dir = NULL,
    summary_stats = NULL,
    info_scores = NULL
  )

  # Load data when button is clicked
  observeEvent(input$load_data, {
    base_dir <- input$data_dir

    if (!dir.exists(base_dir)) {
      showNotification("Directory not found!", type = "error", duration = 5)
      return()
    }

    # Store base_dir for use in value boxes
    data_store$base_dir <- base_dir

    withProgress(message = 'Loading data...', value = 0, {
      # Load MultiQC data (primary source)
      incProgress(0.2, detail = "Loading MultiQC data")
      data_store$multiqc_data <- load_multiqc_data(base_dir)

      # Load parameters
      incProgress(0.2, detail = "Loading parameters")
      data_store$params <- load_params(base_dir)

      # Detect tools used - can be multiple
      # First try from params
      if (!is.null(data_store$params)) {
        tools_param <- data_store$params$tools
        # If it's a string with comma-separated tools, split it
        if (is.character(tools_param) && grepl(",", tools_param)) {
          data_store$tool_used <- strsplit(tools_param, ",\\s*")[[1]]
        } else {
          data_store$tool_used <- tools_param
        }
      }

      # Also detect from actual data (in case multiple tools were run)
      if (!is.null(data_store$multiqc_data)) {
        tryCatch({
          bcftools_data <- data_store$multiqc_data$report_general_stats_data$bcftools
          samples <- names(bcftools_data)
          # Extract unique tools from sample names
          tools_found <- unique(gsub(".*\\.(beagle5|minimac4|glimpse1|glimpse2|quilt|stitch)$", "\\1",
                                     samples[grepl("\\.(beagle5|minimac4|glimpse1|glimpse2|quilt|stitch)$", samples)]))
          if (length(tools_found) > 0) {
            data_store$tool_used <- tools_found
          }
        }, error = function(e) {})
      }

      # Load project summary
      incProgress(0.2, detail = "Loading project summary")
      data_store$project_summary <- load_project_summary(base_dir)

      # Load summary stats (panel/input counts)
      data_store$summary_stats <- load_summary_stats(base_dir)

      # Load MAF spectrum
      incProgress(0.2, detail = "Loading MAF spectrum")
      data_store$maf_spectrum <- load_maf_spectrum(base_dir)

      # Load tool accuracy
      incProgress(0.1, detail = "Loading accuracy metrics")
      data_store$tool_accuracy <- load_tool_accuracy(base_dir)

      # Load GLIMPSE accuracy if tool_accuracy not found
      if (is.null(data_store$tool_accuracy)) {
        data_store$tool_accuracy <- load_glimpse_accuracy(base_dir)
      }

      # Load validation data
      incProgress(0.05, detail = "Loading validation data")
      data_store$validation_data <- load_validation_data(base_dir)

      # Load MAF-stratified accuracy (replaces MAF spectrum if not found)
      incProgress(0.05, detail = "Loading MAF accuracy")
      if (is.null(data_store$maf_spectrum)) {
        data_store$maf_spectrum <- load_maf_accuracy(base_dir)
      }

      # Chromosome accuracy and INFO scores removed - not feasible/useful
    })

    # Show diagnostic message
    msg_parts <- c()
    if (!is.null(data_store$multiqc_data)) msg_parts <- c(msg_parts, "MultiQC")
    if (!is.null(data_store$params)) msg_parts <- c(msg_parts, "Params")
    if (!is.null(data_store$summary_stats)) msg_parts <- c(msg_parts, "Summary Stats")
    if (!is.null(data_store$tool_accuracy)) msg_parts <- c(msg_parts, paste0("Accuracy (", nrow(data_store$tool_accuracy), " rows)"))
    if (!is.null(data_store$maf_spectrum)) msg_parts <- c(msg_parts, paste0("MAF (", nrow(data_store$maf_spectrum), " rows)"))

    notification_msg <- if (length(msg_parts) > 0) {
      paste("Loaded:", paste(msg_parts, collapse = ", "))
    } else {
      "Data loaded (no files found)"
    }

    showNotification(notification_msg, type = "message", duration = 5)
  })

  # Extract sample stats from MultiQC
  sample_stats <- reactive({
    if (is.null(data_store$multiqc_data)) return(NULL)

    tryCatch({
      bcftools_data <- data_store$multiqc_data$report_general_stats_data$bcftools
      if (is.null(bcftools_data)) return(NULL)

      samples <- names(bcftools_data)

      # Filter to imputed samples only (exclude .truth and REF_PANEL)
      imputed_samples <- samples[!grepl("\\.truth$", samples) & samples != "REF_PANEL"]
      bcftools_imputed <- bcftools_data[imputed_samples]

      df <- data.frame(
        Sample = gsub("\\.(beagle5|minimac4|glimpse1|glimpse2|quilt|stitch)$", "", imputed_samples),  # Remove tool suffix
        Tool = gsub(".*\\.", "", imputed_samples),  # Extract tool name
        SNPs = sapply(bcftools_imputed, function(x) x$number_of_SNPs),
        Records = sapply(bcftools_imputed, function(x) x$number_of_records),
        Ts_Tv = sapply(bcftools_imputed, function(x) round(x$tstv, 2)),
        stringsAsFactors = FALSE
      )
      return(df)
    }, error = function(e) {
      return(NULL)
    })
  })

  # Helper function to format numbers in K
  format_k <- function(num) {
    if (is.na(num) || num == "NA") return("N/A")
    num <- as.numeric(num)
    if (is.na(num)) return("N/A")
    paste0(round(num / 1000, 1), "K")
  }

  # Sidebar navigation - updates based on active tab
  output$sidebar_nav <- renderUI({
    # Get active tab from input (defaults to overview)
    active_tab <- input$active_tab
    if (is.null(active_tab)) active_tab <- "overview"

    if (active_tab == "overview") {
      tagList(
        tags$a(href = "#overview-metrics", class = "sidebar-nav-link",
              onclick = "scrollToSection('overview-metrics'); return false;",
              "Overview Metrics"),
        tags$a(href = "#pipeline-params", class = "sidebar-nav-link",
              onclick = "scrollToSection('pipeline-params'); return false;",
              "Pipeline Parameters"),
        tags$a(href = "#sample-overview", class = "sidebar-nav-link",
              onclick = "scrollToSection('sample-overview'); return false;",
              "Sample Overview")
      )
    } else {
      # Quality & Accuracy tab
      tagList(
        tags$a(href = "#maf-plot", class = "sidebar-nav-link",
              onclick = "scrollToSection('maf-plot'); return false;",
              "MAF Quality Plot"),
        tags$a(href = "#overall-accuracy", class = "sidebar-nav-link",
              onclick = "scrollToSection('overall-accuracy'); return false;",
              "Overall Accuracy"),
        tags$a(href = "#quality-summary", class = "sidebar-nav-link",
              onclick = "scrollToSection('quality-summary'); return false;",
              "Accuracy Plots"),
        tags$a(href = "#per-sample-accuracy", class = "sidebar-nav-link",
              onclick = "scrollToSection('per-sample-accuracy'); return false;",
              "Per-Sample Accuracy")
      )
    }
  })

  # Overview value boxes
  output$overview_value_boxes <- renderUI({
    target_text <- if (!is.null(data_store$summary_stats)) {
      target_samples <- data_store$summary_stats$Value[data_store$summary_stats$Metric == "Target_Samples"]
      target_snps <- data_store$summary_stats$Value[data_store$summary_stats$Metric == "Target_SNPs"]
      if (length(target_samples) > 0 && !is.na(target_samples) && target_samples != "NA" &&
          length(target_snps) > 0 && !is.na(target_snps) && target_snps != "NA") {
        paste0(target_samples, " samples\n", format_k(target_snps), " SNPs")
      } else { "N/A" }
    } else { "N/A" }

    truth_text <- if (!is.null(data_store$summary_stats)) {
      truth_samples <- data_store$summary_stats$Value[data_store$summary_stats$Metric == "Truth_Samples"]
      truth_snps <- data_store$summary_stats$Value[data_store$summary_stats$Metric == "Truth_SNPs"]
      if (length(truth_samples) > 0 && !is.na(truth_samples) && truth_samples != "NA" &&
          length(truth_snps) > 0 && !is.na(truth_snps) && truth_snps != "NA") {
        paste0(truth_samples, " samples\n", format_k(truth_snps), " SNPs")
      } else { "N/A" }
    } else { "N/A" }

    panel_text <- if (!is.null(data_store$summary_stats)) {
      panel_samples <- data_store$summary_stats$Value[data_store$summary_stats$Metric == "Panel_Samples"]
      panel_snps <- data_store$summary_stats$Value[data_store$summary_stats$Metric == "Panel_SNPs"]
      if (length(panel_samples) > 0 && !is.na(panel_samples) && panel_samples != "NA" &&
          length(panel_snps) > 0 && !is.na(panel_snps) && panel_snps != "NA") {
        paste0(panel_samples, " samples\n", format_k(panel_snps), " SNPs")
      } else { "N/A" }
    } else { "N/A" }

    imputation_text <- if (!is.null(data_store$summary_stats)) {
      target_snps <- data_store$summary_stats$Value[data_store$summary_stats$Metric == "Target_SNPs"]
      panel_snps <- data_store$summary_stats$Value[data_store$summary_stats$Metric == "Panel_SNPs"]
      if (length(target_snps) > 0 && !is.na(target_snps) && target_snps != "NA" &&
          length(panel_snps) > 0 && !is.na(panel_snps) && panel_snps != "NA") {
        paste0(format_k(target_snps), " → ", format_k(panel_snps))
      } else { "N/A" }
    } else { "N/A" }

    tool_text <- if (!is.null(data_store$tool_used)) {
      if (length(data_store$tool_used) > 1) {
        paste(toupper(data_store$tool_used), collapse = " + ")
      } else {
        toupper(data_store$tool_used)
      }
    } else { "N/A" }

    tags$div(
      class = "value-boxes",
      tags$div(
        class = "value-box navy",
        tags$div(class = "value-box-value", target_text),
        tags$div(class = "value-box-subtitle", "Target Data")
      ),
      tags$div(
        class = "value-box olive",
        tags$div(class = "value-box-value", truth_text),
        tags$div(class = "value-box-subtitle", "Validation Data")
      ),
      tags$div(
        class = "value-box light-blue",
        tags$div(class = "value-box-value", panel_text),
        tags$div(class = "value-box-subtitle", "Reference Panel")
      ),
      tags$div(
        class = "value-box navy",
        tags$div(class = "value-box-value", imputation_text),
        tags$div(class = "value-box-subtitle", "Imputation")
      ),
      tags$div(
        class = "value-box light-blue",
        tags$div(class = "value-box-value", tool_text),
        tags$div(class = "value-box-subtitle", "Imputation Tool(s)")
      )
    )
  })

  # Parameters table
  output$params_table <- renderDT({
    if (is.null(data_store$params)) {
      return(data.frame(Info = "No parameters loaded"))
    }

    key_params <- c("tools", "steps", "input", "panel", "input_truth", "fasta", "outdir")

    param_df <- data.frame(
      Parameter = character(),
      Value = character(),
      stringsAsFactors = FALSE
    )

    for (param in key_params) {
      if (!is.null(data_store$params[[param]])) {
        value <- data_store$params[[param]]
        if (is.logical(value)) {
          value <- ifelse(value, "Yes", "No")
        } else if (is.list(value)) {
          value <- paste(unlist(value), collapse = ", ")
        }
        param_df <- rbind(param_df, data.frame(
          Parameter = param,
          Value = as.character(value)
        ))
      }
    }

    datatable(param_df,
              options = list(pageLength = 10, dom = 't', scrollX = TRUE),
              rownames = FALSE)
  })

  # Project summary display - generate from available data

  output$sample_overview_text <- renderUI({
    if (is.null(data_store$summary_stats)) {
      return(tags$p("N/A - No sample data available"))
    }

    # Get values from summary stats
    target_samples <- data_store$summary_stats$Value[data_store$summary_stats$Metric == "Target_Samples"]
    truth_samples <- data_store$summary_stats$Value[data_store$summary_stats$Metric == "Truth_Samples"]
    panel_samples <- data_store$summary_stats$Value[data_store$summary_stats$Metric == "Panel_Samples"]

    # Convert to numeric, handling NA
    n_target <- if (length(target_samples) > 0 && target_samples != "NA") as.numeric(target_samples) else NA
    n_truth <- if (length(truth_samples) > 0 && truth_samples != "NA") as.numeric(truth_samples) else NA
    n_panel <- if (length(panel_samples) > 0 && panel_samples != "NA") as.numeric(panel_samples) else NA

    # Determine overlap from validation data if available
    n_both <- if (!is.null(data_store$tool_accuracy) && "Sample" %in% names(data_store$tool_accuracy)) {
      # Count unique samples that have validation results
      length(unique(gsub("\\.(beagle5|minimac4|glimpse1|glimpse2|quilt|stitch)$", "",
                        data_store$tool_accuracy$Sample)))
    } else { NA }

    tagList(
      tags$p(tags$strong("Target samples (imputed):"), ifelse(!is.na(n_target), n_target, "N/A")),
      tags$p(tags$strong("Validation samples (truth):"), ifelse(!is.na(n_truth), n_truth, "N/A")),
      tags$p(tags$strong("Reference panel samples:"), ifelse(!is.na(n_panel), n_panel, "N/A")),
      tags$p(tags$strong("Samples with validation:"), ifelse(!is.na(n_both), n_both, "N/A")),
      tags$hr(),
      tags$p(tags$em(style = "color: #666; font-size: 90%;",
                    "Samples with validation = target samples that have corresponding truth data for accuracy assessment"))
    )
  })


  # Build comprehensive sample list from all data sources
  build_sample_list <- reactive({
    sample_list <- data.frame(
      Sample = character(),
      In_Target = character(),
      In_Truth = character(),
      In_Panel = character(),
      stringsAsFactors = FALSE
    )

    # Get samples from MultiQC bcftools data
    if (!is.null(data_store$multiqc_data)) {
      tryCatch({
        bcftools_data <- data_store$multiqc_data$report_general_stats_data$bcftools
        if (!is.null(bcftools_data)) {
          all_samples <- names(bcftools_data)

          # Extract base sample names and categorize
          for (sample_name in all_samples) {
            # Skip REF_PANEL
            if (sample_name == "REF_PANEL") next

            # Remove tool suffixes to get base sample name
            base_name <- gsub("\\.(beagle5|minimac4|glimpse1|glimpse2|quilt|stitch|truth)$", "", sample_name)

            # Check if it's a truth sample or target sample
            is_truth <- grepl("\\.truth$", sample_name)
            is_target <- grepl("\\.(beagle5|minimac4|glimpse1|glimpse2|quilt|stitch)$", sample_name)

            # Check if sample already in list
            if (base_name %in% sample_list$Sample) {
              # Update existing entry
              idx <- which(sample_list$Sample == base_name)
              if (is_target && sample_list$In_Target[idx] == "No") {
                sample_list$In_Target[idx] <- "Yes"
              }
              if (is_truth && sample_list$In_Truth[idx] == "No") {
                sample_list$In_Truth[idx] <- "Yes"
              }
            } else {
              # Add new entry
              sample_list <- rbind(sample_list, data.frame(
                Sample = base_name,
                In_Target = ifelse(is_target, "Yes", "No"),
                In_Truth = ifelse(is_truth, "Yes", "No"),
                In_Panel = "No",
                stringsAsFactors = FALSE
              ))
            }
          }
        }
      }, error = function(e) {})
    }

    # Get panel samples from panel VCF files
    if (!is.null(data_store$base_dir)) {
      tryCatch({
        # Look for panel CSV
        panel_csv <- find_file(data_store$base_dir, "panel.csv",
                               c("prep_panel/csv/", "*/prep_panel/csv/", "."))

        if (!is.null(panel_csv) && file.exists(panel_csv)) {
          csv_data <- read.csv(panel_csv, stringsAsFactors = FALSE, header = TRUE)
          if (nrow(csv_data) > 0) {
            # Get first VCF file from the CSV
            vcf_path <- csv_data$vcf[1]

            if (file.exists(vcf_path)) {
              # Use bcftools to query sample names
              sample_cmd <- paste0("bcftools query -l ", shQuote(vcf_path))
              panel_samples <- system(sample_cmd, intern = TRUE)

              # Add panel samples to list
              for (panel_sample in panel_samples) {
                if (panel_sample %in% sample_list$Sample) {
                  # Update existing entry
                  idx <- which(sample_list$Sample == panel_sample)
                  sample_list$In_Panel[idx] <- "Yes"
                } else {
                  # Add new entry
                  sample_list <- rbind(sample_list, data.frame(
                    Sample = panel_sample,
                    In_Target = "No",
                    In_Truth = "No",
                    In_Panel = "Yes",
                    stringsAsFactors = FALSE
                  ))
                }
              }
            }
          }
        }
      }, error = function(e) {
        # Silently fail if panel samples can't be extracted
      })
    }

    # Sort by sample name
    if (nrow(sample_list) > 0) {
      sample_list <- sample_list[order(sample_list$Sample), ]
    }

    return(sample_list)
  })

  # Show full sample list in modal
  observeEvent(input$show_full_sample_list, {
    sample_list <- build_sample_list()

    if (!is.null(sample_list) && nrow(sample_list) > 0) {
      showModal(modalDialog(
        title = "Complete Sample List",
        size = "l",
        DTOutput("full_sample_table"),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    } else {
      showModal(modalDialog(
        title = "Sample List",
        "No sample data available",
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    }
  })

  output$full_sample_table <- renderDT({
    sample_list <- build_sample_list()
    if (is.null(sample_list) || nrow(sample_list) == 0) return(NULL)

    datatable(sample_list,
              options = list(pageLength = 25, scrollX = TRUE),
              rownames = FALSE,
              filter = "top",
              colnames = c("Sample ID", "In Target", "In Validation", "In Panel")) %>%
      formatStyle(
        'In_Target',
        backgroundColor = styleEqual(c('Yes', 'No'), c('#d4edda', '#f8d7da'))
      ) %>%
      formatStyle(
        'In_Truth',
        backgroundColor = styleEqual(c('Yes', 'No'), c('#d4edda', '#f8d7da'))
      ) %>%
      formatStyle(
        'In_Panel',
        backgroundColor = styleEqual(c('Yes', 'No'), c('#d4edda', '#f8d7da'))
      )
  })

  # MAF spectrum plot
  output$maf_spectrum_plot <- renderPlotly({
    if (is.null(data_store$maf_spectrum)) {
      return(plotly_empty() %>%
               layout(title = list(text = "N/A - No MAF spectrum data available",
                                 font = list(color = "#999", size = 14))))
    }

    tryCatch({
      # Ensure data is a proper data frame
      maf_data <- as.data.frame(data_store$maf_spectrum)

      # Check if required columns exist
      if (!all(c("MAF", "Dosage_r2", "Tool") %in% names(maf_data))) {
        return(plotly_empty() %>%
                 layout(title = list(text = paste("Error: Missing columns. Found:", paste(names(maf_data), collapse = ", ")),
                                   font = list(color = "#999", size = 14))))
      }

      # Ensure proper data types
      maf_data$MAF <- as.numeric(maf_data$MAF)
      maf_data$Dosage_r2 <- as.numeric(maf_data$Dosage_r2)
      maf_data$Tool <- as.character(maf_data$Tool)

      # Remove any rows with NA values
      maf_data <- maf_data[complete.cases(maf_data), ]

      if (nrow(maf_data) == 0) {
        return(plotly_empty() %>%
                 layout(title = list(text = "No valid MAF data after filtering",
                                   font = list(color = "#999", size = 14))))
      }

      # Sort by MAF for proper line plotting
      maf_data <- maf_data[order(maf_data$Tool, maf_data$MAF), ]

      # Create plot directly with plotly instead of ggplot
      tools <- unique(maf_data$Tool)
      # Use distinct colors for tools
      tool_colors <- c(
        "beagle5" = "#2091ae",    # Light blue
        "glimpse2" = "#9cbe2b",   # Green
        "glimpse1" = "#16315b",   # Navy
        "minimac4" = "#e76f51",   # Orange
        "quilt" = "#f4a261",      # Light orange
        "stitch" = "#8338ec"      # Purple
      )

      p <- plot_ly()

      for (i in seq_along(tools)) {
        tool_data <- maf_data[maf_data$Tool == tools[i], ]

        # Sort by MAF to ensure proper line connection
        tool_data <- tool_data[order(tool_data$MAF), ]

        # Remove duplicates - keep only one point per MAF value (shouldn't happen but safety check)
        tool_data <- tool_data[!duplicated(tool_data$MAF), ]

        tool_color <- if (tools[i] %in% names(tool_colors)) tool_colors[[tools[i]]] else "#666666"

        p <- p %>%
          add_trace(
            data = tool_data,
            x = ~MAF,
            y = ~Dosage_r2,
            name = tools[i],
            type = 'scatter',
            mode = 'lines+markers',
            line = list(color = tool_color, width = 2, shape = 'linear'),
            marker = list(
              color = tool_color,
              size = 8,
              line = list(color = 'white', width = 1.5)
            ),
            hovertemplate = paste0(
              '<b>', tools[i], '</b><br>',
              'MAF: %{x:.3f}<br>',
              'Dosage r²: %{y:.3f}<br>',
              '<extra></extra>'
            )
          )
      }

      p %>%
        layout(
          title = list(text = "Imputation Quality by MAF Bin", font = list(size = 16)),
          xaxis = list(
            title = "Minor Allele Frequency (MAF)",
            gridcolor = '#e8e8e8',
            showgrid = TRUE
          ),
          yaxis = list(
            title = "Mean Dosage r²",
            range = c(0, 1),
            gridcolor = '#e8e8e8',
            showgrid = TRUE
          ),
          hovermode = "x unified",
          legend = list(title = list(text = "Tool")),
          plot_bgcolor = 'white',
          paper_bgcolor = 'white'
        )
    }, error = function(e) {
      plotly_empty() %>%
        layout(title = list(text = paste("Error creating plot:", e$message),
                          font = list(color = "#999", size = 14)))
    })
  })

  # Tool accuracy display
  output$tool_accuracy_display <- renderUI({
    if (is.null(data_store$tool_accuracy)) {
      return(div(class = "no-data-message",
                icon("info-circle", class = "fa-3x"),
                br(), br(),
                tags$b("N/A - No validation data available"),
                br(),
                tags$small("Run pipeline with --input_truth to enable validation")))
    }

    # Show data structure info if in debug mode
    if (nrow(data_store$tool_accuracy) == 0) {
      return(div(class = "no-data-message",
                tags$b("Tool accuracy data loaded but empty")))
    }

    DTOutput("tool_accuracy_table")
  })

  output$tool_accuracy_table <- renderDT({
    if (is.null(data_store$tool_accuracy)) return(NULL)

    tryCatch({
      # Calculate summary statistics if we have per-sample data
      if ("Sample" %in% names(data_store$tool_accuracy)) {
        # This is GLIMPSE per-sample data, summarize it by tool
        # Extract tool from sample name
        tool_data <- as.data.frame(data_store$tool_accuracy)

        # Check if we have the expected columns (N_Variants is optional)
        required_cols <- c("Sample", "Dosage_r2", "Best_GT_r2", "NonRef_Discordance_Pct")
        missing_cols <- setdiff(required_cols, names(tool_data))

        if (length(missing_cols) > 0) {
          return(data.frame(
            Error = paste("Missing columns:", paste(missing_cols, collapse = ", "),
                         ". Available:", paste(names(tool_data), collapse = ", "))
          ))
        }

        # Extract tool name from sample
        tool_data$Tool <- sapply(tool_data$Sample, function(s) {
          if (grepl("\\.beagle5$", s)) return("beagle5")
          if (grepl("\\.glimpse2$", s)) return("glimpse2")
          if (grepl("\\.glimpse1$", s)) return("glimpse1")
          if (grepl("\\.minimac4$", s)) return("minimac4")
          if (grepl("\\.quilt$", s)) return("quilt")
          if (grepl("\\.stitch$", s)) return("stitch")
          return("unknown")
        })

        # Convert metrics to numeric
        tool_data$Dosage_r2 <- as.numeric(tool_data$Dosage_r2)
        tool_data$Best_GT_r2 <- as.numeric(tool_data$Best_GT_r2)
        tool_data$NonRef_Discordance_Pct <- as.numeric(tool_data$NonRef_Discordance_Pct)

        # Aggregate by tool using base R
        tools <- unique(tool_data$Tool)
        summary_list <- lapply(tools, function(t) {
          tool_subset <- tool_data[tool_data$Tool == t, ]
          data.frame(
            Tool = t,
            N_Samples = nrow(tool_subset),
            Mean_Dosage_r2 = round(mean(tool_subset$Dosage_r2, na.rm = TRUE), 4),
            Mean_Best_GT_r2 = round(mean(tool_subset$Best_GT_r2, na.rm = TRUE), 4),
            Mean_General_Concordance_Pct = round(mean(tool_subset$General_Concordance_Pct, na.rm = TRUE), 2),
            Mean_NonRef_Concordance_Pct = round(mean(tool_subset$NonRef_Concordance_Pct, na.rm = TRUE), 2),
            stringsAsFactors = FALSE
          )
        })
        summary_data <- do.call(rbind, summary_list)

        datatable(as.data.frame(summary_data),
                  options = list(pageLength = 10, dom = 't', scrollX = TRUE),
                  rownames = FALSE,
                  colnames = c("Tool", "N Samples", "Mean Dosage r²", "Mean Best GT r²",
                              "Mean General Concordance %", "Mean Non-Ref Concordance %"))
      } else {
        # Original format
        datatable(as.data.frame(data_store$tool_accuracy),
                  options = list(pageLength = 10, dom = 't', scrollX = TRUE),
                  rownames = FALSE) %>%
          formatRound(columns = grep("r2|concordance|pct", names(data_store$tool_accuracy), ignore.case = TRUE),
                      digits = 4)
      }
    }, error = function(e) {
      return(data.frame(Error = paste("Could not parse accuracy data:", e$message)))
    })
  })

  # Quality summary plot
  output$quality_summary_plot <- renderPlotly({
    if (is.null(data_store$tool_accuracy)) {
      return(plotly_empty() %>%
               layout(title = list(text = "N/A - No validation data available",
                                 font = list(color = "#999", size = 14))))
    }

    tryCatch({
      # Handle both formats
      if ("Sample" %in% names(data_store$tool_accuracy)) {
        # GLIMPSE format - calculate means by tool
        tool_data <- as.data.frame(data_store$tool_accuracy)

        # Extract tool name from sample
        tool_data$Tool <- sapply(tool_data$Sample, function(s) {
          if (grepl("\\.beagle5$", s)) return("beagle5")
          if (grepl("\\.glimpse2$", s)) return("glimpse2")
          if (grepl("\\.glimpse1$", s)) return("glimpse1")
          if (grepl("\\.minimac4$", s)) return("minimac4")
          if (grepl("\\.quilt$", s)) return("quilt")
          if (grepl("\\.stitch$", s)) return("stitch")
          return("unknown")
        })

        # Convert metrics to numeric
        tool_data$Dosage_r2 <- as.numeric(tool_data$Dosage_r2)
        tool_data$Best_GT_r2 <- as.numeric(tool_data$Best_GT_r2)
        tool_data$General_Concordance_Pct <- as.numeric(tool_data$General_Concordance_Pct)
        tool_data$NonRef_Concordance_Pct <- as.numeric(tool_data$NonRef_Concordance_Pct)

        # Aggregate by tool using base R
        tools <- unique(tool_data$Tool)
        summary_list <- lapply(tools, function(t) {
          tool_subset <- tool_data[tool_data$Tool == t, ]
          data.frame(
            Tool = t,
            Dosage_r2 = mean(tool_subset$Dosage_r2, na.rm = TRUE),
            Best_GT_r2 = mean(tool_subset$Best_GT_r2, na.rm = TRUE),
            General_Concordance = mean(tool_subset$General_Concordance_Pct, na.rm = TRUE),
            NonRef_Concordance = mean(tool_subset$NonRef_Concordance_Pct, na.rm = TRUE),
            stringsAsFactors = FALSE
          )
        })
        summary_data <- do.call(rbind, summary_list)

        # Reshape to long format for plotting using base R
        # Convert r² values to percentages (multiply by 100)
        metrics <- data.frame(
          Tool = rep(summary_data$Tool, 4),
          Metric = rep(c("Dosage r²", "Best GT r²", "General Concordance", "Non-Ref Concordance"), each = nrow(summary_data)),
          Value = c(summary_data$Dosage_r2 * 100, summary_data$Best_GT_r2 * 100, summary_data$General_Concordance, summary_data$NonRef_Concordance),
          stringsAsFactors = FALSE
        )

        # Ensure numeric values
        metrics$Value <- as.numeric(metrics$Value)

        # Create plot directly with plotly
        tools <- unique(metrics$Tool)
        metric_names <- unique(metrics$Metric)
        colors <- RColorBrewer::brewer.pal(max(3, length(tools)), "Set2")[1:length(tools)]

        p <- plot_ly()

        for (i in seq_along(tools)) {
          tool_data <- metrics[metrics$Tool == tools[i], ]
          p <- p %>%
            add_trace(
              data = tool_data,
              x = ~Metric,
              y = ~Value,
              name = tools[i],
              type = 'bar',
              marker = list(color = colors[i])
            )
        }

        p <- p %>%
          layout(
            xaxis = list(
              title = "",
              tickangle = 0,
              tickmode = "array",
              tickvals = 0:3,
              ticktext = c("Dosage r²", "Best GT r²", "General<br>Concordance", "Non-Ref<br>Concordance")
            ),
            yaxis = list(title = "Value (%)", range = c(0, 100)),
            barmode = 'group',
            legend = list(title = list(text = "Tool")),
            margin = list(t = 10, b = 70)
          )

      } else {
        # Original format
        metrics_data <- data_store$tool_accuracy[1, ]
        metric_names <- names(metrics_data)[grepl("Mean_", names(metrics_data))]

        metrics <- data.frame(
          Metric = gsub("Mean_", "", gsub("_", " ", metric_names)),
          Value = as.numeric(metrics_data[metric_names]),
          stringsAsFactors = FALSE
        )

        colors <- RColorBrewer::brewer.pal(max(3, nrow(metrics)), "Set2")[1:nrow(metrics)]

        p <- plot_ly(
          data = metrics,
          x = ~Value,
          y = ~reorder(Metric, Value),
          type = 'bar',
          orientation = 'h',
          marker = list(color = colors)
        ) %>%
          layout(
            xaxis = list(title = "Value"),
            yaxis = list(title = ""),
            showlegend = FALSE,
            margin = list(t = 10)
          )
      }

      return(p)
    }, error = function(e) {
      plotly_empty() %>%
        layout(title = list(text = paste("Error:", e$message),
                          font = list(color = "#999")))
    })
  })

  # Per-sample accuracy display
  output$sample_accuracy_display <- renderUI({
    if (is.null(data_store$validation_data)) {
      return(div(class = "no-data-message",
                icon("info-circle", class = "fa-3x"),
                br(), br(),
                tags$b("N/A - No validation data available"),
                br(),
                tags$small("Per-sample accuracy requires validation data")))
    }

    DTOutput("sample_accuracy_table")
  })

  output$sample_accuracy_table <- renderDT({
    # Try GLIMPSE accuracy first, then validation_data
    if (!is.null(data_store$tool_accuracy) && "Sample" %in% names(data_store$tool_accuracy)) {
      # Use GLIMPSE per-sample accuracy
      dt <- datatable(data_store$tool_accuracy,
                options = list(pageLength = 15, scrollX = TRUE),
                rownames = FALSE,
                colnames = c('Sample', 'N Variants', 'Best GT r²', 'Dosage r²',
                            'General Concordance %', 'Non-Ref Concordance %', 'Non-Ref Discordance %'))

      # Format numeric columns
      if ("N_Variants" %in% names(data_store$tool_accuracy)) {
        dt <- dt %>% formatCurrency("N_Variants", currency = "", digits = 0, mark = ",")
      }
      if ("Dosage_r2" %in% names(data_store$tool_accuracy)) {
        dt <- dt %>% formatRound("Dosage_r2", digits = 4)
      }
      if ("Best_GT_r2" %in% names(data_store$tool_accuracy)) {
        dt <- dt %>% formatRound("Best_GT_r2", digits = 4)
      }
      if ("General_Concordance_Pct" %in% names(data_store$tool_accuracy)) {
        dt <- dt %>% formatRound("General_Concordance_Pct", digits = 2)
      }
      if ("NonRef_Concordance_Pct" %in% names(data_store$tool_accuracy)) {
        dt <- dt %>% formatRound("NonRef_Concordance_Pct", digits = 2)
      }
      if ("NonRef_Discordance_Pct" %in% names(data_store$tool_accuracy)) {
        dt <- dt %>% formatRound("NonRef_Discordance_Pct", digits = 2)
      }

      return(dt)
    } else if (!is.null(data_store$validation_data)) {
      # Use old validation format
      tryCatch({
        sample_summary <- data_store$validation_data %>%
          group_by(ID) %>%
          summarise(
            Total_Variants = sum(n_genotypes, na.rm = TRUE),
            Mean_AF = mean(mean_AF, na.rm = TRUE),
            Total_Matches = sum(RR_hom_matches + RA_het_matches + AA_hom_matches, na.rm = TRUE),
            Total_Mismatches = sum(RR_hom_mismatches + RA_het_mismatches + AA_hom_mismatches, na.rm = TRUE),
            .groups = 'drop'
          ) %>%
          mutate(
            Accuracy_Pct = round((Total_Matches / (Total_Matches + Total_Mismatches)) * 100, 2)
          )

        datatable(sample_summary,
                  options = list(pageLength = 10, scrollX = TRUE),
                  rownames = FALSE) %>%
          formatRound(columns = c("Mean_AF"), digits = 4)
      }, error = function(e) {
        return(data.frame(Error = "Could not parse validation data"))
      })
    } else {
      return(NULL)
    }
  })

  # Chromosome accuracy and INFO scores removed - not feasible with current pipeline
}

# Run the application
shinyApp(ui = ui, server = server)
