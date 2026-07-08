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
      # Order: General Concordance, Non-Ref Concordance, Dosage r², Best GT r²
      result <- data.frame(
        Sample = data$Sample,  # Keep full sample name (e.g., "DUR_10004940.beagle5")
        N_Variants = sapply(parsed_data, function(x) x$n_variants),
        General_Concordance_Pct = sapply(parsed_data, function(x) x$general_concordance),
        NonRef_Concordance_Pct = sapply(parsed_data, function(x) x$non_ref_concordance),
        Dosage_r2 = sapply(parsed_data, function(x) x$imputed_ds_rsquared),
        Best_GT_r2 = sapply(parsed_data, function(x) x$best_gt_rsquared),
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
      tool <- if (grepl("Tbeagle5$|_beagle5$", sample_name)) {
        "beagle5"
      } else if (grepl("Tglimpse2$|_glimpse2$", sample_name)) {
        "glimpse2"
      } else if (grepl("Tglimpse1$|_glimpse1$", sample_name)) {
        "glimpse1"
      } else if (grepl("Tminimac4$|_minimac4$", sample_name)) {
        "minimac4"
      } else if (grepl("Tquilt$|_quilt$", sample_name)) {
        "quilt"
      } else if (grepl("Tstitch$|_stitch$", sample_name)) {
        "stitch"
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
              Sample = sample_name,
              stringsAsFactors = FALSE
            )
          }
        }
      }
    }

    if (length(result_list) == 0) return(NULL)

    # Combine all results - keep per-sample data
    all_data <- do.call(rbind, result_list)
    all_data$MAF <- as.numeric(all_data$MAF)
    all_data$Dosage_r2 <- as.numeric(all_data$Dosage_r2)

    # Return the per-sample data (aggregation done in plot if needed)
    return(all_data)
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

# Load per-chromosome stats from pipeline
load_per_chr_stats <- function(base_dir) {
  stats_file <- find_file(base_dir, "per_chromosome_stats.csv",
                          c("qc_stats/", "*/qc_stats/", "."))
  if (is.null(stats_file)) return(NULL)

  tryCatch({
    data <- read.csv(stats_file, stringsAsFactors = FALSE)
    # Expected columns: chr, tool, snps_before_imputation, snps_after_imputation,
    #                   general_concordance, dosage_r2, best_gt_r2, non_ref_concordance
    return(data)
  }, error = function(e) {
    return(NULL)
  })
}

# Helper function to clean sample names
# Converts "ISA_0555957_PREF_PANEL_Tbeagle5" to "ISA_0555957.beagle5"
clean_sample_name <- function(sample_name) {
  # Extract tool suffix
  tool <- ""
  if (grepl("Tbeagle5$|_beagle5$", sample_name)) tool <- "beagle5"
  else if (grepl("Tglimpse2$|_glimpse2$", sample_name)) tool <- "glimpse2"
  else if (grepl("Tglimpse1$|_glimpse1$", sample_name)) tool <- "glimpse1"
  else if (grepl("Tminimac4$|_minimac4$", sample_name)) tool <- "minimac4"
  else if (grepl("Tquilt$|_quilt$", sample_name)) tool <- "quilt"
  else if (grepl("Tstitch$|_stitch$", sample_name)) tool <- "stitch"

  # Remove tool suffix and intermediate parts like _PREF_PANEL_
  base_name <- gsub("(_PREF_PANEL)?_?T?(beagle5|glimpse2|glimpse1|minimac4|quilt|stitch)$", "", sample_name)

  # Return cleaned name
  if (tool != "") {
    paste0(base_name, ".", tool)
  } else {
    sample_name
  }
}

# Helper function to sort chromosomes naturally (1, 2, ..., 22, X, Y, MT)
# Handles both "chr1" and "1" formats
sort_chromosomes <- function(chroms) {
  # Remove chr prefix for sorting, keep track of original format
  has_chr_prefix <- any(grepl("^chr", chroms, ignore.case = TRUE))
  chroms_clean <- gsub("^chr", "", chroms, ignore.case = TRUE)

  # Extract numeric and non-numeric chromosomes
  numeric_chroms <- chroms_clean[grepl("^[0-9]+$", chroms_clean)]
  non_numeric <- chroms_clean[!grepl("^[0-9]+$", chroms_clean)]

  # Sort numeric chromosomes by value
  sorted_numeric <- as.character(sort(as.numeric(unique(numeric_chroms))))

  # Sort non-numeric alphabetically but put common ones first
  priority_order <- c("X", "Y", "MT", "M")
  priority_chroms <- unique(non_numeric[non_numeric %in% priority_order])
  other_chroms <- unique(non_numeric[!non_numeric %in% priority_order])

  # Order priority chroms by their position in priority_order
  priority_chroms <- priority_chroms[order(match(priority_chroms, priority_order))]

  sorted_clean <- c(sorted_numeric, priority_chroms, sort(other_chroms))

  # Add chr prefix back if original had it
  if (has_chr_prefix) {
    sorted_clean <- paste0("chr", sorted_clean)
  }

  sorted_clean
}

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
                 onclick = "switchTab('quality')"),
      tags$button("Per-Chromosome Stats", class = "navbar-tab",
                 `data-tab` = "chromosome",
                 onclick = "switchTab('chromosome')")
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
          style = "flex: 0 0 50%; min-width: 0;",
          tags$div(
            class = "section",
            id = "maf-plot",
            tags$h3("Imputation Quality by MAF"),
            tags$div(
              class = "content-box",
              tags$div(
                style = "margin-bottom: 10px;",
                checkboxInput("show_by_sample", "Show by sample", value = FALSE)
              ),
              plotlyOutput("maf_spectrum_plot", height = "450px"),
              uiOutput("maf_legend_box")
            )
          )
        ),

        # Right column - Accuracy metrics stacked
        tags$div(
          style = "flex: 0 0 calc(50% - 10px); min-width: 0; display: flex; flex-direction: column; gap: 10px;",
          tags$div(
            class = "section",
            id = "overall-accuracy",
            tags$h3("Overall Accuracy Metrics"),
            tags$div(
              class = "content-box",
              uiOutput("tool_accuracy_display")
            )
          ),
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
    ),

    # Per-Chromosome Stats tab content
    tags$div(
      id = "chromosome-content",
      class = "tab-content",
      style = "display: none;",

      # Two columns: Summary table (left) + Charts stacked (right)
      tags$div(
        style = "display: flex; gap: 10px; width: 100%; overflow: hidden;",

        # Left column - Summary table (flexible, shrinks to fit)
        tags$div(
          style = "flex: 1; min-width: 0; max-width: 55%;",
          tags$div(
            class = "section",
            id = "chr-summary",
            tags$h3("Per-Chromosome Summary"),
            tags$div(
              class = "content-box",
              style = "overflow-x: auto;",
              DTOutput("chr_stats_table")
            )
          )
        ),

        # Right column - Charts stacked
        tags$div(
          style = "flex: 1; min-width: 0; display: flex; flex-direction: column; gap: 10px;",
          tags$div(
            class = "section",
            id = "chr-snp-plot",
            tags$h3("SNP Counts by Chromosome"),
            tags$div(
              class = "content-box",
              plotlyOutput("chr_snp_plot", height = "300px")
            )
          ),
          tags$div(
            class = "section",
            id = "chr-accuracy",
            tags$h3("Accuracy by Chromosome"),
            tags$div(
              class = "content-box",
              plotlyOutput("chr_accuracy_plot", height = "300px")
            )
          )
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
    info_scores = NULL,
    per_chr_stats = NULL
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

      # Load per-chromosome stats
      incProgress(0.05, detail = "Loading per-chromosome stats")
      data_store$per_chr_stats <- load_per_chr_stats(base_dir)
    })

    # Show diagnostic message
    msg_parts <- c()
    if (!is.null(data_store$multiqc_data)) msg_parts <- c(msg_parts, "MultiQC")
    if (!is.null(data_store$params)) msg_parts <- c(msg_parts, "Params")
    if (!is.null(data_store$summary_stats)) msg_parts <- c(msg_parts, "Summary Stats")
    if (!is.null(data_store$tool_accuracy)) msg_parts <- c(msg_parts, paste0("Accuracy (", nrow(data_store$tool_accuracy), " rows)"))
    if (!is.null(data_store$maf_spectrum)) msg_parts <- c(msg_parts, paste0("MAF (", nrow(data_store$maf_spectrum), " rows)"))
    if (!is.null(data_store$per_chr_stats)) msg_parts <- c(msg_parts, paste0("Per-Chr (", nrow(data_store$per_chr_stats), " rows)"))

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
    } else if (active_tab == "quality") {
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
    } else {
      # Per-Chromosome Stats tab
      tagList(
        tags$a(href = "#chr-summary", class = "sidebar-nav-link",
              onclick = "scrollToSection('chr-summary'); return false;",
              "Per-Chromosome Summary"),
        tags$a(href = "#chr-snp-plot", class = "sidebar-nav-link",
              onclick = "scrollToSection('chr-snp-plot'); return false;",
              "SNP Count Chart"),
        tags$a(href = "#chr-accuracy", class = "sidebar-nav-link",
              onclick = "scrollToSection('chr-accuracy'); return false;",
              "Accuracy by Chromosome")
      )
    }
  })

  # Overview value boxes
  output$overview_value_boxes <- renderUI({

    panel_text <- if (!is.null(data_store$summary_stats)) {
      panel_samples <- data_store$summary_stats$Value[data_store$summary_stats$Metric == "Panel_Samples"]
      panel_snps <- data_store$summary_stats$Value[data_store$summary_stats$Metric == "Panel_SNPs"]
      if (length(panel_samples) > 0 && !is.na(panel_samples) && panel_samples != "NA" &&
          length(panel_snps) > 0 && !is.na(panel_snps) && panel_snps != "NA") {
        paste0(panel_samples, " samples\n", format_k(panel_snps), " SNPs")
      } else { "N/A" }
    } else { "N/A" }

    target_text <- if (!is.null(data_store$summary_stats)) {
      target_samples <- data_store$summary_stats$Value[data_store$summary_stats$Metric == "Target_Samples"]
      target_snps <- data_store$summary_stats$Value[data_store$summary_stats$Metric == "Target_SNPs"]
      if (length(target_samples) > 0 && !is.na(target_samples) && target_samples != "NA" &&
          length(target_snps) > 0 && !is.na(target_snps) && target_snps != "NA") {
        paste0(target_samples, " samples\n", format_k(target_snps), " SNPs")
      } else { "N/A" }
    } else { "N/A" }

    imputed_text <- if (!is.null(data_store$summary_stats)) {
      imputed_samples <- data_store$summary_stats$Value[data_store$summary_stats$Metric == "Imputed_Samples"]
      imputed_snps <- data_store$summary_stats$Value[data_store$summary_stats$Metric == "Imputed_SNPs"]
      if (length(imputed_samples) > 0 && !is.na(imputed_samples) && imputed_samples != "NA" &&
          length(imputed_snps) > 0 && !is.na(imputed_snps) && imputed_snps != "NA") {
        paste0(imputed_samples, " samples\n", format_k(imputed_snps), " SNPs")
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
        class = "value-box light-blue",
        tags$div(class = "value-box-value", panel_text),
        tags$div(class = "value-box-subtitle", "Reference Panel")
      ),      
      tags$div(
        class = "value-box navy",
        tags$div(class = "value-box-value", target_text),
        tags$div(class = "value-box-subtitle", "Target Data")
      ),
      tags$div(
        class = "value-box olive",
        tags$div(class = "value-box-value", imputed_text),
        tags$div(class = "value-box-subtitle", "Imputed data")
      ),
      tags$div(
        class = "value-box navy",
        tags$div(class = "value-box-value", truth_text),
        tags$div(class = "value-box-subtitle", "Validation Data")
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

    # Get panel samples from various sources
    if (!is.null(data_store$base_dir)) {
      panel_samples <- c()

      # Method 1: Try to get from REF_PANEL bcftools stats (when prep_panel was run)
      tryCatch({
        panel_stats_file <- find_file(data_store$base_dir, "REF_PANEL.panel.bcftools_stats.txt",
                                       c("prep_panel/stats/", "*/prep_panel/stats/", "."))
        if (!is.null(panel_stats_file) && file.exists(panel_stats_file)) {
          stats_lines <- readLines(panel_stats_file)
          # Look for sample lines - format: "PSC	0	sample_name	..."
          psc_lines <- grep("^PSC\t0\t", stats_lines, value = TRUE)
          if (length(psc_lines) > 0) {
            # Extract sample names from PSC lines
            panel_samples <- sapply(strsplit(psc_lines, "\t"), function(x) x[3])
          }
        }
      }, error = function(e) {})

      # Method 2: Try the original input panel CSV from params (when panel was pre-prepared)
      # Panel CSV has columns: panel, chromosome, file, index - need to query VCF for sample names
      if (length(panel_samples) == 0 && !is.null(data_store$params)) {
        tryCatch({
          panel_param <- data_store$params$panel
          if (!is.null(panel_param) && file.exists(panel_param)) {
            panel_csv_data <- read.csv(panel_param, stringsAsFactors = FALSE, header = TRUE)
            if (nrow(panel_csv_data) > 0) {
              # Get the VCF file path (column might be 'file' or 'vcf')
              vcf_col <- if ("file" %in% names(panel_csv_data)) "file" else if ("vcf" %in% names(panel_csv_data)) "vcf" else NULL
              if (!is.null(vcf_col)) {
                vcf_path <- panel_csv_data[[vcf_col]][1]
                if (file.exists(vcf_path)) {
                  # Try bcftools if available
                  bcftools_result <- tryCatch({
                    sample_cmd <- paste0("bcftools query -l ", shQuote(vcf_path))
                    system(sample_cmd, intern = TRUE, ignore.stderr = TRUE)
                  }, error = function(e) NULL, warning = function(w) NULL)

                  if (!is.null(bcftools_result) && length(bcftools_result) > 0) {
                    panel_samples <- bcftools_result
                  }
                }
              }
            }
          }
        }, error = function(e) {})
      }

      # Method 3: Try output panel.csv + bcftools (if available)
      if (length(panel_samples) == 0) {
        tryCatch({
          panel_csv <- find_file(data_store$base_dir, "panel.csv",
                                 c("prep_panel/csv/", "*/prep_panel/csv/", "."))

          if (!is.null(panel_csv) && file.exists(panel_csv)) {
            csv_data <- read.csv(panel_csv, stringsAsFactors = FALSE, header = TRUE)
            if (nrow(csv_data) > 0 && "vcf" %in% names(csv_data)) {
              vcf_path <- csv_data$vcf[1]

              if (file.exists(vcf_path)) {
                # Try bcftools if available
                bcftools_result <- tryCatch({
                  sample_cmd <- paste0("bcftools query -l ", shQuote(vcf_path))
                  system(sample_cmd, intern = TRUE, ignore.stderr = TRUE)
                }, error = function(e) NULL, warning = function(w) NULL)

                if (!is.null(bcftools_result) && length(bcftools_result) > 0) {
                  panel_samples <- bcftools_result
                }
              }
            }
          }
        }, error = function(e) {})
      }

      # Add panel samples to the list
      if (length(panel_samples) > 0) {
        for (panel_sample in panel_samples) {
          if (panel_sample %in% sample_list$Sample) {
            idx <- which(sample_list$Sample == panel_sample)
            sample_list$In_Panel[idx] <- "Yes"
          } else {
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
      show_by_sample <- input$show_by_sample

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
      maf_data <- maf_data[complete.cases(maf_data[, c("MAF", "Dosage_r2", "Tool")]), ]

      if (nrow(maf_data) == 0) {
        return(plotly_empty() %>%
                 layout(title = list(text = "No valid MAF data after filtering",
                                   font = list(color = "#999", size = 14))))
      }

      # Tool colors
      tool_colors <- c(
        "beagle5" = "#2091ae",    # Light blue
        "glimpse2" = "#9cbe2b",   # Green
        "glimpse1" = "#16315b",   # Navy
        "minimac4" = "#e76f51",   # Orange
        "quilt" = "#f4a261",      # Light orange
        "stitch" = "#8338ec",     # Purple
        "unknown" = "#666666"
      )

      p <- plot_ly()

      if (show_by_sample && "Sample" %in% names(maf_data)) {
        # Show individual sample lines
        samples <- unique(maf_data$Sample)
        for (samp in samples) {
          sample_data <- maf_data[maf_data$Sample == samp, ]
          sample_data <- sample_data[order(sample_data$MAF), ]
          tool <- sample_data$Tool[1]
          tool_color <- if (tool %in% names(tool_colors)) tool_colors[[tool]] else "#666666"
          clean_name <- clean_sample_name(samp)

          p <- p %>%
            add_trace(
              data = sample_data,
              x = ~MAF,
              y = ~Dosage_r2,
              name = clean_name,
              type = 'scatter',
              mode = 'lines',
              line = list(color = tool_color, width = 1),
              opacity = 0.5,
              hovertemplate = paste0(
                '<b>', clean_name, '</b><br>',
                'MAF: %{x:.3f}<br>',
                'Dosage r²: %{y:.2f}<br>',
                '<extra></extra>'
              )
            )
        }
      } else {
        # Aggregate by tool and MAF bin
        maf_data$MAF_rounded <- round(maf_data$MAF, 3)
        tools <- unique(maf_data$Tool)

        for (tool in tools) {
          tool_subset <- maf_data[maf_data$Tool == tool, ]
          unique_mafs <- unique(tool_subset$MAF_rounded)

          agg_list <- list()
          for (maf in unique_mafs) {
            maf_subset <- tool_subset[tool_subset$MAF_rounded == maf, ]
            agg_list[[length(agg_list) + 1]] <- data.frame(
              MAF = maf,
              Dosage_r2 = mean(maf_subset$Dosage_r2, na.rm = TRUE),
              stringsAsFactors = FALSE
            )
          }
          agg_data <- do.call(rbind, agg_list)
          agg_data <- agg_data[order(agg_data$MAF), ]

          tool_color <- if (tool %in% names(tool_colors)) tool_colors[[tool]] else "#666666"

          p <- p %>%
            add_trace(
              data = agg_data,
              x = ~MAF,
              y = ~Dosage_r2,
              name = tool,
              type = 'scatter',
              mode = 'lines+markers',
              line = list(color = tool_color, width = 2),
              marker = list(color = tool_color, size = 8, line = list(color = 'white', width = 1.5)),
              hovertemplate = paste0(
                '<b>', tool, '</b><br>',
                'MAF: %{x:.3f}<br>',
                'Dosage r²: %{y:.2f}<br>',
                '<extra></extra>'
              )
            )
        }
      }

      p %>%
        layout(
          xaxis = list(
            title = "Minor Allele Frequency (MAF)",
            range = c(0, 0.5),
            tickmode = "linear",
            tick0 = 0,
            dtick = 0.05,
            gridcolor = '#e8e8e8',
            showgrid = TRUE
          ),
          yaxis = list(
            title = "Dosage r²",
            range = c(0, 1),
            tickmode = "linear",
            tick0 = 0,
            dtick = 0.1,
            gridcolor = '#e8e8e8',
            showgrid = TRUE
          ),
          hovermode = "x unified",
          showlegend = FALSE,
          plot_bgcolor = 'white',
          paper_bgcolor = 'white',
          margin = list(t = 10, b = 50, l = 60, r = 10)
        )
    }, error = function(e) {
      plotly_empty() %>%
        layout(title = list(text = paste("Error creating plot:", e$message),
                          font = list(color = "#999", size = 14)))
    })
  })

  # MAF plot legend box (separate from plot)
  output$maf_legend_box <- renderUI({
    if (is.null(data_store$maf_spectrum)) return(NULL)

    maf_data <- as.data.frame(data_store$maf_spectrum)
    show_by_sample <- input$show_by_sample

    # Tool colors
    tool_colors <- c(
      "beagle5" = "#2091ae",
      "glimpse2" = "#9cbe2b",
      "glimpse1" = "#16315b",
      "minimac4" = "#e76f51",
      "quilt" = "#f4a261",
      "stitch" = "#8338ec",
      "unknown" = "#666666"
    )

    if (show_by_sample && "Sample" %in% names(maf_data)) {
      # Show sample legend with tool colors
      samples <- unique(maf_data$Sample)
      legend_items <- lapply(samples, function(samp) {
        sample_data <- maf_data[maf_data$Sample == samp, ]
        tool <- sample_data$Tool[1]
        color <- if (tool %in% names(tool_colors)) tool_colors[[tool]] else "#666666"
        clean_name <- clean_sample_name(samp)
        tags$span(
          style = paste0("display: inline-block; margin-right: 12px; margin-bottom: 4px; font-size: 11px;"),
          tags$span(style = paste0("display: inline-block; width: 12px; height: 3px; background-color: ", color, "; margin-right: 4px; vertical-align: middle;")),
          clean_name
        )
      })

      tags$div(
        style = "margin-top: 10px; padding: 10px; background-color: #f8f8f8; border: 1px solid #ddd; border-radius: 4px;",
        tags$strong(style = "font-size: 11px; display: block; margin-bottom: 6px;", "Samples:"),
        tags$div(style = "line-height: 1.8;", legend_items)
      )
    } else {
      # Show tool legend
      tools <- unique(maf_data$Tool)
      legend_items <- lapply(tools, function(tool) {
        color <- if (tool %in% names(tool_colors)) tool_colors[[tool]] else "#666666"
        tags$span(
          style = paste0("display: inline-block; margin-right: 15px; font-size: 12px;"),
          tags$span(style = paste0("display: inline-block; width: 20px; height: 3px; background-color: ", color, "; margin-right: 6px; vertical-align: middle;")),
          tool
        )
      })

      tags$div(
        style = "margin-top: 10px; padding: 10px; background-color: #f8f8f8; border: 1px solid #ddd; border-radius: 4px;",
        tags$strong(style = "font-size: 12px; margin-right: 10px;", "Tool:"),
        legend_items
      )
    }
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

        # Check if we have the expected columns
        required_cols <- c("Sample", "Dosage_r2", "Best_GT_r2", "General_Concordance_Pct", "NonRef_Concordance_Pct")
        missing_cols <- setdiff(required_cols, names(tool_data))

        if (length(missing_cols) > 0) {
          return(data.frame(
            Error = paste("Missing columns:", paste(missing_cols, collapse = ", "),
                         ". Available:", paste(names(tool_data), collapse = ", "))
          ))
        }

        # Extract tool name from sample
        tool_data$Tool <- sapply(tool_data$Sample, function(s) {
          if (grepl("Tbeagle5$|_beagle5$|\\.beagle5$", s)) return("beagle5")
          if (grepl("Tglimpse2$|_glimpse2$|\\.glimpse2$", s)) return("glimpse2")
          if (grepl("Tglimpse1$|_glimpse1$|\\.glimpse1$", s)) return("glimpse1")
          if (grepl("Tminimac4$|_minimac4$|\\.minimac4$", s)) return("minimac4")
          if (grepl("Tquilt$|_quilt$|\\.quilt$", s)) return("quilt")
          if (grepl("Tstitch$|_stitch$|\\.stitch$", s)) return("stitch")
          return("unknown")
        })

        # Convert metrics to numeric
        tool_data$Dosage_r2 <- as.numeric(tool_data$Dosage_r2)
        tool_data$Best_GT_r2 <- as.numeric(tool_data$Best_GT_r2)
        tool_data$General_Concordance_Pct <- as.numeric(tool_data$General_Concordance_Pct)
        tool_data$NonRef_Concordance_Pct <- as.numeric(tool_data$NonRef_Concordance_Pct)

        # Aggregate by tool using base R
        # Order: General Concordance, Non-Ref Concordance, Dosage r², Best GT r²
        tools <- unique(tool_data$Tool)
        summary_list <- lapply(tools, function(t) {
          tool_subset <- tool_data[tool_data$Tool == t, ]
          data.frame(
            Tool = t,
            N_Samples = nrow(tool_subset),
            Mean_General_Concordance_Pct = round(mean(tool_subset$General_Concordance_Pct, na.rm = TRUE), 1),
            Mean_NonRef_Concordance_Pct = round(mean(tool_subset$NonRef_Concordance_Pct, na.rm = TRUE), 1),
            Mean_Dosage_r2 = round(mean(tool_subset$Dosage_r2, na.rm = TRUE), 2),
            Mean_Best_GT_r2 = round(mean(tool_subset$Best_GT_r2, na.rm = TRUE), 2),
            stringsAsFactors = FALSE
          )
        })
        summary_data <- do.call(rbind, summary_list)

        datatable(as.data.frame(summary_data),
                  options = list(pageLength = 10, dom = 't', scrollX = TRUE),
                  rownames = FALSE,
                  colnames = c("Tool", "N Samples", "General Concordance %", "Non-Ref Concordance %",
                              "Dosage r²", "Best GT r²"))
      } else {
        # Original format
        datatable(as.data.frame(data_store$tool_accuracy),
                  options = list(pageLength = 10, dom = 't', scrollX = TRUE),
                  rownames = FALSE) %>%
          formatRound(columns = grep("r2|concordance|pct", names(data_store$tool_accuracy), ignore.case = TRUE),
                      digits = 2)
      }
    }, error = function(e) {
      return(data.frame(Error = paste("Could not parse accuracy data:", e$message)))
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
      # Order: General Concordance, Non-Ref Concordance, Dosage r², Best GT r²
      dt <- datatable(data_store$tool_accuracy,
                options = list(pageLength = 15, scrollX = TRUE),
                rownames = FALSE,
                colnames = c('Sample', 'N Variants', 'General Concordance %', 'Non-Ref Concordance %',
                            'Dosage r²', 'Best GT r²'))

      # Format numeric columns: 1 decimal for percentages, 2 decimals for r²
      if ("N_Variants" %in% names(data_store$tool_accuracy)) {
        dt <- dt %>% formatCurrency("N_Variants", currency = "", digits = 0, mark = ",")
      }
      if ("General_Concordance_Pct" %in% names(data_store$tool_accuracy)) {
        dt <- dt %>% formatRound("General_Concordance_Pct", digits = 1)
      }
      if ("NonRef_Concordance_Pct" %in% names(data_store$tool_accuracy)) {
        dt <- dt %>% formatRound("NonRef_Concordance_Pct", digits = 1)
      }
      if ("Dosage_r2" %in% names(data_store$tool_accuracy)) {
        dt <- dt %>% formatRound("Dosage_r2", digits = 2)
      }
      if ("Best_GT_r2" %in% names(data_store$tool_accuracy)) {
        dt <- dt %>% formatRound("Best_GT_r2", digits = 2)
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
          formatRound(columns = c("Mean_AF"), digits = 2) %>%
          formatRound(columns = c("Accuracy_Pct"), digits = 1)
      }, error = function(e) {
        return(data.frame(Error = "Could not parse validation data"))
      })
    } else {
      return(NULL)
    }
  })

  # Per-chromosome stats table
  output$chr_stats_table <- renderDT({
    if (is.null(data_store$per_chr_stats)) {
      return(data.frame(Info = "No per-chromosome data available. Run the pipeline with --steps impute or validate to generate per_chromosome_stats.csv"))
    }

    df <- data_store$per_chr_stats

    # Sort chromosomes naturally
    df$chr <- factor(df$chr, levels = sort_chromosomes(unique(df$chr)))
    df <- df[order(df$chr, df$tool), ]

    # Convert concordance from decimal (0-1) to percentage (0-100)
    df$general_concordance <- as.numeric(df$general_concordance) * 100
    df$non_ref_concordance <- as.numeric(df$non_ref_concordance) * 100

    # Reorder columns: chr, tool, snps_before, snps_after, general_conc, non_ref_conc, dosage_r2, best_gt_r2
    df <- df[, c("chr", "tool", "snps_before_imputation", "snps_after_imputation",
                 "general_concordance", "non_ref_concordance", "dosage_r2", "best_gt_r2")]

    # Format column names for display (with line breaks for long headers)
    colnames(df) <- c("Chr", "Tool", "SNPs<br>Before", "SNPs<br>After",
                      "General<br>Conc %", "Non-Ref<br>Conc %", "Dosage<br>r²", "Best GT<br>r²")

    datatable(df,
              options = list(
                pageLength = 25,
                scrollX = TRUE,
                dom = 'frtip',
                autoWidth = FALSE
              ),
              rownames = FALSE,
              escape = FALSE) %>%  # escape=FALSE allows HTML in column names
      formatRound(columns = c("General<br>Conc %", "Non-Ref<br>Conc %"), digits = 1) %>%
      formatRound(columns = c("Dosage<br>r²", "Best GT<br>r²"), digits = 2) %>%
      formatCurrency(columns = c("SNPs<br>Before", "SNPs<br>After"),
                    currency = "", digits = 0, mark = ",")
  })

  # Per-chromosome SNP counts plot
  output$chr_snp_plot <- renderPlotly({
    if (is.null(data_store$per_chr_stats)) {
      return(plotly_empty() %>%
               layout(title = list(text = "N/A - No per-chromosome data available",
                                 font = list(color = "#999", size = 14))))
    }

    df <- data_store$per_chr_stats

    # Convert SNP columns to numeric, handling "NA" strings
    df$snps_before_imputation <- as.numeric(df$snps_before_imputation)
    df$snps_after_imputation <- as.numeric(df$snps_after_imputation)

    # Sort chromosomes naturally
    chr_order <- sort_chromosomes(unique(df$chr))
    df$chr <- factor(df$chr, levels = chr_order)

    # Get unique tools
    tools <- unique(df$tool)

    # Create grouped bar chart
    p <- plot_ly()

    # Add bars for each tool
    tool_colors <- c(
      "beagle5" = "#2091ae",
      "glimpse2" = "#9cbe2b",
      "glimpse1" = "#16315b",
      "minimac4" = "#e76f51",
      "quilt" = "#f4a261",
      "stitch" = "#8338ec"
    )

    for (tool in tools) {
      tool_data <- df[df$tool == tool, ]
      tool_color <- if (tool %in% names(tool_colors)) tool_colors[[tool]] else "#666666"

      # Before imputation bars (lighter shade)
      p <- p %>%
        add_trace(
          data = tool_data,
          x = ~chr,
          y = ~snps_before_imputation,
          name = paste(tool, "- Before"),
          type = 'bar',
          marker = list(color = tool_color, opacity = 0.5),
          hovertemplate = paste0(
            '<b>', tool, ' - Before</b><br>',
            'Chr: %{x}<br>',
            'SNPs: %{y:,.0f}<br>',
            '<extra></extra>'
          )
        )

      # After imputation bars (solid)
      p <- p %>%
        add_trace(
          data = tool_data,
          x = ~chr,
          y = ~snps_after_imputation,
          name = paste(tool, "- After"),
          type = 'bar',
          marker = list(color = tool_color),
          hovertemplate = paste0(
            '<b>', tool, ' - After</b><br>',
            'Chr: %{x}<br>',
            'SNPs: %{y:,.0f}<br>',
            '<extra></extra>'
          )
        )
    }

    p %>%
      layout(
        xaxis = list(title = "", tickangle = 45),
        yaxis = list(title = "SNP Count"),
        barmode = 'group',
        legend = list(orientation = "h", y = -0.3, x = 0.5, xanchor = "center"),
        hovermode = "x unified",
        margin = list(b = 80)
      )
  })

  # Per-chromosome accuracy plot
  output$chr_accuracy_plot <- renderPlotly({
    if (is.null(data_store$per_chr_stats)) {
      return(plotly_empty() %>%
               layout(title = list(text = "N/A - No per-chromosome data available",
                                 font = list(color = "#999", size = 14))))
    }

    df <- data_store$per_chr_stats

    # Convert accuracy columns to numeric, handling "NA" strings
    df$dosage_r2 <- as.numeric(df$dosage_r2)

    # Check if all accuracy values are NA (no validation data)
    if (all(is.na(df$dosage_r2))) {
      return(plotly_empty() %>%
               layout(title = list(text = "N/A - No validation data available (run with --input_truth)",
                                 font = list(color = "#999", size = 14))))
    }

    # Sort chromosomes naturally
    chr_order <- sort_chromosomes(unique(df$chr))
    df$chr <- factor(df$chr, levels = chr_order)

    # Get unique tools
    tools <- unique(df$tool)

    tool_colors <- c(
      "beagle5" = "#2091ae",
      "glimpse2" = "#9cbe2b",
      "glimpse1" = "#16315b",
      "minimac4" = "#e76f51",
      "quilt" = "#f4a261",
      "stitch" = "#8338ec"
    )

    p <- plot_ly()

    for (tool in tools) {
      tool_data <- df[df$tool == tool, ]
      tool_data <- tool_data[order(tool_data$chr), ]
      tool_color <- if (tool %in% names(tool_colors)) tool_colors[[tool]] else "#666666"

      # Dosage r² line only
      p <- p %>%
        add_trace(
          data = tool_data,
          x = ~chr,
          y = ~dosage_r2,
          name = tool,
          type = 'scatter',
          mode = 'lines+markers',
          line = list(color = tool_color, width = 2),
          marker = list(color = tool_color, size = 6),
          hovertemplate = paste0(
            '<b>', tool, '</b><br>',
            'Chr: %{x}<br>',
            'Dosage r²: %{y:.2f}<br>',
            '<extra></extra>'
          )
        )
    }

    p %>%
      layout(
        xaxis = list(title = "", tickangle = 45),
        yaxis = list(title = "Dosage r²", range = c(
          max(0, min(df$dosage_r2, na.rm = TRUE) - 0.05),
          1
        )),
        legend = list(orientation = "h", y = -0.3, x = 0.5, xanchor = "center"),
        hovermode = "x unified",
        margin = list(b = 80)
      )
  })
}

# Run the application
shinyApp(ui = ui, server = server)
