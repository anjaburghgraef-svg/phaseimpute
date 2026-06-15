#!/usr/bin/env Rscript

# Install dependencies for phaseimpute Shiny QC Report

cat("Installing required R packages for phaseimpute Shiny QC Report...\n\n")

# Ensure user library exists (for Windows permission issues)
user_lib <- Sys.getenv("R_LIBS_USER")
if (!dir.exists(user_lib)) {
  cat(sprintf("Creating user library at: %s\n", user_lib))
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}

# Make sure we're using the user library
if (!user_lib %in% .libPaths()) {
  .libPaths(c(user_lib, .libPaths()))
}

cat(sprintf("Installing packages to: %s\n\n", .libPaths()[1]))

# List of required packages
required_packages <- c(
  "shiny",
  "shinydashboard",
  "DT",
  "ggplot2",
  "plotly",
  "jsonlite",
  "dplyr",
  "tidyr"
)

# Function to check and install packages
install_if_missing <- function(package) {
  if (!require(package, character.only = TRUE, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", package))
    tryCatch({
      install.packages(package,
                      repos = "https://cran.rstudio.com/",
                      dependencies = TRUE,
                      lib = .libPaths()[1])
      if (require(package, character.only = TRUE, quietly = TRUE)) {
        cat(sprintf("✓ %s installed successfully\n", package))
        return(TRUE)
      } else {
        cat(sprintf("✗ Failed to load %s after installation\n", package))
        return(FALSE)
      }
    }, error = function(e) {
      cat(sprintf("✗ Error installing %s: %s\n", package, e$message))
      return(FALSE)
    })
  } else {
    cat(sprintf("✓ %s already installed\n", package))
    return(TRUE)
  }
}

# Install all packages
cat("Checking and installing packages...\n")
cat("=====================================\n\n")

success <- sapply(required_packages, install_if_missing)

cat("\n=====================================\n")
if (all(success)) {
  cat("All dependencies installed successfully!\n")
  cat("\nTo launch the app, run:\n")
  cat("  Rscript launch_report.R\n")
  cat("\nOr in R:\n")
  cat("  shiny::runApp()\n")
} else {
  cat("Some packages failed to install. Please install them manually:\n")
  cat(paste(required_packages[!success], collapse = ", "), "\n")
}
