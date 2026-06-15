#!/bin/bash

echo "=========================================="
echo "phaseimpute QC Report"
echo "=========================================="
echo ""

# Try to load R module (on HPC systems)
if command -v module &> /dev/null; then
    echo "Loading R module..."
    # Try R 4.4.2 first (compatible with ICU 75), then fallback to older versions
    module load 2024 2>/dev/null && module load R/4.4.2-gfbf-2024a 2>/dev/null || \
    module load R/4.3.1 2>/dev/null || \
    module load R/4.2 2>/dev/null || \
    module load R 2>/dev/null
fi

# Check if R is available
if ! command -v R &> /dev/null; then
    echo "ERROR: R is not available"
    echo "Please load R module manually: module load R"
    exit 1
fi

# Get R library path
R_VERSION=$(R --version | head -n1 | sed 's/R version //' | cut -d' ' -f1 | cut -d'.' -f1-2)
R_LIB_DIR="$HOME/R/x86_64-pc-linux-gnu-library/$R_VERSION"

# Create R library directory if it doesn't exist
mkdir -p "$R_LIB_DIR"

# Check and install required packages
echo "Installing R dependencies (first time only)..."
Rscript -e "
lib_path <- '$R_LIB_DIR'
.libPaths(c(lib_path, .libPaths()))

required_packages <- c('shiny', 'shinydashboard', 'DT', 'ggplot2', 'plotly',
                       'jsonlite', 'dplyr', 'tidyr', 'RColorBrewer')

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat('Installing', pkg, '...\n')
    install.packages(pkg, lib = lib_path, repos = 'https://cloud.r-project.org/', quiet = TRUE)
  }
}
" 2>&1 | grep -v "^>"

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to install R packages"
    echo "Try running manually: R"
    echo "Then in R console: install.packages(c('shiny', 'shinydashboard', 'DT', 'ggplot2', 'plotly', 'jsonlite', 'dplyr', 'tidyr', 'RColorBrewer'))"
    exit 1
fi

echo ""
echo "Starting Shiny app..."
echo "Access at: http://localhost:3838"
echo ""
echo "In the app:"
echo "  1. Enter your output directory path"
echo "  2. Click 'Load Data'"
echo ""
echo "Press Ctrl+C to stop"
echo "=========================================="
echo ""

# Set R library path and run app
export R_LIBS_USER="$R_LIB_DIR"

# Run Rscript with the app
Rscript -e "
.libPaths(c('$R_LIB_DIR', .libPaths()))
shiny::runApp('app.R', host = '0.0.0.0', port = 3838, launch.browser = FALSE)
"
