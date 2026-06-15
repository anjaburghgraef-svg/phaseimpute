#!/bin/bash

echo "=========================================="
echo "phaseimpute QC Report (using conda)"
echo "=========================================="
echo ""

# Load conda if available
if command -v module &> /dev/null; then
    module load Anaconda3 2>/dev/null || module load conda 2>/dev/null || module load miniconda 2>/dev/null
fi

# Check if conda environment exists, create if not
if [ ! -d "$HOME/conda_envs/r_shiny" ]; then
    echo "Creating conda environment (one-time setup)..."
    conda create -y -p $HOME/conda_envs/r_shiny -c conda-forge \
        r-base=4.3 r-shiny r-shinydashboard r-dt r-ggplot2 r-plotly \
        r-jsonlite r-dplyr r-tidyr r-rcolorbrewer
fi

echo "Activating conda environment..."
source activate $HOME/conda_envs/r_shiny

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

# Run the app
Rscript -e "shiny::runApp('app.R', host = '0.0.0.0', port = 3838, launch.browser = FALSE)"
