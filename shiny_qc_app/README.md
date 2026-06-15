# phaseimpute QC App

Interactive quality control dashboard for phaseimpute pipeline results.

---

## Quick Start

### 1. Launch the App

(Might be preferable to start the app in a new terminal)
Make sure you are in the pipeline directory

**On Linux/Mac/Cluster:**
```bash
cd shiny_qc_app
bash launch_app.sh
```

**On Windows:**
```powershell
cd shiny_qc_app
.\launch_app.bat
```

The first time the R dependencies will be downloaded, might take a while
The app will open in your browser at `http://localhost:3838`
  (you might get a notification asking if the app can be opened)

### 2. Load Your Data

In the app interface:
1. Enter your pipeline output directory path (e.g., `/path/to/output` or `results/run1`)
2. Click **"Load Data"**
3. View your results!

---

## What You'll See

### Overview Tab
- **Sample Counts**: Imputed samples, validation samples, ref panel samples
- **SNP Counts**: Before and after imputation
- **Tool Used**: Imputation tool (BEAGLE5, GLIMPSE, etc.)
- **Sample Overview**: Detailed breakdown with full sample list
- **Pipeline Parameters**: All parameters used in the run

### Quality & Accuracy Tab
- **MAF Spectrum**: Imputation quality across minor allele frequencies
- **Overall Accuracy**: Mean metrics across all samples
- **Quality Summary**: Visual comparison of key metrics
- **Per-Sample Accuracy**: Detailed accuracy for each validated sample

---

## Requirements

**R (version 4.0+)** with the following packages (installed automatically on first run):
- shiny
- shinydashboard
- DT
- ggplot2
- plotly
- jsonlite
- dplyr
- tidyr

---

## Troubleshooting

### "No data available" message
✓ Check that your output directory path is correct  
✓ Ensure the pipeline has completed successfully  
✓ Verify these folders exist in your output:
  - `multiqc/multiqc_data/`
  - `pipeline_info/`

### "See params" in value boxes
This is normal if:
- `qc_stats/summary_stats.txt` doesn't exist (panelprep step not run)
- Input/panel CSV files are not accessible from where you're running the app

### R packages installation fails
Try installing manually:
```r
Rscript install_dependencies.R
```

### Port 3838 already in use
Stop the existing app or change the port in the launch script:
```bash
# Change 3838 to another port like 3839
Rscript -e "shiny::runApp('app.R', port=3839, launch.browser=TRUE)"
```

### App won't start - missing R
**On cluster:**
```bash
module avail R          # Check available R versions
module load R/4.3.1     # Load R (adjust version)
./launch_app.sh
```

**On Windows:**  
Install R from: https://cran.r-project.org/

---

## Example Usage

```bash
# After running your pipeline
nextflow run nf-core/phaseimpute \
  --input samples.csv \
  --panel refpanel.csv \
  --outdir results/my_run \
  ... other params ...

# Launch QC app
cd shiny_qc_app
./launch_app.sh

# In the app, enter: results/my_run
# Click "Load Data"
```

---

## Tips

💡 **Multiple Runs**: You can point the app to different output directories without restarting - just change the path and click "Load Data" again

## Support

For issues or questions:
- Check pipeline documentation: https://nf-co.re/phaseimpute
- Verify your pipeline completed without errors
- Check `.nextflow.log` for pipeline issues
- Ensure all required output files exist

---

**App Version:** 1.0  
**Compatible with:** nf-core/phaseimpute v1.0+  
**Last Updated:** June 2026
