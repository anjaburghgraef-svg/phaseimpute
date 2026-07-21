# Testing phaseimpute on Azure Batch

## Prerequisites

1. Azure credentials filled in `nextflow.config` (ask IT for these)
2. Access to the test datasets in this folder
3. Nextflow installed

## Step 1: Verify Config

First check that the Azure config loads correctly:

```bash
cd /path/to/phaseimpute
nextflow config -profile azure --species layer --sample_scale low
```

You should see `executor = 'azurebatch'` and `machineType` values in the output.

## Step 2: Run a Test

### Small test (recommended first)

```bash
nextflow run main.nf \
    -profile azure \
    --species layer \
    --sample_scale low \
    --input <path-to-test-input.csv> \
    --panel <path-to-test-panel.vcf.gz> \
    --outdir results_test \
    --tools beagle5
```

### What to check

- Jobs appear in Azure Batch portal
- No "out of memory" errors
- Pipeline completes without errors
- Output files are created in `results_test/`

## Test Datasets

<!-- TODO: Fill in test dataset locations and descriptions -->

| Dataset | Location | Species | Samples | Description |
|---------|----------|---------|---------|-------------|
| | | | | |

## Expected Runtime

| Test | Approximate Time |
|------|------------------|
| Small test (layer, low) | ~X hours |
| Medium test | ~X hours |

## Common Issues

### "Insufficient quota"
The 350 CPU limit was exceeded. Wait for other jobs to finish or reduce parallelism.

### "Pool not found"
The Azure pools need to be pre-created by IT. Contact them.

### Jobs stuck in queue
Check Azure Batch portal for error messages.

## Verifying Results

After the run completes:

1. Check `results_test/` folder exists with output files
2. Look at `results_test/pipeline_info/execution_report_*.html` for timing/resource info
3. Compare imputation accuracy if you have truth data

## Notes

- The pipeline does NOT create Azure pools (they must exist already)
- CPU limit: 350 total across all jobs
- If something fails, check `.nextflow.log` for details
