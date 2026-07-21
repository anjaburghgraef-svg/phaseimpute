# Azure Batch Configuration for phaseimpute

This folder contains Azure-specific configuration files for running the phaseimpute pipeline on Azure Batch.

## Quick Start

```bash
nextflow run main.nf -profile azure --species layer --sample_scale medium --outdir results
```

## How It Works

The Azure VM allocation works the same way as the HPC (slurm) resource allocation:

1. **`--species`** determines VMs for genome-size dependent processes (BEAGLE5, VCFFIXUP, etc.)
2. **`--sample_scale`** determines VMs for sample-count dependent processes (PLUGINSPLIT, PLINK_VCF, etc.)

### Config Loading Order

```
base.config              → Default VM mappings for all processes
    ↓
species/{species}.config → Overrides for genome-dependent processes
    ↓
sample_scale/{scale}.config → Overrides for sample-dependent processes
```

Later configs override earlier ones, so sample_scale settings take priority for processes like PLUGINSPLIT.

## Parameters

| Parameter | Options | Description |
|-----------|---------|-------------|
| `--species` | layer, broiler, turkey, pig, salmon, trout, shrimp | Species (determines genome-size resources) |
| `--sample_scale` | low, medium, high | Sample count category (default: medium) |

### Sample Scale Guide

| Scale | Sample Count | Example Use Case |
|-------|--------------|------------------|
| low | < 500 | Small test runs, pilot studies |
| medium | 500 - 2000 | Standard production runs |
| high | > 2000 | Large-scale imputation (e.g., 5000 pigs) |

## Available Azure VMs

Hendrix Azure has the following D-series VMs available:

| VM Size | vCPU | Memory | Used For |
|---------|------|--------|----------|
| Standard_D2ds_v4 | 2 | 8 GB | Light tasks (GAWK, INDEX, TABIX) |
| Standard_D4ds_v4 | 4 | 16 GB | Medium tasks (BCFTOOLS_NORM) |
| Standard_D8ds_v4 | 8 | 32 GB | Heavier tasks (BEAGLE5 for chicken) |
| Standard_D16ds_v4 | 16 | 64 GB | High memory (BEAGLE5 for pig) |
| Standard_D32ds_v4 | 32 | 128 GB | Very high memory (PLUGINSPLIT high) |
| Standard_D64ds_v4 | 64 | 256 GB | Maximum (fallback for extreme cases) |

**CPU Limit:** 350 CPUs total across all running jobs

## Before Running

### 1. Fill in Azure Credentials

Edit `nextflow.config` and replace the placeholder values in the azure profile:

```groovy
azure {
    batch {
        accountName = '<batch-account>'   // ← Replace with actual account name
        accountKey = '<batch-key>'        // ← Replace with actual key
    }
    storage {
        accountName = '<storage-account>' // ← Replace with actual account name
        accountKey = '<storage-key>'      // ← Replace with actual key
    }
}
```

Ask the IT/cloud team for these credentials if you don't have them.

### 2. Verify Config

Check that the config loads correctly:

```bash
nextflow config -profile azure --species layer --sample_scale medium
```

This prints the merged configuration — verify the `machineType` values look correct.

## File Structure

```
conf/azure/
├── README.md                    # This file
├── base.config                  # Default VM mappings
├── species/
│   ├── layer.config             # Layer chicken (~1 GB genome)
│   ├── broiler.config           # Broiler chicken (~1 GB genome)
│   ├── turkey.config            # Turkey (~1.1 GB genome)
│   ├── pig.config               # Pig (~2.5 GB genome, 65M variants)
│   ├── salmon.config            # Atlantic salmon (~3 GB genome)
│   ├── trout.config             # Rainbow trout (~2.4 GB genome)
│   └── shrimp.config            # Pacific white shrimp (~2 GB genome)
└── sample_scale/
    ├── low.config               # < 500 samples
    ├── medium.config            # 500-2000 samples (default)
    └── high.config              # > 2000 samples
```

## Example Commands

<!-- TODO: Add example commands with test datasets -->

## Troubleshooting

### Job fails with "out of memory"
Increase the VM size in the relevant config file. Check which process failed and find it in the species or sample_scale config.

### Jobs queue but don't start
Check the CPU quota (350 limit). Too many parallel jobs may exceed the quota.

### Config not loading correctly
Make sure `--species` and `--sample_scale` parameters are spelled correctly (lowercase).

## Cost Bottlenecks

These processes are the most expensive — monitor them:

| Process | Why Expensive |
|---------|---------------|
| BCFTOOLS_PLUGINSPLIT | 5+ days for large sample counts, uses D32/D64 |
| BEAGLE5_BEAGLE | High memory, long runtime per chromosome |
| PLINK_VCF | Scales badly with sample count |
| BCFTOOLS_STATS | Runs once per sample (thousands of tasks) |

## Contact

For questions about the pipeline: [Your contact info]
For Azure/infrastructure issues: [IT contact]
