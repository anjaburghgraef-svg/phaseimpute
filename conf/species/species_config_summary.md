# Hendrix Genetics Species Configuration Summary for Phaseimpute

## Genomic Characteristics and Configuration Parameters

| Species | Scientific Name                         | Genome Size | Autosomes (chr_set) | Avg Chr Size |
|---------|-----------------                        |-------------|---------------------|--------------|
| **Layer (Chicken)** | *Gallus gallus*             | 1.05 Gb     | 38                  | ~27 Mb |
| **Broiler (Chicken)** | *Gallus gallus*           | 1.05 Gb     | 38                  | ~27 Mb | 
| **Turkey** | *Meleagris gallopavo*                | 1.1 Gb      | 39                  | ~28 Mb | 
| **Pig** | *Sus scrofa*                            | 2.5 Gb      | 18                  | ~130 Mb |
| **Atlantic Salmon** | *Salmo salar*               | 3.0 Gb      | 29                  | ~100 Mb | 
| **Rainbow Trout** | *Oncorhynchus mykiss*         | 2.4 Gb      | 29                  | ~80 Mb | 
| **Pacific White Shrimp** | *Litopenaeus vannamei* | 1.9-2.0 Gb  | 44                  | ~42 Mb |

## Configuration Files Created

1. `species_layer.config` - Layer chickens
2. `species_broiler.config` - Broiler chickens
3. `species_turkey.config` - Turkeys
4. `species_pig.config` - Pigs
5. `species_salmon.config` - Atlantic salmon
6. `species_trout.config` - Rainbow trout
7. `species_shrimp.config` - Pacific white shrimp

## Key Parameter Rationale

### chr_set
Number of autosomal chromosomes. Required for PLINK and other tools that need to know the chromosome count.

### MINIMAC4 Chunk/Overlap
- Chunk size scales with average chromosome size (~35-40% of avg chr)
- Overlap is typically 10% of chunk size
- Smaller chromosomes (poultry microchromosomes) need smaller chunks

### GLIMPSE Window Size
- Window size scales with chromosome size
- Poultry: 1 Mb (small microchromosomes)
- Pig: 5 Mb (large chromosomes)
- Fish: 4 Mb (medium-large chromosomes)
- Shrimp: 2 Mb (medium, uniform chromosomes)

### Beagle ne (Effective Population Size)
- Poultry/Pig: ne=50 (highly selected commercial populations)
- Fish: ne=100 (larger breeding programs, more diversity)
- Shrimp: ne=200 (higher diversity, less selection history)

## Species-Specific Notes

### Poultry (Chicken & Turkey)
- Have microchromosomes (3.5-23 Mb) which are gene-rich
- May benefit from chromosome-specific chunk sizes for optimal imputation
- ZW sex determination system

### Salmonids (Salmon & Trout)
- Underwent whole genome duplication ~80 MYA
- Ongoing rediploidization creates duplicated homeologous regions
- May require specialized handling for duplicated regions
- Variable chromosome number between strains

### Shrimp
- Very high repeat content (50-64% transposable elements)
- Relatively uniform chromosome sizes
- ZW sex determination but no distinct sex chromosomes
- Consider increased iterations for imputation due to repeat complexity

**IMPORTANT: Scaffold-level assembly only!**
- The current reference genome (ASM378908v1) is **scaffold-level**, not chromosome-level
- Uses NCBI scaffold IDs: `NW_020868286.1`, `NW_020868287.1`, etc.
- Your reference panel VCF **must use the same scaffold IDs** for imputation to work
- No `_chr` variant available (scaffold names don't follow chr1/1 convention)
- **TODO**: Replace with chromosome-level assembly when one becomes available (check NCBI/Ensembl periodically)

## Reference Genomes

Pre-configured reference genomes are available for each species via the `--genome` parameter:

| Species | `--species` | `--genome` | Assembly | Source |
|---------|-------------|------------|----------|--------|
| Chicken (layer/broiler) | `layer` / `broiler` | `GRCg7b` | bGalGal1.mat.broiler.GRCg7b | Ensembl 114 |
| Turkey | `turkey` | `Turkey_5.1` | Turkey_5.1 | Ensembl 114 |
| Pig | `pig` | `Sscrofa11.1` | Sscrofa11.1 | Ensembl 114 |
| Atlantic Salmon | `salmon` | `Ssal_v3.1` | Ssal_v3.1 | Ensembl 114 |
| Rainbow Trout | `trout` | `USDA_OmykA_1.1` | USDA_OmykA_1.1 | Ensembl 114 |
| Pacific White Shrimp | `shrimp` | `ASM378908v1` | ASM378908v1 (scaffold-level*) | NCBI |

### Automatic Reference Genome Download

Reference genomes are **automatically downloaded and prepared** on first use:

1. Pipeline downloads genome from Ensembl/NCBI
2. Converts from gzip to uncompressed format
3. Optionally adds `chr` prefix (for `_chr` variants)
4. Creates FAI index
5. Caches in `genomes_cache/` directory for future runs

**No manual setup required!** Just use `--genome Sscrofa11.1` and the pipeline handles the rest.

The default cache location is `./genomes_cache`. Override with:

```bash
--genomes_cache /path/to/shared/cache
```

For shared HPC environments, point this to a shared directory so all users benefit from cached genomes.

## Usage

Specify both `--species` (for imputation parameters) and `--genome` (for reference genome):

```bash
nextflow run nf-core/phaseimpute \
  -profile singularity \
  --species pig \
  --genome Sscrofa11.1 \
  --input samplesheet.csv \
  --panel panel.csv \
  --outdir results
```

### Using a custom reference genome

If you need to use a different reference genome build, use `--fasta` instead of `--genome`:

```bash
nextflow run nf-core/phaseimpute \
  -profile singularity \
  --species pig \
  --fasta /path/to/custom_reference.fa \
  --input samplesheet.csv \
  --panel panel.csv \
  --outdir results
```

### Available parameters

**Species options** (`--species`):
- `layer` - Layer chickens
- `broiler` - Broiler chickens  
- `turkey` - Turkeys
- `pig` - Pigs
- `salmon` - Atlantic salmon
- `trout` - Rainbow trout
- `shrimp` - Pacific white shrimp

**Genome options** (`--genome`):
- `GRCg7b` - Chicken
- `Turkey_5.1` - Turkey
- `Sscrofa11.1` - Pig
- `Ssal_v3.1` - Atlantic Salmon
- `USDA_OmykA_1.1` - Rainbow Trout
- `ASM378908v1` - Pacific White Shrimp

The `--species` parameter automatically configures:
- `chr_set` for PLINK compatibility
- MINIMAC4 chunk and overlap sizes
- GLIMPSE/GLIMPSE2 window sizes
- Beagle effective population size (ne)

## Resource Allocation

Memory and CPU requirements depend on **reference panel size**, not species. A separate config file is provided for resource tuning:

```bash
# For large panels (> 30M variants), include the resources config:
nextflow run nf-core/phaseimpute \
  -profile singularity \
  -c conf/resources.config \
  --species pig \
  --genome Sscrofa11.1 \
  ...
```

Edit `conf/resources.config` to adjust memory/CPU based on your panel size:

| Panel Size | Memory Recommendation |
|------------|----------------------|
| < 10M variants | Default settings (no extra config needed) |
| 10-30M variants | Consider increasing memory |
| > 30M variants | Use high memory settings in resources.config |

Key processes that scale with panel size:
- `BCFTOOLS_PLUGINSPLIT` - splits multi-sample VCF
- `VCFLIB_VCFFIXUP` - recalculates AC/AN fields
- `SHAPEIT5_PHASECOMMON` - phasing
- `BEAGLE5_BEAGLE` - imputation
