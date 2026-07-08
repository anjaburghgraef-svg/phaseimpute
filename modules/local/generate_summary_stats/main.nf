process GENERATE_QC_METRICS {
    tag "qc_metrics"
    label 'process_single'

    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    path target_vcfs, stageAs: "target_*.vcf.gz"
    path truth_vcfs, stageAs: "truth_*.vcf.gz"
    path panel_vcfs, stageAs: "panel_*.vcf.gz"
    path imputed_vcfs, stageAs: "imputed_*.vcf.gz"
    path per_chr_stats

    output:
    path "summary_stats.txt", emit: summary_stats
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def has_per_chr = per_chr_stats && !per_chr_stats.name.contains("NO_FILE")
    """
    #!/usr/bin/env bash

    # Function to get sample count and SNP count from multiple VCF files
    # For samples: assumes same samples in all VCFs (per-chromosome splits)
    # For SNPs: sums across all VCF files (handles per-chromosome split panels)
    get_vcf_stats() {
        local pattern=\$1
        local sum_snps=\${2:-false}  # If true, sum SNPs across all files
        local vcf_files=(\$(ls \${pattern} 2>/dev/null | grep -v "NO_FILE"))

        if [ \${#vcf_files[@]} -eq 0 ]; then
            echo "0\t0"
            return
        fi

        local first_vcf="\${vcf_files[0]}"
        local num_files=\${#vcf_files[@]}

        # Get samples from first VCF
        local samples_in_first=\$(bcftools query -l "\$first_vcf" 2>/dev/null | wc -l)
        samples_in_first=\${samples_in_first:-0}

        # If first VCF has 1 sample and we have multiple files, assume per-sample VCFs
        # Sample count = number of files (much faster than querying each file)
        local samples
        if [ "\$samples_in_first" -eq 1 ] && [ "\$num_files" -gt 1 ]; then
            samples=\$num_files
        else
            # Multi-sample VCF(s) - need to count unique samples across all
            local sample_file="samples_\${RANDOM}.txt"
            for vcf in "\${vcf_files[@]}"; do
                bcftools query -l "\$vcf" 2>/dev/null >> "\$sample_file"
            done
            samples=\$(sort -u "\$sample_file" 2>/dev/null | wc -l)
            samples=\${samples:-0}
            rm -f "\$sample_file"
        fi

        # Get SNP count - sum across all files if requested (for panel with per-chr VCFs)
        local total_snps=0
        if [ "\$sum_snps" = "true" ]; then
            for vcf in "\${vcf_files[@]}"; do
                local file_snps=\$(bcftools index --nrecords "\$vcf" 2>/dev/null)
                if [ -z "\$file_snps" ] || [ "\$file_snps" = "0" ]; then
                    file_snps=\$(bcftools view -H "\$vcf" 2>/dev/null | wc -l)
                fi
                file_snps=\${file_snps:-0}
                total_snps=\$((total_snps + file_snps))
            done
        else
            # Just use first VCF (for target/imputed which are typically same structure)
            total_snps=\$(bcftools index --nrecords "\$first_vcf" 2>/dev/null)
            if [ -z "\$total_snps" ] || [ "\$total_snps" = "0" ]; then
                total_snps=\$(bcftools view -H "\$first_vcf" 2>/dev/null | wc -l)
            fi
            total_snps=\${total_snps:-0}
        fi

        echo "\${samples}\t\${total_snps}"
    }

    # Check if we have per-chromosome stats to use for faster SNP counting
    HAS_PER_CHR=${has_per_chr ? 'true' : 'false'}

    if [ "\$HAS_PER_CHR" = "true" ] && [ -f "${per_chr_stats}" ]; then
        echo "Using per-chromosome stats for SNP counts" >&2

        # Get the first tool name from the per-chr stats (to avoid summing across tools)
        FIRST_TOOL=\$(awk -F',' 'NR==2 {print \$2}' "${per_chr_stats}")
        echo "Using tool '\$FIRST_TOOL' for SNP counts" >&2

        # Sum imputed SNPs from per-chr stats for first tool only (column 4: snps_after_imputation)
        IMPUTED_SNPS=\$(awk -F',' -v tool="\$FIRST_TOOL" 'NR>1 && \$2==tool && \$4 != "NA" {sum += \$4} END {print sum+0}' "${per_chr_stats}")

        # Sum target SNPs from per-chr stats for first tool only (column 3: snps_before_imputation)
        TARGET_SNPS_FROM_CHR=\$(awk -F',' -v tool="\$FIRST_TOOL" 'NR>1 && \$2==tool && \$3 != "NA" {sum += \$3} END {print sum+0}' "${per_chr_stats}")

        # Get sample counts from VCFs (still need this)
        read TARGET_SAMPLES _unused <<< \$(get_vcf_stats "target_*.vcf.gz")
        read IMPUTED_SAMPLES _unused <<< \$(get_vcf_stats "imputed_*.vcf.gz")

        # Use per-chr target SNPs if available, otherwise fall back to VCF counting
        if [ "\$TARGET_SNPS_FROM_CHR" -gt 0 ]; then
            TARGET_SNPS=\$TARGET_SNPS_FROM_CHR
        else
            read _unused TARGET_SNPS <<< \$(get_vcf_stats "target_*.vcf.gz")
        fi
    else
        # Fallback to traditional VCF counting
        read TARGET_SAMPLES TARGET_SNPS <<< \$(get_vcf_stats "target_*.vcf.gz")
        read IMPUTED_SAMPLES IMPUTED_SNPS <<< \$(get_vcf_stats "imputed_*.vcf.gz")
    fi

    # Get stats for truth (always from VCF) - sum across all files in case of per-chr split
    read TRUTH_SAMPLES TRUTH_SNPS <<< \$(get_vcf_stats "truth_*.vcf.gz" "true")

    # Get stats for panel (always from VCF) - sum across all files for per-chr panels
    read PANEL_SAMPLES PANEL_SNPS <<< \$(get_vcf_stats "panel_*.vcf.gz" "true")

    # Write summary stats
    cat > summary_stats.txt <<EOF
Metric\tValue
Target_Samples\t\${TARGET_SAMPLES}
Target_SNPs\t\${TARGET_SNPS}
Truth_Samples\t\${TRUTH_SAMPLES}
Truth_SNPs\t\${TRUTH_SNPS}
Panel_Samples\t\${PANEL_SAMPLES}
Panel_SNPs\t\${PANEL_SNPS}
Imputed_Samples\t\${IMPUTED_SAMPLES}
Imputed_SNPs\t\${IMPUTED_SNPS}
EOF

    # Write versions
    cat > versions.yml <<EOF
"${task.process}":
    bcftools: \$(bcftools --version | head -n1 | sed 's/^bcftools //')
EOF
    """
}
