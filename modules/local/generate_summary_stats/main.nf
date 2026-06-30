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

    output:
    path "summary_stats.txt", emit: summary_stats
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env bash

    # Function to get sample count and SNP count from multiple VCF files
    # Optimized: only queries first VCF for SNPs, samples counted efficiently
    get_vcf_stats() {
        local pattern=\$1
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

        # Get SNP count from first VCF only
        local total_snps=\$(bcftools index --nrecords "\$first_vcf" 2>/dev/null)
        if [ -z "\$total_snps" ] || [ "\$total_snps" = "0" ]; then
            total_snps=\$(bcftools view -H "\$first_vcf" 2>/dev/null | wc -l)
        fi
        total_snps=\${total_snps:-0}

        echo "\${samples}\t\${total_snps}"
    }

    # Get stats for target
    read TARGET_SAMPLES TARGET_SNPS <<< \$(get_vcf_stats "target_*.vcf.gz")

    # Get stats for truth
    read TRUTH_SAMPLES TRUTH_SNPS <<< \$(get_vcf_stats "truth_*.vcf.gz")

    # Get stats for panel
    read PANEL_SAMPLES PANEL_SNPS <<< \$(get_vcf_stats "panel_*.vcf.gz")

    # Get stats for imputed
    read IMPUTED_SAMPLES IMPUTED_SNPS <<< \$(get_vcf_stats "imputed_*.vcf.gz")

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
