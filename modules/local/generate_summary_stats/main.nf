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
    get_vcf_stats() {
        local pattern=\$1
        local vcf_files=(\$(ls \${pattern} 2>/dev/null | grep -v "NO_FILE"))

        if [ \${#vcf_files[@]} -eq 0 ]; then
            echo "0\t0"
            return
        fi

        # Get unique sample names across all VCFs
        local sample_file="samples_\${RANDOM}.txt"
        for vcf in "\${vcf_files[@]}"; do
            bcftools query -l "\$vcf" 2>/dev/null >> "\$sample_file"
        done

        local samples=\$(sort -u "\$sample_file" 2>/dev/null | wc -l)
        samples=\${samples:-0}
        rm -f "\$sample_file"

        # Deduplicate VCF files by resolving symlinks to real paths
        declare -A seen_paths
        local unique_vcfs=()

        for vcf in "\${vcf_files[@]}"; do
            local real_path=\$(readlink -f "\$vcf" 2>/dev/null || realpath "\$vcf" 2>/dev/null || echo "\$vcf")
            if [ -z "\${seen_paths[\$real_path]}" ]; then
                seen_paths[\$real_path]=1
                unique_vcfs+=("\$real_path")
            fi
        done

        # Sum SNPs across unique VCFs only
        local total_snps=0
        for vcf in "\${unique_vcfs[@]}"; do
            local snps=\$(bcftools index --nrecords "\$vcf" 2>/dev/null)
            if [ -z "\$snps" ] || [ "\$snps" = "0" ]; then
                snps=\$(bcftools view -H "\$vcf" 2>/dev/null | wc -l)
            fi
            snps=\${snps:-0}
            total_snps=\$((total_snps + snps))
        done

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
