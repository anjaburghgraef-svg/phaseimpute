process PER_CHROMOSOME_STATS {
    tag "${meta.id}_${meta.chr}_${meta.tools}"
    label 'process_medium'

    conda "bioconda::bcftools=1.20"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0' :
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(imputed_vcf), path(imputed_idx)
    tuple val(meta_target), path(target_vcf), path(target_idx)
    tuple val(meta_truth), path(truth_vcf), path(truth_idx)

    output:
    path "chr_stats_*.csv", emit: stats
    path "versions.yml"   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def chr = meta.chr
    def tool = meta.tools ?: "unknown"
    def has_truth = truth_vcf && !truth_vcf.name.contains("NO_FILE")
    def has_target = target_vcf && !target_vcf.name.contains("NO_FILE")
    // Sample every Nth variant for detailed metrics (R2, non-ref concordance)
    def sample_rate = 100
    """
    #!/bin/bash
    set -euo pipefail

    CHR_NAME="${chr}"
    TOOL_NAME="${tool}"
    HAS_TRUTH=${has_truth ? 'true' : 'false'}
    HAS_TARGET=${has_target ? 'true' : 'false'}
    SAMPLE_RATE=${sample_rate}

    echo "Processing chromosome \$CHR_NAME for tool \$TOOL_NAME" >&2

    # Count SNPs in imputed VCF (already per-chromosome)
    SNPS_AFTER=\$(bcftools index --nrecords "${imputed_vcf}" 2>/dev/null || echo "0")
    echo "Imputed SNPs for chr \$CHR_NAME: \$SNPS_AFTER" >&2

    # Count SNPs in target VCF for this chromosome
    SNPS_BEFORE="NA"
    if [ "\$HAS_TARGET" = "true" ]; then
        SNPS_BEFORE=\$(bcftools view -H -r "\$CHR_NAME" "${target_vcf}" 2>/dev/null | wc -l || echo "0")
        if [ "\$SNPS_BEFORE" = "0" ]; then
            if [[ "\$CHR_NAME" == chr* ]]; then
                ALT_CHR="\${CHR_NAME#chr}"
            else
                ALT_CHR="chr\$CHR_NAME"
            fi
            SNPS_BEFORE=\$(bcftools view -H -r "\$ALT_CHR" "${target_vcf}" 2>/dev/null | wc -l || echo "0")
        fi
        echo "Target SNPs for chr \$CHR_NAME: \$SNPS_BEFORE" >&2
    fi

    # Initialize accuracy metrics as NA
    GENERAL_CONCORDANCE="NA"
    DOSAGE_R2="NA"
    BEST_GT_R2="NA"
    NON_REF_CONCORDANCE="NA"

    if [ "\$HAS_TRUTH" = "true" ]; then
        echo "Calculating accuracy metrics for chr \$CHR_NAME..." >&2

        # Determine which chromosome name works for truth VCF
        TRUTH_CHR="\$CHR_NAME"
        TRUTH_COUNT=\$(bcftools view -H -r "\$CHR_NAME" "${truth_vcf}" 2>/dev/null | head -1 | wc -l || echo "0")
        if [ "\$TRUTH_COUNT" = "0" ]; then
            if [[ "\$CHR_NAME" == chr* ]]; then
                TRUTH_CHR="\${CHR_NAME#chr}"
            else
                TRUTH_CHR="chr\$CHR_NAME"
            fi
            echo "Using alternate chromosome name for truth: \$TRUTH_CHR" >&2
        fi

        # Get common samples
        IMPUTED_SAMPLES=\$(bcftools query -l "${imputed_vcf}" | sort)
        TRUTH_SAMPLES=\$(bcftools query -l "${truth_vcf}" | sort)
        COMMON_SAMPLES=\$(comm -12 <(echo "\$IMPUTED_SAMPLES") <(echo "\$TRUTH_SAMPLES"))
        N_COMMON=\$(echo "\$COMMON_SAMPLES" | grep -c . || echo "0")

        if [ "\$N_COMMON" -eq 0 ]; then
            echo "WARNING: No common samples between imputed and truth" >&2
        else
            echo "Common samples: \$N_COMMON" >&2
            echo "\$COMMON_SAMPLES" > common_samples.txt

            # Extract chromosome from truth with renamed chr if needed
            bcftools view -r "\$TRUTH_CHR" -S common_samples.txt -v snps -m2 -M2 -Oz -o truth_chr.vcf.gz "${truth_vcf}"
            bcftools index -t truth_chr.vcf.gz

            # Rename chromosome in truth to match imputed if different
            if [ "\$TRUTH_CHR" != "\$CHR_NAME" ]; then
                echo "\$TRUTH_CHR \$CHR_NAME" > chr_rename.txt
                bcftools annotate --rename-chrs chr_rename.txt -Oz -o truth_chr_renamed.vcf.gz truth_chr.vcf.gz
                mv truth_chr_renamed.vcf.gz truth_chr.vcf.gz
                bcftools index -t truth_chr.vcf.gz
            fi

            # Extract chromosome from imputed
            bcftools view -r "\$CHR_NAME" -S common_samples.txt -v snps -m2 -M2 -Oz -o imputed_chr.vcf.gz "${imputed_vcf}"
            bcftools index -t imputed_chr.vcf.gz

            # Find intersection
            bcftools isec -n=2 -w1 -Oz -o imputed_isec.vcf.gz imputed_chr.vcf.gz truth_chr.vcf.gz
            bcftools isec -n=2 -w2 -Oz -o truth_isec.vcf.gz imputed_chr.vcf.gz truth_chr.vcf.gz

            N_ISEC=\$(bcftools view -H imputed_isec.vcf.gz | wc -l)
            echo "Found \$N_ISEC intersecting variants, sampling every \${SAMPLE_RATE}th for detailed metrics..." >&2

            if [ "\$N_ISEC" -gt 0 ]; then
                # Check if DS field exists
                HAS_DS=\$(bcftools view -h imputed_isec.vcf.gz | grep -c "ID=DS," || echo "0")

                # Extract genotypes with sampling - only every Nth line
                if [ "\$HAS_DS" -gt 0 ]; then
                    bcftools query -f '[%GT:%DS\\t]\\n' imputed_isec.vcf.gz | awk -v rate="\$SAMPLE_RATE" 'NR % rate == 1' > imputed_gt.txt
                else
                    bcftools query -f '[%GT\\t]\\n' imputed_isec.vcf.gz | awk -v rate="\$SAMPLE_RATE" 'NR % rate == 1' > imputed_gt.txt
                fi
                bcftools query -f '[%GT\\t]\\n' truth_isec.vcf.gz | awk -v rate="\$SAMPLE_RATE" 'NR % rate == 1' > truth_gt.txt

                SAMPLED_VARIANTS=\$(wc -l < truth_gt.txt)
                echo "Calculating metrics from \$SAMPLED_VARIANTS sampled variants..." >&2

                # Calculate all metrics using sampled data
                METRICS=\$(paste truth_gt.txt imputed_gt.txt | awk -v has_ds="\$HAS_DS" '
                function parse_gt(gt) {
                    gsub(/\\|/, "/", gt)
                    if (gt == "." || gt == "./.") return -1
                    split(gt, a, "/")
                    if (a[1] == "." || a[2] == ".") return -1
                    g1 = int(a[1]); g2 = int(a[2])
                    if (g1 < 0 || g1 > 1 || g2 < 0 || g2 > 1) return -1
                    return g1 + g2
                }
                BEGIN {
                    total=0; match_total=0
                    rr_total=0; rr_match=0
                    ra_total=0; ra_match=0
                    aa_total=0; aa_match=0
                    gt_n=0; gt_sx=0; gt_sy=0; gt_sx2=0; gt_sy2=0; gt_sxy=0
                    ds_n=0; ds_sx=0; ds_sy=0; ds_sx2=0; ds_sy2=0; ds_sxy=0
                }
                {
                    n = NF / 2
                    for (i = 1; i <= n; i++) {
                        t_gt = parse_gt(\$i)
                        imp_field = \$(i + n)
                        if (has_ds > 0 && index(imp_field, ":") > 0) {
                            split(imp_field, parts, ":")
                            i_gt = parse_gt(parts[1])
                            i_ds = (parts[2] == "." || parts[2] == "") ? -1 : parts[2] + 0
                        } else {
                            i_gt = parse_gt(imp_field)
                            i_ds = -1
                        }
                        if (t_gt < 0 || i_gt < 0) continue

                        # General concordance
                        total++
                        if (t_gt == i_gt) match_total++

                        # Per-genotype concordance
                        if (t_gt == 0) { rr_total++; if (i_gt == 0) rr_match++ }
                        else if (t_gt == 1) { ra_total++; if (i_gt == 1) ra_match++ }
                        else if (t_gt == 2) { aa_total++; if (i_gt == 2) aa_match++ }

                        # R2 accumulators
                        gt_n++; gt_sx += t_gt; gt_sy += i_gt
                        gt_sx2 += t_gt*t_gt; gt_sy2 += i_gt*i_gt; gt_sxy += t_gt*i_gt

                        ds_val = (i_ds >= 0) ? i_ds : i_gt
                        ds_n++; ds_sx += t_gt; ds_sy += ds_val
                        ds_sx2 += t_gt*t_gt; ds_sy2 += ds_val*ds_val; ds_sxy += t_gt*ds_val
                    }
                }
                END {
                    # General concordance
                    if (total > 0) gen_conc = match_total / total
                    else gen_conc = "NA"

                    # Non-ref concordance (RA + AA)
                    non_ref_total = ra_total + aa_total
                    non_ref_match = ra_match + aa_match
                    if (non_ref_total > 0) non_ref_conc = non_ref_match / non_ref_total
                    else non_ref_conc = "NA"

                    # Best GT R2
                    if (gt_n >= 2) {
                        cov = gt_sxy - (gt_sx * gt_sy / gt_n)
                        vx = gt_sx2 - (gt_sx * gt_sx / gt_n)
                        vy = gt_sy2 - (gt_sy * gt_sy / gt_n)
                        if (vx > 0 && vy > 0) gt_r2 = (cov / sqrt(vx * vy))^2
                        else gt_r2 = "NA"
                    } else gt_r2 = "NA"

                    # Dosage R2
                    if (ds_n >= 2) {
                        cov = ds_sxy - (ds_sx * ds_sy / ds_n)
                        vx = ds_sx2 - (ds_sx * ds_sx / ds_n)
                        vy = ds_sy2 - (ds_sy * ds_sy / ds_n)
                        if (vx > 0 && vy > 0) ds_r2 = (cov / sqrt(vx * vy))^2
                        else ds_r2 = "NA"
                    } else ds_r2 = "NA"

                    printf "%.6f\\t%.6f\\t%.6f\\t%.6f", gen_conc, ds_r2, gt_r2, non_ref_conc
                }')

                GENERAL_CONCORDANCE=\$(echo "\$METRICS" | cut -f1)
                DOSAGE_R2=\$(echo "\$METRICS" | cut -f2)
                BEST_GT_R2=\$(echo "\$METRICS" | cut -f3)
                NON_REF_CONCORDANCE=\$(echo "\$METRICS" | cut -f4)

                echo "General concordance: \$GENERAL_CONCORDANCE" >&2
                echo "Dosage R2: \$DOSAGE_R2, Best GT R2: \$BEST_GT_R2" >&2
                echo "Non-ref concordance: \$NON_REF_CONCORDANCE" >&2
            fi

            # Cleanup
            rm -f truth_chr.vcf.gz* imputed_chr.vcf.gz* common_samples.txt chr_rename.txt
            rm -f imputed_isec.vcf.gz truth_isec.vcf.gz imputed_gt.txt truth_gt.txt
        fi
    fi

    # Write output CSV
    echo "chr,tool,snps_before_imputation,snps_after_imputation,general_concordance,dosage_r2,best_gt_r2,non_ref_concordance" > "chr_stats_\${CHR_NAME}_\${TOOL_NAME}.csv"
    echo "\$CHR_NAME,\$TOOL_NAME,\$SNPS_BEFORE,\$SNPS_AFTER,\$GENERAL_CONCORDANCE,\$DOSAGE_R2,\$BEST_GT_R2,\$NON_REF_CONCORDANCE" >> "chr_stats_\${CHR_NAME}_\${TOOL_NAME}.csv"

    echo "Per-chromosome stats written for \$CHR_NAME" >&2

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version | head -n1 | sed 's/^bcftools //')
    END_VERSIONS
    """
}
