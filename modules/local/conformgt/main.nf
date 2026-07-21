process CONFORMGT {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://eclipse-temurin:8-jre' :
        'eclipse-temurin:8-jre' }"

    input:
    tuple val(meta), path(target_vcf), path(target_index)
    path(panel_vcfs)
    path(panel_indices)
    val(chromosomes)
    path(conformgt_jar)

    output:
    tuple val(meta), path("${meta.id}.conformed.vcf.gz"), path("${meta.id}.conformed.vcf.gz.tbi"), emit: vcf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def match_mode = task.ext.args ?: "POS"
    def chrom_list = chromosomes.join(' ')
    """
    #!/bin/bash
    set -euo pipefail

    echo "Starting CONFORMGT process" >&2
    echo "Java version:" >&2
    java -version 2>&1 || { echo "ERROR: Java not found" >&2; exit 1; }
    echo "Working directory: \$(pwd)" >&2
    echo "Panel files found:" >&2
    ls -la *.vcf.gz 2>&1 | head -20 >&2

    # Process each chromosome
    processed=0
    for chr in ${chrom_list}; do
        # Find the panel file for this chromosome
        # Try multiple patterns: *_chrX_*.vcf.gz, *_chrX.vcf.gz, *chrX.vcf.gz
        panel_file=\$(ls -1 *_\${chr}_*.vcf.gz 2>/dev/null | grep -v conformed | grep -v tmp | head -1 || true)
        if [ -z "\$panel_file" ]; then
            panel_file=\$(ls -1 *_\${chr}.vcf.gz 2>/dev/null | grep -v conformed | grep -v tmp | head -1 || true)
        fi
        if [ -z "\$panel_file" ]; then
            panel_file=\$(ls -1 *\${chr}.vcf.gz 2>/dev/null | grep -v conformed | grep -v tmp | grep -v "^\${chr}" | head -1 || true)
        fi

        if [ -n "\$panel_file" ]; then
            echo "Processing \$chr with panel \$panel_file" >&2

            # Remove any existing output to avoid "file exists" error
            rm -f ${meta.id}.\$chr.tmp.vcf.gz ${meta.id}.\$chr.tmp.log

            java -jar ${conformgt_jar} \\
                ref=\$panel_file \\
                gt=${target_vcf} \\
                match=${match_mode} \\
                chrom=\$chr \\
                out=${meta.id}.\$chr.tmp 2>&1 || echo "WARNING: conform-gt failed for \$chr, continuing..." >&2

            if [ -f "${meta.id}.\$chr.tmp.vcf.gz" ]; then
                processed=\$((processed + 1))
            fi
        else
            echo "WARNING: No panel file found for \$chr" >&2
        fi
    done

    echo "Processed \$processed chromosomes" >&2

    # Concatenate all chromosome outputs into ONE file
    if ls *.tmp.vcf.gz 1>/dev/null 2>&1; then
        echo "Validating and indexing output files..." >&2
        # Validate each tmp file is a proper BGZF file before indexing
        valid_files=""
        for f in *.tmp.vcf.gz; do
            if bcftools view -h "\$f" >/dev/null 2>&1; then
                bcftools index -t \$f
                valid_files="\$valid_files \$f"
            else
                echo "WARNING: Corrupted or incomplete file \$f (possibly OOM killed), skipping..." >&2
                rm -f "\$f"
            fi
        done

        if [ -z "\$valid_files" ]; then
            echo "ERROR: No valid output files after validation" >&2
            echo "This usually means the process ran out of memory." >&2
            echo "Try increasing memory allocation for CONFORMGT process." >&2
            exit 1
        fi

        echo "Concatenating valid files:\$valid_files" >&2
        bcftools concat -Oz -o ${meta.id}.conformed.vcf.gz \$valid_files
        bcftools index -t ${meta.id}.conformed.vcf.gz
        rm -f *.tmp.vcf.gz *.tmp.vcf.gz.tbi *.tmp.log
        echo "Successfully created ${meta.id}.conformed.vcf.gz" >&2
    else
        echo "ERROR: No output files created after processing \$processed chromosomes" >&2
        echo "Files in directory:" >&2
        ls -la >&2
        exit 1
    fi

    echo "${task.process}:" > versions.yml
    echo "    conformgt: r1174" >> versions.yml
    echo "    java: 8" >> versions.yml
    """

    stub:
    """
    echo "" | gzip > ${meta.id}.conformed.vcf.gz
    touch ${meta.id}.conformed.vcf.gz.tbi

    echo "${task.process}:" > versions.yml
    echo "    conformgt: r1174" >> versions.yml
    echo "    java: 8" >> versions.yml
    """
}
