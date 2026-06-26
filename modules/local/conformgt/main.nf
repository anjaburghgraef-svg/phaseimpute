process CONFORMGT {
    tag "$meta.id"
    label 'process_medium'

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
    # Process each chromosome
    for chr in ${chrom_list}; do
        # Find the panel file for this chromosome (match _chrX_ exactly to avoid chr1 matching chr10)
        panel_file=\$(ls -1 *_\${chr}_*.vcf.gz 2>/dev/null | grep -v conformed | grep -v tmp | head -1)

        if [ -n "\$panel_file" ]; then
            echo "Processing \$chr with panel \$panel_file" >&2

            java -jar ${conformgt_jar} \\
                ref=\$panel_file \\
                gt=${target_vcf} \\
                match=${match_mode} \\
                chrom=\$chr \\
                out=${meta.id}.\$chr.tmp || true
        else
            echo "WARNING: No panel file found for \$chr" >&2
        fi
    done

    # Concatenate all chromosome outputs into ONE file
    if ls *.tmp.vcf.gz 1>/dev/null 2>&1; then
        # Index each tmp file first
        for f in *.tmp.vcf.gz; do
            bcftools index -t \$f
        done
        bcftools concat -Oz -o ${meta.id}.conformed.vcf.gz *.tmp.vcf.gz
        bcftools index -t ${meta.id}.conformed.vcf.gz
        rm -f *.tmp.vcf.gz *.tmp.vcf.gz.tbi *.tmp.log
    else
        echo "ERROR: No output files created" >&2
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
