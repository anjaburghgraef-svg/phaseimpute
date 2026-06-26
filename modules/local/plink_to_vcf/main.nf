// Module: PLINK_TO_VCF
// Converts PLINK binary files (.bed/.bim/.fam) to VCF format
// To be placed in: modules/local/plink_to_vcf/main.nf

process PLINK_TO_VCF {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink2:2.00a5.12--h4ac6f70_0':
        'quay.io/biocontainers/plink2:2.00a5.12--h4ac6f70_0' }"

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("*.vcf.gz"), emit: vcf
    path "versions.yml"              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chr_set = params.chr_set ? "--chr-set ${params.chr_set}" : ""

    """
    # Convert PLINK to VCF using plink2
    plink2 \\
        --bed ${bed} \\
        --bim ${bim} \\
        --fam ${fam} \\
        ${chr_set} \\
        --allow-extra-chr \\
        --export vcf-4.2 bgz id-paste=iid \\
        --out ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | head -1 | sed 's/PLINK v//; s/ .*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink2: \$(plink2 --version 2>&1 | head -1 | sed 's/PLINK v//; s/ .*//')
    END_VERSIONS
    """
}
