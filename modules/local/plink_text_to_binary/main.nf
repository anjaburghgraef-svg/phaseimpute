// Module: PLINK_TEXT_TO_BINARY
// Converts PLINK text format (.ped/.map) to PLINK binary format (.bed/.bim/.fam)
// To be placed in: modules/local/plink_text_to_binary/main.nf

process PLINK_TEXT_TO_BINARY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1':
        'quay.io/biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
    tuple val(meta), path(ped), path(map)

    output:
    tuple val(meta), path("*.bed"), path("*.bim"), path("*.fam"), emit: plink
    path "versions.yml"                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def chr_set = params.chr_set ? "--chr-set ${params.chr_set}" : ""

    """
    plink \\
        --ped ${ped} \\
        --map ${map} \\
        --make-bed \\
        --allow-extra-chr \\
        ${chr_set} \\
        --out ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(plink --version 2>&1 | head -1 | sed 's/PLINK v//; s/ .*//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed
    touch ${prefix}.bim
    touch ${prefix}.fam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plink: \$(plink --version 2>&1 | head -1 | sed 's/PLINK v//; s/ .*//')
    END_VERSIONS
    """
}
