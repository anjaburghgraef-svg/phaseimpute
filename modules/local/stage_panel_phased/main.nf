process STAGE_PANEL_PHASED {
    tag "${meta.id} ${meta.chr}"

    publishDir "${params.outdir}/prep_panel/panel", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(vcf), path(csi)

    output:
    tuple val(meta),
          path("${vcf.getName()}"),
          path("${csi.getName()}"),
          emit: vcf_tbi
    path "versions.yml", emit: versions

    script:
    """
    cat <<-END_VERSIONS > versions.yml
    "NFCORE_PHASEIMPUTE:PHASEIMPUTE:STAGE_PANEL_PHASED":
        stage_panel_phased: "1"
    END_VERSIONS
    """
}