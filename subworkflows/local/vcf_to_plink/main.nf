/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VCF → PLINK (BED + PED)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PLINK_VCF } from '/lustre/backup/HG/burgh012/phaseimpute/modules/nf-core/plink/vcf'

process PLINK_MAKE_PED {

    tag "$meta.id"

    module 'plink/1.9-180913'

    input:
    tuple val(meta), path(bed), path(bim), path(fam)

    output:
    tuple val(meta), path("*.ped"), path("*.map"), emit: ped_map
    path "versions.yml", emit: versions

    script:
    """
    plink \
        --bed $bed \
        --bim $bim \
        --fam $fam \
        --recode \
        --out ${meta.id}

    plink --version > versions.yml
    """
}

workflow VCF_TO_PLINK {

    take:
    ch_vcf   // [ meta, vcf, index ]

    main:

    // Prepare for nf-core module
    ch_plink_input = ch_vcf.map { meta, vcf, index ->
        [ meta, vcf ]
    }

    // Step 1: VCF → BED/BIM/FAM
    PLINK_VCF(ch_plink_input)

    // Combine outputs
    ch_binary = PLINK_VCF.out.bed
        .join(PLINK_VCF.out.bim)
        .join(PLINK_VCF.out.fam)

    // Step 2: BED → PED/MAP
    PLINK_MAKE_PED(ch_binary)

    emit:
bed = PLINK_VCF.out.bed
bim = PLINK_VCF.out.bim
fam = PLINK_VCF.out.fam

ped = PLINK_MAKE_PED.out.ped_map.map { it[1] }
map = PLINK_MAKE_PED.out.ped_map.map { it[2] }

versions = PLINK_MAKE_PED.out.versions
}
