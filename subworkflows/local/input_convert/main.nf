//
// SUBWORKFLOW: INPUT_CONVERT
// Converts PLINK and Mixblup format inputs to VCF for use in the phaseimpute pipeline
//
// To be placed in: subworkflows/local/input_convert/main.nf
//

include { PLINK_TO_VCF         } from '../../../modules/local/plink_to_vcf/main'
include { PLINK_TEXT_TO_BINARY } from '../../../modules/local/plink_text_to_binary/main'
include { MIXBLUP_TO_PLINK     } from '../../../modules/local/mixblup_to_plink/main'
include { TABIX_TABIX          } from '../../../modules/nf-core/tabix/tabix/main'

workflow INPUT_CONVERT {

    take:
    ch_input_plink_binary  // channel: [ [meta], bed, bim, fam ]
    ch_input_plink_text    // channel: [ [meta], ped, map ]
    ch_input_mixblup       // channel: [ [meta], gtp, ped, map_mix, snp_details ]

    main:
    ch_versions = Channel.empty()

    //
    // Convert PLINK text -> PLINK binary
    //
    PLINK_TEXT_TO_BINARY(ch_input_plink_text)
    ch_versions = ch_versions.mix(PLINK_TEXT_TO_BINARY.out.versions.first())

    //
    // Convert Mixblup -> PLINK binary
    //
    MIXBLUP_TO_PLINK(ch_input_mixblup)
    ch_versions = ch_versions.mix(MIXBLUP_TO_PLINK.out.versions.first())

    //
    // Combine all PLINK binary inputs
    //
    ch_all_plink = ch_input_plink_binary
        .mix(PLINK_TEXT_TO_BINARY.out.plink)
        .mix(MIXBLUP_TO_PLINK.out.plink)

    //
    // Convert PLINK binary -> VCF
    //
    PLINK_TO_VCF(ch_all_plink)
    ch_versions = ch_versions.mix(PLINK_TO_VCF.out.versions.first())

    //
    // Index the VCF with tabix
    //
    TABIX_TABIX(PLINK_TO_VCF.out.vcf)
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions.first())

    // Combine VCF with its index
    ch_vcf_tbi = PLINK_TO_VCF.out.vcf
        .join(TABIX_TABIX.out.tbi)

    emit:
    vcf_tbi  = ch_vcf_tbi    // channel: [ [meta], vcf.gz, vcf.gz.tbi ]
    versions = ch_versions   // channel: [ versions.yml ]
}
