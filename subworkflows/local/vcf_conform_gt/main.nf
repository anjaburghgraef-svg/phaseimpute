include { CONFORMGT } from '../../../modules/local/conformgt'

workflow VCF_CONFORM_GT {
    take:
    ch_input  // channel: [ [id, ...], vcf, tbi ] - target VCF (single file, all chromosomes)
    ch_panel  // channel: [ [id, chr], vcf, index ] - panel per chromosome
    ch_jar    // channel: conform-gt.jar file

    main:
    ch_versions = Channel.empty()

    // Get unique target (without chr-specific metadata)
    ch_target = ch_input
        .map { meta, vcf, index ->
            def base_meta = [id: meta.id, tools: meta.tools ?: [], batch: meta.batch ?: 0]
            [base_meta, vcf, index]
        }
        .unique()

    // Collect all panel VCFs and their chromosomes
    ch_panel_collected = ch_panel
        .map { meta, vcf, index -> [vcf, index, meta.chr] }
        .toList()
        .map { items ->
            def vcfs = items.collect { it[0] }
            def indices = items.collect { it[1] }
            def chroms = items.collect { it[2] }
            [vcfs, indices, chroms]
        }

    // Run conform-gt harmonization (single task, all chromosomes, outputs ONE concatenated file)
    CONFORMGT(
        ch_target,
        ch_panel_collected.map { it[0] },  // panel VCFs
        ch_panel_collected.map { it[1] },  // panel indices
        ch_panel_collected.map { it[2] },  // chromosome list
        ch_jar
    )
    ch_versions = ch_versions.mix(CONFORMGT.out.versions.first())

    emit:
    vcf_tbi   = CONFORMGT.out.vcf  // channel: [ [id], vcf, tbi ] - ONE file with all chromosomes
    versions  = ch_versions         // channel: [ versions.yml ]
}
