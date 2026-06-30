include { BCFTOOLS_MERGE as BCFTOOLS_MERGE_CONFORMGT } from '../../../modules/nf-core/bcftools/merge'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_MERGED    } from '../../../modules/nf-core/bcftools/index'
include { CONFORMGT                                  } from '../../../modules/local/conformgt'

workflow VCF_CONFORM_GT {
    take:
    ch_input  // channel: [ [id, ...], vcf, tbi ] - target VCF(s), may be multiple single-sample files
    ch_panel  // channel: [ [id, chr], vcf, index ] - panel per chromosome
    ch_jar    // channel: conform-gt.jar file

    main:
    ch_versions = Channel.empty()

    // Convert JAR channel to a value channel so it can be reused
    ch_jar_value = ch_jar.first()

    // Collect all input VCFs and count them
    ch_vcfs_list = ch_input
        .map { meta, vcf, idx -> [vcf, idx] }
        .toList()

    // Prepare inputs for merge (needs at least 2 VCFs)
    // If only 1 VCF, we'll skip merge and use it directly
    ch_for_merge = ch_vcfs_list
        .filter { items -> items.size() > 1 }
        .map { items ->
            def vcfs = items.collect { it[0] }
            def indices = items.collect { it[1] }
            [[id: "MERGED_TARGET", tools: [], batch: 0], vcfs, indices, []]
        }

    // Single VCF case - just pass through with standardized metadata
    ch_single = ch_vcfs_list
        .filter { items -> items.size() == 1 }
        .map { items ->
            def vcf = items[0][0]
            def idx = items[0][1]
            [[id: "MERGED_TARGET", tools: [], batch: 0], vcf, idx]
        }

    // Merge multiple VCFs into one
    BCFTOOLS_MERGE_CONFORMGT(
        ch_for_merge,
        Channel.value([[], [], []])  // empty fasta reference as value channel
    )

    // Index the merged VCF (only runs if merge happened)
    BCFTOOLS_INDEX_MERGED(BCFTOOLS_MERGE_CONFORMGT.out.vcf)

    // Combine merged output with its index
    ch_merged = BCFTOOLS_MERGE_CONFORMGT.out.vcf
        .join(BCFTOOLS_INDEX_MERGED.out.tbi.mix(BCFTOOLS_INDEX_MERGED.out.csi))
        .map { meta, vcf, idx ->
            [[id: meta.id, tools: [], batch: 0], vcf, idx]
        }

    // Mix versions from merge/index (only if they ran - empty channels just don't emit)
    ch_versions = ch_versions.mix(BCFTOOLS_MERGE_CONFORMGT.out.versions)
    ch_versions = ch_versions.mix(BCFTOOLS_INDEX_MERGED.out.versions)

    // Final target: either merged multi-sample or single original VCF
    ch_target = ch_merged.mix(ch_single)

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
        ch_jar_value
    )
    ch_versions = ch_versions.mix(CONFORMGT.out.versions)

    emit:
    vcf_tbi   = CONFORMGT.out.vcf  // channel: [ [id: MERGED_TARGET], vcf, tbi ] - ONE merged harmonized file
    versions  = ch_versions         // channel: [ versions.yml ]
}
