include { GLIMPSE2_CHUNK                         } from '../../../modules/nf-core/glimpse2/chunk'
include { SHAPEIT5_PHASECOMMON                   } from '../../../modules/nf-core/shapeit5/phasecommon'
include { SHAPEIT5_LIGATE                        } from '../../../modules/nf-core/shapeit5/ligate'
include { BCFTOOLS_INDEX as VCF_BCFTOOLS_INDEX_1 } from '../../../modules/nf-core/bcftools/index'
include { BCFTOOLS_INDEX as VCF_BCFTOOLS_INDEX_2 } from '../../../modules/nf-core/bcftools/index'

workflow VCF_PHASE_SHAPEIT5 {

    take:
    ch_vcf        // channel (mandatory) : [ [id, chr], vcf, index, pedigree ]
    ch_region     // channel (mandatory) : [ [chr, region], region ]
    ch_ref        // channel (optional)  : [ [id, chr], vcf, index ]
    ch_scaffold   // channel (optional)  : [ [id, chr], vcf, index ]
    ch_map        // channel (optional) : [ [chr], map]
    chunk_model   // channel (mandatory) : [ model ]

    main:

    ch_versions = Channel.empty()

    /*
     * Keep region available by chromosome.
     * Needed for fallback when GLIMPSE2 writes an empty chunk file.
     */
    ch_region_by_chr = ch_region
        .map { metaCR, region -> [ metaCR.subMap("chr"), metaCR, region ] }

    // Chunk with Glimpse2
    ch_input_glimpse2 = ch_vcf
        .map{
            metaIC, vcf, csi, _pedigree -> [metaIC.subMap("chr"), metaIC, vcf, csi]
        }
        .combine(ch_region.map{ metaCR, region -> [metaCR.subMap("chr"), region]}, by:0)
        .join(ch_map)
        .map{
            _metaC, metaIC, vcf, csi, region, gmap -> [metaIC, vcf, csi, region, gmap]
        }

    GLIMPSE2_CHUNK ( ch_input_glimpse2, chunk_model )
    ch_versions = ch_versions.mix( GLIMPSE2_CHUNK.out.versions.first() )

    /*
     * GLIMPSE2 sequential can detect "one chunk" but still write an empty file
     * for some small chromosomes/regions.
     *
     * Split chunk files into:
     *  - non-empty files: use as-is
     *  - empty files: create a fallback one-line chunk file using the full region
     */
    ch_chunk_files_all = GLIMPSE2_CHUNK.out.chunk_chr
        .map { metaIC, chunkfile ->
            [ metaIC.subMap("chr"), metaIC, chunkfile ]
        }

    ch_chunk_files_nonempty = ch_chunk_files_all
        .filter { _metaC, _metaIC, chunkfile -> chunkfile.size() > 0 }
        .map { _metaC, metaIC, chunkfile -> [ metaIC, chunkfile ] }

    ch_chunk_files_empty = ch_chunk_files_all
        .filter { _metaC, _metaIC, chunkfile -> chunkfile.size() == 0 }
        .join(ch_region_by_chr, by:0)
        .map { _metaC, metaIC, chunkfile, _metaCR, region ->
            [ metaIC, region ]
        }

    MAKE_SINGLE_CHUNK_FILE(ch_chunk_files_empty)

    ch_chunk_files_used = ch_chunk_files_nonempty
        .mix(MAKE_SINGLE_CHUNK_FILE.out.chunkfile)

    // Rearrange channels
    ch_chunks_glimpse2 = ch_chunk_files_used
        .splitCsv(
            header: [
                'ID', 'Chr', 'RegionBuf', 'RegionCnk', 'WindowCm',
                'WindowMb', 'NbTotVariants', 'NbComVariants'
            ], sep: "\t", skip: 0
        )
        .map { metaIC, it -> [metaIC, it["RegionBuf"], it["RegionCnk"]]}

    ch_phase_input = ch_vcf
        .combine(ch_chunks_glimpse2, by:0)
        .map{
            metaIC, vcf, csi, pedigree, regionbuf, regioncnk -> [metaIC.subMap("chr"), metaIC, vcf, csi, pedigree, regionbuf, regioncnk]
        }
        .combine(ch_map, by:0)
        .map { _metaC, metaIC, vcf, index, pedigree, regionbuf, regioncnk, gmap ->
            [metaIC + [chunk: regioncnk], vcf, index, pedigree, regionbuf, gmap]
        }

    SHAPEIT5_PHASECOMMON (
        ch_phase_input, ch_ref,
        ch_scaffold
    )
    ch_versions = ch_versions.mix(SHAPEIT5_PHASECOMMON.out.versions.first())

    VCF_BCFTOOLS_INDEX_1(SHAPEIT5_PHASECOMMON.out.phased_variant)
    ch_versions = ch_versions.mix(VCF_BCFTOOLS_INDEX_1.out.versions.first())

    /*
     * Group phased chunk outputs by chromosome.
     * We keep the phased VCF/BCF files sorted by chunk start.
     */
    ch_grouped_phased = SHAPEIT5_PHASECOMMON.out.phased_variant
        .join(VCF_BCFTOOLS_INDEX_1.out.csi, failOnMismatch:true, failOnDuplicate:true)
        .map{ meta, vcf, csi -> [meta.subMap("id", "chr"), [vcf, meta.chunk], csi]}
        .groupTuple()
        .map{ meta, vcf, csi ->
                [ meta,
                vcf
                    .sort { a, b ->
                        def aStart = a.last().split("-")[-1].toInteger()
                        def bStart = b.last().split("-")[-1].toInteger()
                        aStart <=> bStart
                    }
                    .collect{it.first()},
                csi]}

    /*
     * Split chromosomes into:
     *  - multi-chunk  -> need ligation
     *  - single-chunk -> bypass ligation, but still must be finalized into
     *                    the normal chromosome-level *_phased.vcf.gz output
     */
    ch_multi_chunk = ch_grouped_phased
        .filter { meta, vcfs, csis -> vcfs.size() > 1 }

    ch_single_chunk = ch_grouped_phased
        .filter { meta, vcfs, csis -> vcfs.size() == 1 }
        .map { meta, vcfs, csis ->
            [ meta, vcfs[0] ]
        }

    SHAPEIT5_LIGATE(ch_multi_chunk)
    ch_versions = ch_versions.mix(SHAPEIT5_LIGATE.out.versions.first())

    FINALIZE_SINGLE_PHASED(ch_single_chunk)

    /*
     * Final phased chromosome panel files = ligated multi-chunk outputs
     *                                       + finalized single-chunk outputs
     */
    ch_final_panel_vcf = SHAPEIT5_LIGATE.out.merged_variants
        .mix(FINALIZE_SINGLE_PHASED.out.phased_vcf)

    VCF_BCFTOOLS_INDEX_2(ch_final_panel_vcf)
    ch_versions = ch_versions.mix(VCF_BCFTOOLS_INDEX_2.out.versions.first())

    ch_vcf_tbi_join = ch_final_panel_vcf
        .join(VCF_BCFTOOLS_INDEX_2.out.csi, failOnMismatch:true, failOnDuplicate:true)

    emit:
    vcf_tbi             = ch_vcf_tbi_join         // channel: [ [id, chr], vcf, csi ]
    versions            = ch_versions             // channel: [ versions.yml ]
}

/*
 * Fallback helper for small chromosomes where GLIMPSE2 sequential writes
 * an empty chunk file even though one chunk is effectively implied.
 *
 * We create one synthetic chunk covering the full region so the chromosome
 * can still enter SHAPEIT5 phasing.
 */
process MAKE_SINGLE_CHUNK_FILE {
    tag "${meta.id ?: 'REF_PANEL'} ${meta.chr}"

    input:
    tuple val(meta), val(region)

    output:
    tuple val(meta), path("REF_PANEL_${meta.chr}_fallback_chunks.txt"), emit: chunkfile

    script:
    """
    echo -e "0\t${meta.chr}\t${region}\t${region}\t0\t0\t0\t0" > REF_PANEL_${meta.chr}_fallback_chunks.txt
    """
}


/*
 * Finalize a single phased chunk into the same final chromosome-level file
 * shape as ligated chromosomes.
 *
 * Example:
 *   REF_PANEL_chr25_538-3066828_chunks.bcf
 * becomes:
 *   REF_PANEL_chr25_phased.vcf.gz
 */
process FINALIZE_SINGLE_PHASED {
    tag "${meta.id} ${meta.chr}"

    input:
    tuple val(meta), path(bcf)

    output:
    tuple val(meta), path("${meta.id}_${meta.chr}_phased.vcf.gz"), emit: phased_vcf

    script:
    """
    bcftools view -Oz -o ${meta.id}_${meta.chr}_phased.vcf.gz ${bcf}
    """
}