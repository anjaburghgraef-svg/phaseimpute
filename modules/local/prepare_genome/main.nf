process PREPARE_GENOME {
    tag "$genome_name"
    label 'process_medium'
    storeDir "${params.genomes_cache ?: 'genomes_cache'}/${genome_name}"

    conda "bioconda::samtools=1.21"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0' :
        'biocontainers/samtools:1.21--h50ea8bc_0' }"

    input:
    tuple val(genome_name), val(fasta_url), val(add_chr)

    output:
    tuple val(genome_name), path("${genome_name}.fa"), path("${genome_name}.fa.fai"), emit: genome
    path "versions.yml", emit: versions

    script:
    def add_chr_cmd = add_chr ? "| sed 's/^>\\([0-9XYMTWZUn][0-9a-zA-Z_]*\\)/>chr\\1/'" : ""
    """
    # Download genome (gzip compressed from Ensembl/NCBI)
    wget -q "${fasta_url}" -O genome.fa.gz

    # Decompress and optionally add chr prefix, then output as uncompressed fasta
    # (samtools faidx works with uncompressed fasta)
    gunzip -c genome.fa.gz ${add_chr_cmd} > ${genome_name}.fa

    # Remove downloaded gzip to save space
    rm genome.fa.gz

    # Index with samtools
    samtools faidx ${genome_name}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version 2>&1 | head -1 | sed 's/samtools //')
    END_VERSIONS
    """
}
