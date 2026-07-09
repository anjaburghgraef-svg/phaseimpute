// Module: MIXBLUP_TO_PLINK
// Converts Mixblup format files to PLINK binary format
//
// map.mix format: index snp_name alleles(2chars) chromosome position
// manifest format: snp_name Y/N value (Y=include marker)
// gtp.mix format: sample_id followed by 0/1/2 genotypes (no spaces between genotypes)

process MIXBLUP_TO_PLINK {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plink:1.90b6.21--h779adbc_1':
        'quay.io/biocontainers/plink:1.90b6.21--h779adbc_1' }"

    input:
    tuple val(meta), path(gtp), path(ped), path(map_mix), path(manifest)

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
    # Step 1: Build list of markers to include (Y in manifest)
    awk '\$2 == "Y" {print \$1}' ${manifest} > included_snps.txt
    echo "Markers to include: \$(wc -l < included_snps.txt)"

    # Step 2: Create filtered map file and allele lookup
    # map.mix format: index snp_name alleles chr position
    # Output map format: chr snp_name 0 position
    # Also create allele file: snp_name allele1 allele2
    awk '
    BEGIN { OFS="\\t" }
    NR==FNR { include[\$1]=1; next }
    \$2 in include {
        # Output to map file
        print \$4, \$2, 0, \$5 > "${prefix}.map"
        # Output alleles (first char = allele1, second char = allele2)
        a1 = substr(\$3, 1, 1)
        a2 = substr(\$3, 2, 1)
        print \$2, a1, a2 > "alleles.txt"
        # Store marker index (1-based from map.mix) for filtering genotypes
        print NR > "marker_indices.txt"
    }
    ' included_snps.txt ${map_mix}

    echo "Created map file with \$(wc -l < ${prefix}.map) markers"

    # Step 3: Build sample info lookup from ped.mix
    # ped.mix format: sample_id sire dam sex ...
    awk '{ print \$1, \$2, \$3, \$4, \$5, \$6 }' ${ped} > sample_info.txt

    # Step 4: Convert gtp.mix to PLINK ped format using awk
    # Read alleles into array, then process genotypes
    awk '
    BEGIN { OFS="\\t" }

    # First file: read marker indices to include (1-based line numbers from map.mix)
    FILENAME == "marker_indices.txt" {
        marker_idx[NR] = \$1 - 1  # Convert to 0-based index into genotype string
        n_markers = NR
        next
    }

    # Second file: read alleles
    FILENAME == "alleles.txt" {
        a1[NR] = \$2
        a2[NR] = \$3
        next
    }

    # Third file: read sample info into lookup
    FILENAME == "sample_info.txt" {
        sample_sire[\$1] = \$2
        sample_dam[\$1] = \$3
        sample_sex[\$1] = \$4
        sample_pheno[\$1] = (\$5 != "") ? \$5 : "0"
        next
    }

    # Fourth file: process genotypes
    {
        # Split line: first field is sample_id, rest is genotypes
        sample_id = \$1
        gsub(/^[^ ]+ /, "", \$0)  # Remove sample_id and space
        gsub(/ /, "", \$0)        # Remove any remaining spaces
        geno_str = \$0

        # Get sample info (use defaults if not found)
        sire = (sample_id in sample_sire) ? sample_sire[sample_id] : "0"
        dam = (sample_id in sample_dam) ? sample_dam[sample_id] : "0"
        sex = (sample_id in sample_sex) ? sample_sex[sample_id] : "0"
        pheno = (sample_id in sample_pheno) ? sample_pheno[sample_id] : "0"

        # Build output line: FID IID father mother sex phenotype genotypes...
        printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s", sample_id, sample_id, sire, dam, sex, pheno

        # Convert each included marker
        for (i = 1; i <= n_markers; i++) {
            idx = marker_idx[i]
            g = substr(geno_str, idx + 1, 1)  # +1 because substr is 1-based

            if (g == "0") {
                printf "\\t%s %s", a1[i], a1[i]
            } else if (g == "1") {
                printf "\\t%s %s", a1[i], a2[i]
            } else if (g == "2") {
                printf "\\t%s %s", a2[i], a2[i]
            } else {
                printf "\\t0 0"
            }
        }
        printf "\\n"
    }
    ' marker_indices.txt alleles.txt sample_info.txt ${gtp} > ${prefix}.ped

    echo "Created ped file with \$(wc -l < ${prefix}.ped) samples"

    # Step 5: Convert to binary PLINK format
    plink \\
        --ped ${prefix}.ped \\
        --map ${prefix}.map \\
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
