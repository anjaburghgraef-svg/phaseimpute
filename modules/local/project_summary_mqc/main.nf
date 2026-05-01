process PROJECT_SUMMARY_MQC {
    tag "project_summary"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'biocontainers/gawk:5.3.0' }"

    input:
    path regions_txt,           stageAs: 'regions.txt'
    path imputed_samples_txt,   stageAs: 'samples_imputed.txt'
    path validated_samples_txt, stageAs: 'samples_validated.txt'
    path truth_samples_txt,     stageAs: 'samples_truth.txt'

    output:
    path "project_summary_mqc.txt", emit: mqc

    script:
    '''
    # Join regions into one readable string
    regions_str=$(awk '{for (i=1; i<=NF; i++) {printf (c++ ? " | " : "") $i}} END { print "" }' regions.txt)    n_imputed=$(wc -l < samples_imputed.txt | tr -d ' ')
    n_validated=$(wc -l < samples_validated.txt | tr -d ' ')
    n_truth=$(wc -l < samples_truth.txt | tr -d ' ')

    cat > project_summary_mqc.txt <<'EOF'
# id: "project_summary"
# section_name: "Project summary"
# description: "High-level summary metrics for this run."
# plot_type: "table"
# pconfig:
#   id: "project_summary_table"
#   title: "Project summary"
#   col1_header: "Metric"
Metric\tValue
EOF

    printf "Regions processed\t%s\n" "${regions_str}" >> project_summary_mqc.txt
    printf "Number of samples imputed\t%s\n" "${n_imputed}" >> project_summary_mqc.txt
    printf "Number of samples validated\t%s\n" "${n_validated}" >> project_summary_mqc.txt
    printf "Number of truth samples used\t%s\n" "${n_truth}" >> project_summary_mqc.txt
    '''
}