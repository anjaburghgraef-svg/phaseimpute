process MAF_SPECTRUM_TOOL {

    tag "maf_spectrum"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' &&
        !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'biocontainers/gawk:5.3.0' }"

    input:
    path err_grp_files

    output:
    path "maf_spectrum_tool_mqc.json", emit: mqc

    script:
"""
set -euo pipefail
command -v gzip >/dev/null 2>&1 || { echo "ERROR: gzip not found" >&2; exit 1; }

cat > maf_spectrum_tool_mqc.json <<'HDR'
{
  "id": "maf_spectrum_tool",
  "section_name": "Imputed MAF spectrum by tool",
  "description": "Mean imputed dosage r² stratified by minor allele frequency (MAF) bins.",
  "plot_type": "linegraph",
  "pconfig": {
    "id": "maf_spectrum_tool_plot",
    "title": "Imputed MAF spectrum by tool",
    "xlab": "Minor allele frequency",
    "ylab": "Mean dosage r²"
  },
  "data": {
HDR

gawk '
BEGIN {
  FS = "[[:space:]]+"
  first_tool = 1

  # collect filenames from ARGV; prevent awk from reading them directly
  n = 0
  for (i = 1; i < ARGC; i++) { files[++n] = ARGV[i]; ARGV[i] = "" }
}

function tool_from_file(fn, m) {
  if (match(fn, /_T([^\\.]+)/, m)) return m[1]
  return "unknown"
}

END {
  # accumulate sums
  for (i = 1; i <= n; i++) {
    fn = files[i]
    if (index(fn, "concordance.rsquare.grp.txt.gz") == 0) continue

    tool = tool_from_file(fn)
    cmd = "gzip -dc " fn

    while ((cmd | getline line) > 0) {
      if (line == "" || line ~ /^#/) continue
      split(line, f, FS)

      b   = f[1]      # bin index
      maf = f[3]      # mean AF
      r2  = f[5]      # imputed dosage r2

      if (r2 == "" || r2 == "NA") continue

      maf_sum[b] += maf; maf_n[b]++
      r2_sum[tool,b] += r2; r2_n[tool,b]++
      tools[tool] = 1
    }
    close(cmd)
  }

  # mean maf per bin (x)
  for (b = 0; b <= 4; b++) maf_mean[b] = (maf_n[b] ? maf_sum[b]/maf_n[b] : 0)

  # write JSON traces
  for (tool in tools) {
    if (!first_tool) printf(",\\n") >> "maf_spectrum_tool_mqc.json"
    first_tool = 0

    printf("    \\"%s\\": [", tool) >> "maf_spectrum_tool_mqc.json"
    first_pt = 1
    for (b = 0; b <= 4; b++) {
      if (r2_n[tool,b]) {
        mean_r2 = r2_sum[tool,b] / r2_n[tool,b]
        if (!first_pt) printf(", ") >> "maf_spectrum_tool_mqc.json"
        first_pt = 0
        printf("[%0.7f, %0.6f]", maf_mean[b], mean_r2) >> "maf_spectrum_tool_mqc.json"
      }
    }
    printf("]") >> "maf_spectrum_tool_mqc.json"
  }

  printf("\\n  }\\n}\\n") >> "maf_spectrum_tool_mqc.json"
}
' ${err_grp_files}

"""
}