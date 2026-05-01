process TOOL_ACCURACY_MQC {
    tag "tool_accuracy"
    label 'process_single'

    // Same container selection style as LISTTOFILE
    // Needs: gawk + gzip (or zcat). If gzip is missing in this image, see note below.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gawk:5.3.0' :
        'biocontainers/gawk:5.3.0' }"

    input:
    path conc_files

    output:
    path "tool_accuracy_mqc.txt", emit: mqc

    when:
    task.ext.when == null || task.ext.when

    script:
'''
# MultiQC custom content header (auto-detected because *_mqc.txt) [2](https://docs.galaxyproject.org/en/latest/admin/special_topics/mulled_containers.html)
cat > tool_accuracy_mqc.txt <<'HDR'
# id: "tool_accuracy"
# section_name: "Imputation accuracy by tool"
# description: "Mean accuracy across truth samples (SNPs), grouped by imputation tool."
# plot_type: "table"
Tool\tN_samples\tValidation_variants\tMean_dosage_r2\tMean_best_gt_r2\tMean_nonref_discordance_pct
HDR

# Only the per-sample concordance summaries (these already contain r2 + discordance)
ls -1 *.concordance.renamed.error.spl.txt.gz > _err_spl_files.txt

gawk -v OFS='\\t' '
  BEGIN { FS=" " }   # These files are space-separated (not tab) once decompressed

  function tool_from_sample(s, a) {
    split(s, a, "\\.")
    return (length(a) >= 2 ? a[2] : "unknown")
  }

  # Read list of gz files
  FNR==NR { files[++n] = $0; next }

  END {
    for (i=1; i<=n; i++) {
      fn = files[i]
      cmd = "gzip -dc " fn

      # Reset column index mapping for each file
      delete idx
      header_ok = 0

      while ((cmd | getline line) > 0) {

        # Skip empty lines
        if (line == "") continue

        # Header line for SNP section looks like: "#GCsS id sample_name #val_gt_RR ..."
        if (line ~ /^#GCsS[[:space:]]/) {
          # Build column-name -> index map
          split(line, h, FS)
          for (j=1; j<=length(h); j++) {
            name = h[j]
            sub(/^#/, "", name)     # remove leading #
            idx[name] = j
          }
          header_ok = 1
          continue
        }

        # Ignore other comment lines
        if (line ~ /^#/) continue

        # Only process SNP rows (GCsS ...)
        split(line, f, FS)
        if (f[1] != "GCsS") continue
        if (!header_ok) continue

        sample = f[idx["sample_name"]]
        tool   = tool_from_sample(sample)

        # Validation variants: val_gt_RR + val_gt_RA + val_gt_AA
        rr = f[idx["val_gt_RR"]]
        ra = f[idx["val_gt_RA"]]
        aa = f[idx["val_gt_AA"]]
        vv = rr + ra + aa

        # Non-ref discordance column has a typo in your file: "non_reference_discordanc_rate_percent"
        nrd = f[idx["non_reference_discordanc_rate_percent"]]

        # r2 columns
        best = f[idx["best_gt_rsquared"]]
        dos  = f[idx["imputed_ds_rsquared"]]

        sum_best[tool] += best;  n_best[tool]++
        sum_dos[tool]  += dos;   n_dos[tool]++
        sum_nrd[tool]  += nrd;   n_nrd[tool]++

        # keep vv (should be consistent per tool); store first seen
        if (!(tool in valvar)) valvar[tool] = vv

        # count unique samples per tool (not rows)
        key = tool "|" sample
        if (!(key in seen)) { seen[key] = 1; n_samples[tool]++ }
      }

      close(cmd)
    }

    # Print final table
    for (t in n_samples) {
      mdos  = (n_dos[t]  ? sum_dos[t]/n_dos[t]   : 0)
      mbest = (n_best[t] ? sum_best[t]/n_best[t] : 0)
      mnrd  = (n_nrd[t]  ? sum_nrd[t]/n_nrd[t]   : 0)

      printf "%s\\t%d\\t%s\\t%.6f\\t%.6f\\t%.6f\\n", \\
        t, n_samples[t], (t in valvar ? valvar[t] : ""), mdos, mbest, mnrd \\
        >> "tool_accuracy_mqc.txt"
    }
  }
' _err_spl_files.txt _err_spl_files.txt
'''
}