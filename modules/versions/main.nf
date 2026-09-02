process DUMP_VERSIONS {

    label 'process_low'

    publishDir "pipeline_info", mode: 'copy', overwrite: true

    input:
        path(versions, stageAs: "versions_?.yml")

    output:
        path "software_versions.yml"     , emit: yml
        path "software_versions_mqc.yml" , emit: mqc_yml

    script:
    """
    cat versions_*.yml | awk '!seen[\$0]++' > software_versions.yml

    # MultiQC custom content: anything named *_mqc.yml is picked up as a
    # section. Built with awk rather than python so this process needs nothing
    # beyond coreutils - it runs on the host, not in an image.
    {
      echo "id: pipeline_software_versions"
      echo "section_name: Software Versions"
      echo "description: Versions of the software this run used, as reported by the tools themselves."
      echo "section_href: https://github.com/bixBeta/10x"
      echo "plot_type: html"
      echo "data: |"
      echo "  <table class=\\"table table-condensed\\">"
      echo "  <thead><tr><th>Software</th><th>Version</th></tr></thead>"
      echo "  <tbody>"
      awk '
        /^[[:space:]]+[A-Za-z0-9_.+-]+:[[:space:]]*/ {
            line = \$0
            sub(/^[[:space:]]+/, "", line)
            idx  = index(line, ":")
            tool = substr(line, 1, idx - 1)
            ver  = substr(line, idx + 1)
            gsub(/^[[:space:]]+|[[:space:]]+\$/, "", ver)
            if ( ver != "" && !(tool in seen) ) { seen[tool] = ver; order[++n] = tool }
        }
        END {
            for ( i = 1; i <= n; i++ )
                printf "  <tr><td><strong>%s</strong></td><td><samp>%s</samp></td></tr>\\n", order[i], seen[order[i]]
        }
      ' software_versions.yml
      echo "  </tbody></table>"
    } > software_versions_mqc.yml
    """

    stub:
    """
    cat versions_*.yml > software_versions.yml
    printf 'id: pipeline_software_versions\\nsection_name: Software Versions\\nplot_type: html\\ndata: |\\n  <table></table>\\n' > software_versions_mqc.yml
    """
}
