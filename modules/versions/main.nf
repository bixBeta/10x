process DUMP_VERSIONS {

    label 'process_low'

    publishDir "pipeline_info", mode: 'copy', overwrite: true

    input:
        path(versions, stageAs: "versions_?.yml")

    output:
        path "software_versions.yml",     emit: yml

    script:
    """
    cat versions_*.yml | awk '!seen[\$0]++' > software_versions.yml
    """

    stub:
    """
    cat versions_*.yml > software_versions.yml
    """
}
