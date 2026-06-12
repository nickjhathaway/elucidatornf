process WAIT_FOR_ASSEMBLIES {
    label 'process_single'

    input:
    path expected_names
    path assemblies_root

    output:
    path "assemblies_ready.txt", emit: ready

    script:
    """
    # Completion guard. The per-sample assemblies are consumed from a published
    # directory (`assemblies_root`) that publishDir populates ASYNCHRONOUSLY, so a
    # freshly-extracted sample's results may not be on disk yet when this fires.
    # Wait until every expected sample's result files are present (the same sentinel
    # files the skip-existing check uses) before any downstream step reads the
    # directory, and fail loudly if they never appear.
    deadline=\$(( SECONDS + 3600 ))
    while IFS= read -r name; do
        [ -z "\$name" ] && continue
        until [ -s "${assemblies_root}/\$name/final/allFinal.fasta" ] \\
           && [ -s "${assemblies_root}/\$name/final/basicInfoPerRegion.tab.txt" ] \\
           && [ -s "${assemblies_root}/\$name/partial/allPartial.fasta" ]; do
            if [ \$SECONDS -gt \$deadline ]; then
                echo "ERROR: timed out waiting for assembly results: \$name" >&2
                exit 1
            fi
            sleep 5
        done
    done < ${expected_names}

    echo ok > assemblies_ready.txt
    """
}
