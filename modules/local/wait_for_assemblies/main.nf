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
    # Wait until every expected sample's result files are PRESENT before any
    # downstream step reads the directory, and fail loudly if they never appear.
    # Use -f (exists), NOT -s (non-empty): a sample can legitimately have an empty
    # allPartial/allFinal.fasta, and this must accept exactly what the skip-existing
    # check (pwAssemblyComplete -> .exists()) accepts, or already-done runs hang.
    deadline=\$(( SECONDS + 3600 ))
    while :; do
        missing=0
        example=""
        while IFS= read -r name; do
            [ -z "\$name" ] && continue
            if [ ! -f "${assemblies_root}/\$name/final/allFinal.fasta" ] \\
            || [ ! -f "${assemblies_root}/\$name/final/basicInfoPerRegion.tab.txt" ] \\
            || [ ! -f "${assemblies_root}/\$name/partial/allPartial.fasta" ]; then
                missing=\$(( missing + 1 ))
                [ -z "\$example" ] && example="\$name"
            fi
        done < ${expected_names}

        [ "\$missing" -eq 0 ] && break

        if [ \$SECONDS -gt \$deadline ]; then
            echo "ERROR: timed out with \$missing assemblies still missing (e.g. \$example)" >&2
            exit 1
        fi
        echo "[wait_for_assemblies] still waiting on \$missing assemblies (e.g. \$example)" >&2
        sleep 10
    done

    echo ok > assemblies_ready.txt
    """
}
