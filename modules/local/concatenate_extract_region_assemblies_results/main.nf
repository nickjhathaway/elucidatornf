process CONCATENATE_EXTRACT_REGION_ASSEMBLIES_RESULTS {
    label 'process_low'

    publishDir { "${pub_results_dir}" }, mode: 'copy', overwrite: true

    input:
    path res_folders
    val pub_results_dir

    output:
    path "allBasicInfo.tsv.gz", emit: all_basic_info
    path "allPartial.fasta.gz", emit: all_partial
    path "allFinal.fasta.gz", emit: all_final

    script:
    """
    rm -f dirs.txt allBasicInfoFiles.txt allFinalFastaFnps.txt allPartialFastaFnps.txt

    # List the staged assembly dirs WITHOUT expanding tens of thousands of paths
    # onto a command line. printf is a bash builtin, so the glob is expanded
    # in-shell and is not subject to ARG_MAX (unlike `ls */` or `for x in \${...}`).
    shopt -s nullglob
    printf '%s\\n' */ | sed 's#/\$##' > dirs.txt

    # Build the per-type file manifests (one path per line).
    sed 's#\$#/final/basicInfoPerRegion.tab.txt#' dirs.txt > allBasicInfoFiles.txt
    sed 's#\$#/final/allFinal.fasta#'             dirs.txt > allFinalFastaFnps.txt
    sed 's#\$#/partial/allPartial.fasta#'         dirs.txt > allPartialFastaFnps.txt

    elucidator rBind --files allBasicInfoFiles.txt --delim tab --header --overWrite --out allBasicInfo.tsv.gz

    # xargs batches the file list under ARG_MAX; `cat \$(cat ...)` would overflow.
    # NUL-delimited so it is robust to any path and to GNU/BSD xargs.
    tr '\\n' '\\0' < allFinalFastaFnps.txt   | xargs -0 cat | pigz > allFinal.fasta.gz
    tr '\\n' '\\0' < allPartialFastaFnps.txt | xargs -0 cat | pigz > allPartial.fasta.gz
    """
}
