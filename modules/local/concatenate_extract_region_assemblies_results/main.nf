process CONCATENATE_EXTRACT_REGION_ASSEMBLIES_RESULTS {
    label 'process_low'

    publishDir { "${pub_results_dir}" }, mode: 'copy', overwrite: true

    input:
    path assemblies_root
    val pub_results_dir

    output:
    path "allBasicInfo.tsv.gz", emit: all_basic_info
    path "allPartial.fasta.gz", emit: all_partial
    path "allFinal.fasta.gz", emit: all_final

    script:
    """
    rm -f dirs.txt allBasicInfoFiles.txt allFinalFastaFnps.txt allPartialFastaFnps.txt

    # `assemblies_root` is the SINGLE staged parent dir holding every per-sample
    # assembly dir. Staging one dir (vs tens of thousands of child paths) keeps
    # nxf_stage tiny -- a huge nxf_stage makes bash segfault -- and Nextflow
    # bind-mounts this one dir into containers automatically. printf is a bash
    # builtin, so the glob is expanded in-shell (no ARG_MAX issue).
    shopt -s nullglob
    printf '%s\\n' ${assemblies_root}/*/ | sed 's#/\$##' > dirs.txt

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
