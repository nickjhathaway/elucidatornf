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
    rm -f allBasicInfoFiles.txt
    rm -f allFinalFastaFnps.txt
    rm -f allPartialFastaFnps.txt
    for x in ${res_folders}; do
        echo \${x}/final/basicInfoPerRegion.tab.txt >> allBasicInfoFiles.txt
        echo \${x}/final/allFinal.fasta >> allFinalFastaFnps.txt
        echo \${x}/partial/allPartial.fasta >> allPartialFastaFnps.txt
    done

    elucidator rBind --files allBasicInfoFiles.txt --delim tab --header --overWrite --out allBasicInfo.tsv.gz

    cat allFinalFastaFnps.txt
    cat \$(cat allFinalFastaFnps.txt) | pigz > allFinal.fasta.gz
    cat \$(cat allPartialFastaFnps.txt) | pigz > allPartial.fasta.gz

    """
}
