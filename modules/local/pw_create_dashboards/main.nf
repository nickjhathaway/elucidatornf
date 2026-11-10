process PW_CREATE_DASHBOARDS {
    label 'process_single'


    input:
    val outdir // the output directory as a string, as a path its modified time ends up in the hash and publishing changes that
    path basic_reports_qmd
    val render_report
    path pop_clustering_res_targets_with_results
    path bed_withGeneInfo_tsv
    path variant_reports_path

    output:
    path "*.html", optional: true

    script:

    if (render_report){
        """
        cp ${basic_reports_qmd} "${outdir}/"
        cd "${outdir}"
        quarto render ${basic_reports_qmd} --no-cache
        """
    } else {
        """
        cp ${basic_reports_qmd} "${outdir}/"
        """
    }
}
