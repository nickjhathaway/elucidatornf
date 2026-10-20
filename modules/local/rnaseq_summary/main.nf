process RNASEQ_SAMPLE_SUMMARY {

    tag "${sample_id}_rnaseq_summary"
    label 'process_single'
    publishDir { "${pub_dir}" }, mode: 'copy', overwrite: true, pattern: "*_rnaseq_summary.tab.txt"

    input:
    // the initial and filtered salmon files share names, so stage them apart
    tuple val(sample_id),
          path(fastp_json),
          path(host_total_read_counts),
          path(kmer_filter_counts),
          path(kmer_sort_counts),
          path(initial_salmon_meta, stageAs: 'initial_salmon/*'),
          path(initial_salmon_lib_format, stageAs: 'initial_salmon/*'),
          path(filtered_salmon_meta, stageAs: 'filtered_salmon/*'),
          val(pub_dir)

    output:
    tuple val(sample_id), path("${sample_id}_rnaseq_summary.tab.txt")

    script:
    """
    rnaseq_sample_summary.py \\
        --sample ${sample_id} \\
        --fastp_json ${fastp_json} \\
        --host_counts ${host_total_read_counts} \\
        --kmer_filter_counts ${kmer_filter_counts} \\
        --kmer_sort_counts ${kmer_sort_counts} \\
        --initial_salmon_meta ${initial_salmon_meta} \\
        --initial_salmon_lib_format ${initial_salmon_lib_format} \\
        --filtered_salmon_meta ${filtered_salmon_meta} \\
        --pf_sets ${params.rnaseq_summary_pf_sets} \\
        --sort_exclude_sets ${params.rnaseq_kmer_sort_exclude_sets} \\
        --gene_sets ${params.rnaseq_summary_gene_sets} \\
        --out ${sample_id}_rnaseq_summary.tab.txt
    """

    stub:
    """
    echo ${sample_id} > ${sample_id}_rnaseq_summary.tab.txt
    """
}


process RNASEQ_RUN_SUMMARY {

    tag "rnaseq_run_summary"
    label 'process_single'
    publishDir { "${pub_dir}" }, mode: 'copy', overwrite: true, pattern: "*.tsv.gz"

    input:
    tuple path(sample_summaries), val(pub_dir)

    output:
    path("rnaseq_summary.tsv.gz")

    script:
    // --fill: the filter-pass columns follow the sets of the k-mer table
    """
    elucidator rBind \\
        --contains _rnaseq_summary.tab.txt \\
        --delim tab \\
        --header \\
        --fill \\
        --overWrite \\
        --out rnaseq_summary.tsv.gz
    """

    stub:
    """
    cat *_rnaseq_summary.tab.txt | gzip > rnaseq_summary.tsv.gz
    """
}


process RNASEQ_RESULTS_README {

    tag "rnaseq_results_readme"
    label 'process_single'
    publishDir { "${pub_dir}" }, mode: 'copy', overwrite: true, pattern: "README.md"

    input:
    // the run summary is staged so the README is written last and can count the samples;
    // refs: the reference paths as resolved, since each may come from the directory or a flag
    tuple path(run_summary), path(readme_template), val(refs), val(pub_dir)

    output:
    path("README.md")

    script:
    """
    n_samples=\$(gzip -dc ${run_summary} | tail -n +2 | wc -l | tr -d ' ')
    sed \\
        -e "s|@PIPELINE_VERSION@|${workflow.manifest.version}|" \\
        -e "s|@RUN_DATE@|\$(date +%Y-%m-%d)|" \\
        -e "s|@N_SAMPLES@|\${n_samples}|" \\
        -e "s|@PF_SETS@|${params.rnaseq_summary_pf_sets}|" \\
        -e "s|@GENE_SETS@|${params.rnaseq_summary_gene_sets}|" \\
        -e "s|@PATHWEAVER_SETS@|${params.rnaseq_pathweaver_gene_sets ?: 'none'}|" \\
        -e "s|@HOST1_GENOME@|${refs.host1_genome}|" \\
        -e "s|@HOST1_CONTIGS@|${refs.host1_contigs}|" \\
        -e "s|@HOST2_GENOME@|${refs.host2_genome}|" \\
        -e "s|@HOST2_CONTIGS@|${refs.host2_contigs}|" \\
        -e "s|@SALMON_INDEX@|${refs.salmon_index}|" \\
        -e "s|@KMER_FILTER_TABLE@|${refs.kmer_filter_table}|" \\
        -e "s|@KMER_SORT_TABLE@|${refs.kmer_sort_table}|" \\
        -e "s|@SINGLE_CELL@|${refs.single_cell}|" \\
        -e "s|@NON_UNIQUE_KMER_TABLES@|${refs.non_unique_kmer_tables}|" \\
        -e "s|@REF_FILES_DIR@|${refs.ref_files_dir}|" \\
        ${readme_template} > README.md

    if grep -q '@[A-Z_]*@' README.md; then
        echo "unfilled placeholders in README.md:" >&2
        grep -o '@[A-Z_]*@' README.md | sort -u >&2
        exit 1
    fi
    """

    stub:
    """
    cp ${readme_template} README.md
    """
}
