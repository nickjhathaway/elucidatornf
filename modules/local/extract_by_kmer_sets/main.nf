// Read sorting by unique k-mer sets with elucidator extractByCountingUniqKmersFromSets.
// Used by the RNA-seq pipeline as two passes over the host-depleted reads:
//   KMER_SET_FILTER  removes reads matching unwanted sets (e.g. mycoplasma, Pf rRNA)
//   KMER_SET_SORT    sorts the remaining reads into the genes/gene families of the table
// Each pass writes three count tables (extractionCounts, multihitCounts,
// inversePairsCounts), already carrying a sample column.

process KMER_SET_FILTER {

    tag "${sample_id} kmer set filter"
    label 'process_medium_low_memory'
    publishDir { "${pub_dir}/counts" }, mode: 'copy', overwrite: true, pattern: "*.tab.txt"

    input:
    tuple val(sample_id), path(r1), path(r2), path(kmer_table), val(pub_dir)

    output:
    tuple val(sample_id),
          path("${sample_id}_kmer_filter_kept_R1.fastq.gz"),
          path("${sample_id}_kmer_filter_kept_R2.fastq.gz"),
          emit: kept
    tuple val(sample_id),
          path("${sample_id}_kmer_filter_extractionCounts.tab.txt"),
          path("${sample_id}_kmer_filter_multihitCounts.tab.txt"),
          path("${sample_id}_kmer_filter_inversePairsCounts.tab.txt"),
          emit: counts

    script:
    // the kept sets are passed as --excludeSetNames: an excluded set is not competed
    // for, and with --writeOutExclude its reads are written out so they can be kept
    def kept_sets  = params.rnaseq_kmer_filter_kept_sets.toString().tokenize(',')*.trim()
    def extra_args = params.rnaseq_kmer_filter_extra_args ?: ''
    """
    elucidator extractByCountingUniqKmersFromSets \\
        --fastq1gz ${r1} \\
        --fastq2gz ${r2} \\
        --sampleName ${sample_id} \\
        --kmerTable ${kmer_table} \\
        --dout kmer_filter \\
        --overWriteDir \\
        --numThreads ${task.cpus} \\
        --maxIterations 0 \\
        --doNotDoFinalExtract \\
        --readKmerFracCutOff 0.01 \\
        --hardCountOff 15 \\
        --noCheckingQFrac \\
        --smallLenCutOff 0 \\
        --excludeSetNames ${kept_sets.join(',')} \\
        --writeOutExclude \\
        --filterMultiHitReads \\
        --writeMultiHitSeparate \\
        ${extra_args}

    # kept reads = the kept sets plus the reads that matched no set. A set with no
    # reads gets no file, so build the outputs up from empty gzip members
    gzip -c </dev/null > ${sample_id}_kmer_filter_kept_R1.fastq.gz
    gzip -c </dev/null > ${sample_id}_kmer_filter_kept_R2.fastq.gz
    for set_name in ${kept_sets.join(' ')} undetermined; do
        if [ -f "kmer_filter/\${set_name}_R1.fastq.gz" ]; then
            cat "kmer_filter/\${set_name}_R1.fastq.gz" >> ${sample_id}_kmer_filter_kept_R1.fastq.gz
            cat "kmer_filter/\${set_name}_R2.fastq.gz" >> ${sample_id}_kmer_filter_kept_R2.fastq.gz
        fi
    done

    for tab in extractionCounts multihitCounts inversePairsCounts; do
        cp "kmer_filter/\${tab}.tab.txt" "${sample_id}_kmer_filter_\${tab}.tab.txt"
    done
    """

    stub:
    """
    echo | gzip > ${sample_id}_kmer_filter_kept_R1.fastq.gz
    echo | gzip > ${sample_id}_kmer_filter_kept_R2.fastq.gz
    for tab in extractionCounts multihitCounts inversePairsCounts; do touch ${sample_id}_kmer_filter_\${tab}.tab.txt; done
    """
}


process KMER_SET_SORT {

    tag "${sample_id} kmer set sort"
    label 'process_medium_low_memory'
    publishDir { "${pub_dir}/sorted_fastqs" }, mode: 'copy', overwrite: true, pattern: "*/*.fastq.gz"
    publishDir { "${pub_dir}/counts" }, mode: 'copy', overwrite: true, pattern: "*.tab.txt"

    input:
    tuple val(sample_id), path(r1), path(r2), path(kmer_table), val(pub_dir)

    output:
    tuple val(sample_id), path("${sample_id}/*.fastq.gz"), emit: fastqs, optional: true
    tuple val(sample_id),
          path("${sample_id}_kmer_sort_extractionCounts.tab.txt"),
          path("${sample_id}_kmer_sort_multihitCounts.tab.txt"),
          path("${sample_id}_kmer_sort_inversePairsCounts.tab.txt"),
          emit: counts

    script:
    // CLI params arrive as strings, so a plain truthiness test would treat "false" as on
    def freq_correction = params.rnaseq_kmer_sort_frequency_correction.toString().toBoolean()
    def freq_args = freq_correction
        ? "--kmerFrequencyCorrection --kmerFrequencyCorrectionMaxCandidateDepth ${params.rnaseq_kmer_sort_frequency_max_candidate_depth}"
        : ''
    def extra_args = params.rnaseq_kmer_sort_extra_args ?: ''
    """
    elucidator extractByCountingUniqKmersFromSets \\
        --fastq1gz ${r1} \\
        --fastq2gz ${r2} \\
        --sampleName ${sample_id} \\
        --kmerTable ${kmer_table} \\
        --dout kmer_sort \\
        --overWriteDir \\
        --numThreads ${task.cpus} \\
        --smallLenCutOff 75 \\
        --maxIterations 10 \\
        --excludeSetNames ${params.rnaseq_kmer_sort_exclude_sets} \\
        --filterMultiHitReads \\
        ${freq_args} \\
        ${extra_args}

    # only the sorted gene fastqs go on; the excluded sets are never written here
    mkdir ${sample_id}
    for fq in kmer_sort/finalExtraction/*.fastq.gz; do
        [ -e "\${fq}" ] || continue
        case "\$(basename "\${fq}")" in
            undetermined*) continue ;;
        esac
        mv "\${fq}" ${sample_id}/
    done

    for tab in extractionCounts multihitCounts inversePairsCounts; do
        cp "kmer_sort/\${tab}.tab.txt" "${sample_id}_kmer_sort_\${tab}.tab.txt"
    done
    """

    stub:
    """
    mkdir ${sample_id}
    echo | gzip > ${sample_id}/vars_R1.fastq.gz
    for tab in extractionCounts multihitCounts inversePairsCounts; do touch ${sample_id}_kmer_sort_\${tab}.tab.txt; done
    """
}


process COMBINE_KMER_EXTRACTION_COUNTS {

    tag "combine ${prefix} counts"
    label 'process_single'
    publishDir { "${pub_dir}" }, mode: 'copy', overwrite: true, pattern: "*.tsv.gz"

    input:
    tuple val(prefix), path(count_tabs), val(pub_dir)

    output:
    tuple path("${prefix}_extractionCounts.tsv.gz"),
          path("${prefix}_multihitCounts.tsv.gz"),
          path("${prefix}_inversePairsCounts.tsv.gz")

    script:
    """
    for tab in extractionCounts multihitCounts inversePairsCounts; do
        elucidator rBind --contains _\${tab}.tab.txt --delim tab --header --overWrite --out ${prefix}_\${tab}.tsv.gz
    done
    """

    stub:
    // the staged count tables, so a -stub run shows which samples were combined
    """
    for tab in extractionCounts multihitCounts inversePairsCounts; do
        ls *_\${tab}.tab.txt | gzip > ${prefix}_\${tab}.tsv.gz
    done
    """
}
