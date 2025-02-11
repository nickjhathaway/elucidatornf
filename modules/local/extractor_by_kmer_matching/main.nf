process EXTRACTOR_BY_KMER_MATCHING {
    label 'process_single'


    input:
    tuple path(fastq_fnp), val(sample_name), path(primers_fnp), path(kmers_sets), path(length_cut_offs_per_target), val (min_len_cut_off)
    // path fastq_fnp
    // val sample_name
    // path primers_fnp
    // path kmers_sets
    // path length_cut_offs_per_target
    // val min_len_cut_off

    output:
    path "extraction/*.fastq.gz", emit: fastqs_per_target
    path "extraction/extractionProfile.tab.txt", emit: extraction_profile_per_target
    path "extraction/extractionStats.tab.txt", emit: extraction_stats

    script:
    """
    SeekDeep extractorByKmerMatching \
            --dout extraction \
            --uniqueKmersPerTarget ${kmers_sets} \
            --rename \
            --fastqgz ${fastq_fnp} \
            --id ${primers_fnp} \
            --sampleName ${sample_name} \
            --lenCutOffs ${length_cut_offs_per_target} \
            --minLenCutOff ${min_len_cut_off}
    """
}
