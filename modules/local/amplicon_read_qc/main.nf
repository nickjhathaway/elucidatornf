process AMPLICON_READ_QC {
    // all the samples go in one task so a sequence that's in many samples only gets BLASTed once
    tag "amplicon_read_qc"
    cpus   { ncpus }

    publishDir { "${pub_dir}" }, mode: 'copy', overwrite: true, pattern: "*.tsv.gz"
    publishDir { "${pub_dir}" }, mode: 'copy', overwrite: true, pattern: "parametersUsed.txt"

    input:
    tuple val(sample_names), path(fastq1_fnps, arity: '1..*'), path(fastq2_fnps, arity: '1..*')
    path primers_fnp, name: "primers.tsv"
    path for_seek_deep_info, name: "target_info/forSeekDeep"
    path locations_by_genome, name: "target_info/locationsByGenome"
    path ref_seqs_dir, name: "qc_ref_seqs"
    path genome_dir, name: "target_genomes"
    val primary_genome
    path host_genome_dir, name: "host_genomes" // the directory of the host fasta so its blast database next to it gets used, params.empty_file_fnp for no host
    val host_genome_name // the name of the fasta in host_genome_dir, blank for no host
    path auto_flags_fnp
    val ncpus
    val pub_dir

    output:
    path "readCategoryCounts.tsv.gz", emit: read_category_counts
    path "primerDimerCounts.tsv.gz", emit: primer_dimer_counts
    path "targetCounts.tsv.gz", emit: target_counts
    path "otherPrimerPairCounts.tsv.gz", emit: other_primer_pair_counts
    path "shortProducts.tsv.gz", emit: short_products
    path "kmerMismatchProducts.tsv.gz", emit: kmer_mismatch_products
    path "mismatchedPrimerProducts.tsv.gz", emit: mismatched_primer_products
    path "placementCounts.tsv.gz", emit: placement_counts
    path "parametersUsed.txt", emit: parameters_used

    script:
    def extra_args = task.ext.args ? task.ext.args : ''
    def host_arg = host_genome_name ? "--hostGenome ${host_genome_dir}/${host_genome_name}" : ''
    def sample_table_args = [sample_names, fastq1_fnps, fastq2_fnps].transpose().collect { sample_name, fastq1, fastq2 ->
        "'${sample_name}' '${fastq1}' '${fastq2}'"
    }.join(" ")

    """
    printf 'sample\\tr1\\tr2\\n' > samples.tsv
    printf '%s\\t%s\\t%s\\n' ${sample_table_args} >> samples.tsv
    SeekDeep ampliconReadQC \
            --sampleTable samples.tsv \
            --primers ${primers_fnp} \
            --targetInfoDir target_info \
            --refSeqsDir ${ref_seqs_dir} \
            --genomeDir ${genome_dir} \
            --primaryGenome ${primary_genome} \
            ${host_arg} \
            --numThreads ${ncpus} \
            --dout read_qc \
            ${extra_args} \$(cat ${auto_flags_fnp})
    ln -s read_qc/*.tsv.gz .
    ln -s read_qc/parametersUsed.txt parametersUsed.txt
    """
}
