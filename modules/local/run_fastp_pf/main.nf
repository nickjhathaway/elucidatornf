
process FASTP_TRIM {

    tag "${sample_id}"
    label 'process_medium_low_memory'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id),
          path ("trimmed_${sample_id}_R1.fastq.gz"),
          path("trimmed_${sample_id}_R2.fastq.gz"),
          path("${sample_id}.fastp.json")

    script:
    def (r1, r2) = reads
    """

    fastp \\
      -i "${r1}" \\
      -I "${r2}" \\
      -o trimmed_${sample_id}_R1.fastq.gz \\
      -O trimmed_${sample_id}_R2.fastq.gz \\
      --detect_adapter_for_pe \\
      --cut_front \\
      --cut_tail \\
      --cut_window_size 2 \\
      --trim_tail1 1 \\
      --trim_tail2 1 \\
      --cut_mean_quality 20 \\
      --qualified_quality_phred 15 \\
      --unqualified_percent_limit 40 \\
      --n_base_limit 0 \\
      --length_required 50 \\
      --dont_eval_duplication \\
      --trim_poly_g \\
      --trim_poly_x \\
      --low_complexity_filter \\
      --complexity_threshold 15 \\
      --thread ${task.cpus} \\
      --json "${sample_id}.fastp.json" \\
      --html "${sample_id}.fastp.html"

    """
}
