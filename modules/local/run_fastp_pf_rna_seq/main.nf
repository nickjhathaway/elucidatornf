process FASTP_TRIM {

  tag "${sample_id}_fastp_trimming"
  label 'process_medium_low_memory'
  publishDir { "${pub_dir}" }, mode: 'copy', overwrite: true, pattern: "*.json.gz"

  input:
  tuple val(sample_id), path(reads), val(pub_dir)

  output:
  tuple val(sample_id),
        path("trimmed_${sample_id}_R1.fastq.gz"),
        path("trimmed_${sample_id}_R2.fastq.gz"),
        path("${sample_id}.fastp.json.gz"),
        emit: trimmed

  script:
  def (r1, r2) = reads
  // CLI params arrive as strings, so a plain truthiness test would treat "false" as on
  def end_trim_bases = params.rnaseq_fastp_post_adapter_end_trim.toString().toBoolean()
    ? params.rnaseq_fastp_post_adapter_end_trim_bases.toString().toInteger()
    : 0
  // fastp's --trim_front/--trim_tail act on the raw read, before adapter trimming. The
  // front trim is still right in the main pass (the 5' end is never adapter), but the
  // 3' trim has to come after the adapter is gone, so the main pass streams interleaved
  // reads into a second fastp that only trims the tail (-A -G -Q -L: nothing else runs,
  // no reads are dropped). The main pass's json stays the report, as its counts are final.
  def outputs = end_trim_bases > 0
    ? """--trim_front1 ${end_trim_bases} --trim_front2 ${end_trim_bases} \\
    --stdout \\"""
    : """-o trimmed_${sample_id}_R1.fastq.gz \\
    -O trimmed_${sample_id}_R2.fastq.gz \\"""
  def tail_trim_pass = end_trim_bases > 0
    ? """| fastp \\
    --stdin --interleaved_in \\
    -o trimmed_${sample_id}_R1.fastq.gz \\
    -O trimmed_${sample_id}_R2.fastq.gz \\
    --trim_tail1 ${end_trim_bases} --trim_tail2 ${end_trim_bases} \\
    -A -G -Q -L \\
    --dont_eval_duplication \\
    --thread ${task.cpus} \\
    --json tail_trim.fastp.json \\
    --html tail_trim.fastp.html"""
    : ''
  """

  fastp \\
    -i "${r1}" -I "${r2}" \\
    ${outputs}
    --detect_adapter_for_pe \\
    --cut_mean_quality 20 \\
    --qualified_quality_phred 15 \\
    --unqualified_percent_limit 40 \\
    --n_base_limit 0 \\
    --length_required 50 \\
    --dont_eval_duplication \\
    --trim_poly_g --trim_poly_x \\
    --low_complexity_filter --complexity_threshold 15 \\
    --thread ${task.cpus} \\
    --json "${sample_id}.fastp.json" \\
    --html "${sample_id}.fastp.html" \\
  ${tail_trim_pass}

  pigz -p ${task.cpus} "${sample_id}.fastp.json"
  """

  stub:
  """
  echo | gzip > trimmed_${sample_id}_R1.fastq.gz
  echo | gzip > trimmed_${sample_id}_R2.fastq.gz
  echo '{}' | gzip > ${sample_id}.fastp.json.gz
  """
}
