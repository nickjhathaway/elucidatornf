process FASTP_TRIM {

  tag "${sample_id}_fastp_trimming"
  label 'process_medium_low_memory'
  publishDir "${pub_dir}", mode: 'copy', overwrite: true, pattern: "*.json.gz"

  maxRetries 2

  input:
  tuple val(sample_id), path(reads), val(pub_dir)

  output:
  tuple val(sample_id),
        path("trimmed_${sample_id}_R1.fastq.gz"),
        path("trimmed_${sample_id}_R2.fastq.gz"),
        path("${sample_id}.fastp.json.gz"),
        emit: trimmed, optional: true

  tuple val(sample_id),
        path("${sample_id}.fastp.status.txt"),
        emit: status

  script:
  def (r1, r2) = reads
  """
  set +e

  fastp \\
    -i "${r1}" -I "${r2}" \\
    -o trimmed_${sample_id}_R1.fastq.gz \\
    -O trimmed_${sample_id}_R2.fastq.gz \\
    --detect_adapter_for_pe \\
    --cut_front --cut_tail --cut_window_size 2 \\
    --trim_tail1 1 --trim_tail2 1 \\
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
    --html "${sample_id}.fastp.html"

  ec=\$?

  if [ \$ec -eq 0 ]; then
    echo "OK" > ${sample_id}.fastp.status.txt
    pigz -p ${task.cpus} "${sample_id}.fastp.json"
    exit 0
  fi

  rm -f trimmed_${sample_id}_R1.fastq.gz trimmed_${sample_id}_R2.fastq.gz \\
        ${sample_id}.fastp.json ${sample_id}.fastp.html

  if [ ${task.attempt} -lt ${task.maxRetries} ]; then
    echo "RETRY\\t${task.attempt}\\t\$ec" > ${sample_id}.fastp.status.txt
    exit \$ec
  else
    echo "FAIL\\t${task.attempt}\\t\$ec" > ${sample_id}.fastp.status.txt
    exit 0
  fi
  """
}


// process FASTP_TRIM {

//     tag "${sample_id}_fastp_trimming"
//     label 'process_medium_low_memory'

//     publishDir "${pub_dir}", mode: 'copy', overwrite: true, pattern: "*.json.gz"

//     input:
//     tuple val(sample_id), path(reads), val(pub_dir)

//     output:
//     tuple val(sample_id),
//           path ("trimmed_${sample_id}_R1.fastq.gz"),
//           path("trimmed_${sample_id}_R2.fastq.gz"),
//           path("${sample_id}.fastp.json.gz")

//     script:
//     def (r1, r2) = reads
//     """

//     fastp \\
//       -i "${r1}" \\
//       -I "${r2}" \\
//       -o trimmed_${sample_id}_R1.fastq.gz \\
//       -O trimmed_${sample_id}_R2.fastq.gz \\
//       --detect_adapter_for_pe \\
//       --cut_front \\
//       --cut_tail \\
//       --cut_window_size 2 \\
//       --trim_tail1 1 \\
//       --trim_tail2 1 \\
//       --cut_mean_quality 20 \\
//       --qualified_quality_phred 15 \\
//       --unqualified_percent_limit 40 \\
//       --n_base_limit 0 \\
//       --length_required 50 \\
//       --dont_eval_duplication \\
//       --trim_poly_g \\
//       --trim_poly_x \\
//       --low_complexity_filter \\
//       --complexity_threshold 15 \\
//       --thread ${task.cpus} \\
//       --json "${sample_id}.fastp.json" \\
//       --html "${sample_id}.fastp.html"
//     pigz -p ${task.cpus} "${sample_id}.fastp.json"
//     """
// }
//
