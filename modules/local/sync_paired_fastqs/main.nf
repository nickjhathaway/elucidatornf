process SYNC_PAIRED_FASTQS {

    tag "${sample_id} sync pairs"

    label 'process_low'

    publishDir { "${pub_dir}" }, mode: 'copy', overwrite: true, pattern: "*_R{1,2}.fastq.gz"
    publishDir { "${unpaired_dir}" }, mode: 'copy', overwrite: true, pattern: "*_unpaired.fastq.gz"

    input:
    tuple val(sample_id), val(pub_dir), val(unpaired_dir), path(kept_r1), path(kept_r2)

    output:
    tuple val(sample_id),
          path("${sample_id}_R1.fastq.gz"),
          path("${sample_id}_R2.fastq.gz"),                     emit: paired
    tuple val(sample_id), path("${sample_id}_pairing.tab.txt"), emit: stats
    path("${sample_id}_unpaired.fastq.gz"), optional: true,     emit: unpaired

    script:
    // seqkit pair needs both mates in the same order, which is how BamFilterByChroms
    // writes them, so this stays a cheap streaming pass rather than a full name index
    def save_unpaired = params.wgs_depletion_save_unpaired ? '--save-unpaired' : ''
    """
    count_reads () {
        echo \$(( \$(gzip -dc "\$1" | wc -l) / 4 ))
    }

    in_r1=\$(count_reads "${kept_r1}")
    in_r2=\$(count_reads "${kept_r2}")

    mkdir -p paired
    seqkit pair \\
        --read1 "${kept_r1}" \\
        --read2 "${kept_r2}" \\
        --out-dir paired \\
        ${save_unpaired} \\
        --force \\
        2> ${sample_id}.seqkit_pair.log.txt

    # seqkit keeps the input basenames in --out-dir, so pull them back out by mate
    mv paired/\$(basename "${kept_r1}") ${sample_id}_R1.fastq.gz
    mv paired/\$(basename "${kept_r2}") ${sample_id}_R2.fastq.gz

    # discarded singletons land in paired/*.unpaired.fastq.gz, one file per mate
    if [ -n "${save_unpaired}" ]; then
        unpaired_files=\$(find paired -name '*.unpaired.*' -type f | sort)
        if [ -n "\$unpaired_files" ]; then
            cat \$unpaired_files > ${sample_id}_unpaired.fastq.gz
        else
            gzip -c < /dev/null > ${sample_id}_unpaired.fastq.gz
        fi
    fi

    out_r1=\$(count_reads "${sample_id}_R1.fastq.gz")
    out_r2=\$(count_reads "${sample_id}_R2.fastq.gz")

    if [ "\$out_r1" -ne "\$out_r2" ]; then
        echo "ERROR: ${sample_id} mate counts differ after pairing (\$out_r1 vs \$out_r2)" >&2
        exit 1
    fi

    dropped=\$(( (in_r1 - out_r1) + (in_r2 - out_r2) ))

    {
        printf "sample_id\\tkept_r1_reads\\tkept_r2_reads\\tfinal_pairs\\tfinal_reads\\tsingletons_discarded\\n"
        printf "%s\\t%d\\t%d\\t%d\\t%d\\t%d\\n" \\
            "${sample_id}" "\$in_r1" "\$in_r2" "\$out_r1" "\$(( out_r1 + out_r2 ))" "\$dropped"
    } > ${sample_id}_pairing.tab.txt
    """
}
