#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { FASTP_TRIM } from '../../../modules/local/run_fastp_pf'



workflow RUN_WGS_PREPROCESS_PF{
    take:
    input_fastq_dir //a directory with all fastqs to be analyzed
    output_dir // the output directory
    main:

    /********************************
    * Helper: build rename map
    ********************************/

    def rename_map = [:]
    if (params.wgs_rename_tsv) {
        new File(params.wgs_rename_tsv).eachLine { line ->
            if (!line.trim() || line.startsWith('#')) return
            def toks = line.split('\t')
            assert toks.size() >= 2 : "wgs_rename_tsv must have at least 2 columns: orig_basename\tnew_name"
            rename_map[toks[0]] = toks[1]
        }
    }

    /********************************
    * Input channel: paired fastqs
    * supports *_1.fastq.gz/_2.fastq.gz and *_R1.fastq.gz/_R2.fastq.gz
    ********************************/

    Channel
        .fromFilePairs("${input_fastq_dir}/*_{1,2}.fastq.gz", flat: true)
        .mix(
            channel.fromFilePairs("${input_fastq_dir}/*_R{1,2}.fastq.gz", flat: true)
        )
        .mix(
            channel.fromFilePairs("${input_fastq_dir}/*_R{1,2}.fq.gz", flat: true)
        )
        .mix(
            channel.fromFilePairs("${input_fastq_dir}/*_{1,2}.fq.gz", flat: true)
        )
        .map { sid, r1, r2 ->
            /*
            * sid is the common basename from fromFilePairs, e.g. ST131 or ST131_R
            * clean it a bit so:
            *   ST131_1.fastq.gz    -> ST131
            *   ST131_R1.fastq.gz   -> ST131
            */
            def base = sid
                    .replaceAll(/_R?1$/, '')
                    .replaceAll(/_R?2$/, '')
            def sample_id = base
            def final_id  = rename_map.get(sample_id, sample_id)
            tuple(sample_id, final_id, [r1, r2])
        }
        .unique()   // in case patterns overlap
        .map { sample_id, final_id, reads ->
            // Skip if final BAM exists and skip_existing is true
            def final_bam = file("${output_dir}/${final_id}.sorted.bam")
            def skip = params.wgs_skip_existing && final_bam.exists()
            tuple(sample_id, final_id, reads, skip)
        }
        .filter { sample_id, final_id, reads, skip -> !skip }
        .map { sample_id, final_id, reads, skip ->
            tuple(sample_id, reads)
        }
        .set { READS_CH }

    FASTP_TRIM(READS_CH)
}
