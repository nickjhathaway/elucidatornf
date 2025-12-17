#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { FASTP_TRIM } from '../../../modules/local/run_fastp_pf'
include { MAP_BWA as HOST1_FILT_MAP_BWA} from '../../../modules/local/map_bwa'
include { MAP_BWA as HOST2_FILT_MAP_BWA} from '../../../modules/local/map_bwa'
include { MAP_BWA as FINAL_MAP_BWA } from '../../../modules/local/map_bwa'
include { MAP_MINIMAP2 as HOST1_FILT_MAP_MINIMAP2 } from '../../../modules/local/map_minimap2'
include { MAP_MINIMAP2 as HOST2_FILT_MAP_MINIMAP2} from '../../../modules/local/map_minimap2'

include { BAM_FILTER_BY_CHROMS as HOST1_MINIMAP2_BAM_FILTER_BY_CHROMS } from '../../../modules/local/bam_filter_by_chroms'
include { BAM_FILTER_BY_CHROMS as HOST1_BWA_BAM_FILTER_BY_CHROMS} from '../../../modules/local/bam_filter_by_chroms'
include { BAM_FILTER_BY_CHROMS as HOST2_MINIMAP2_BAM_FILTER_BY_CHROMS} from '../../../modules/local/bam_filter_by_chroms'
include { BAM_FILTER_BY_CHROMS as HOST2_BWA_BAM_FILTER_BY_CHROMS} from '../../../modules/local/bam_filter_by_chroms'

include { COMBINE_KEPT_FILTERED_FASTQS } from '../../../modules/local/combine_kept_filtered_fastqs'



workflow RUN_WGS_PREPROCESS_PF{
    take:
    input_fastq_dir //a directory with all fastqs to be analyzed
    output_dir // the output directory
    wgs_host1_filter_genome_fasta_fnp
    wgs_host1_filter_contigs_fnp
    wgs_host2_filter_genome_fasta_fnp
    wgs_host2_filter_contigs_fnp
    wgs_final_genome_fasta_fnp
    main:

    /********************************
    * Helper: build rename map
    ********************************/

    def rename_map = SampleRename.loadRenameKey(params.wgs_rename_tsv)


    /********************************
    * Input channel: paired fastqs
    * supports *_1.fastq.gz/_2.fastq.gz and *_R1.fastq.gz/_R2.fastq.gz
    ********************************/

    channel
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
        .filter { _sample_id, _final_id, _reads, skip -> !skip }
        .map { _sample_id, final_id, reads, _skip ->
            tuple(final_id, reads)
        }
        .set { READS_CH }

    def host1_mmi = wgs_host1_filter_genome_fasta_fnp
        .toString()
        .replaceFirst(/\.(fa|fasta|fna)(\.gz)?$/, '.mmi')
    def host1_mmi_path = file(host1_mmi)

    def host2_mmi = wgs_host2_filter_genome_fasta_fnp
        .toString()
        .replaceFirst(/\.(fa|fasta|fna)(\.gz)?$/, '.mmi')
    def host2_mmi_path = file(host2_mmi)


    trimmed_ch = FASTP_TRIM(READS_CH)
    host1_minimap2_in = trimmed_ch.map { sid, r1_trim, r2_trim, _json ->
        tuple(sid, host1_mmi_path, r1_trim, r2_trim)
    }

    host1_minimap2_out = HOST1_FILT_MAP_MINIMAP2(host1_minimap2_in)
    // map minimap2 to host 1 filter
    host1_minimap2_filt_in = host1_minimap2_out.map{sid, bam_fnp, bam_bai_fnp ->
        tuple(sid, bam_fnp, bam_bai_fnp, wgs_host1_filter_contigs_fnp, true)
    }
    // filter by host 1 filter on minimap2
    host1_minimap2_filter_out = HOST1_MINIMAP2_BAM_FILTER_BY_CHROMS(host1_minimap2_filt_in)

    // map bwa to host 1 filter
    host1_bwa_map_in = host1_minimap2_filter_out.map{sid, _kept_r1, _kept_r2, unmapped_r1, unmapped_r2, _chrom_tab, _total_tab ->
        tuple(sid, wgs_host1_filter_genome_fasta_fnp, file("${wgs_host1_filter_genome_fasta_fnp}.*"), unmapped_r1, unmapped_r2)
    }
    host1_bwa_filter_out = HOST1_FILT_MAP_BWA(host1_bwa_map_in)
    // filter by host 1 filter on bwa
    host1_bwa_filt_in = host1_bwa_filter_out.map{sid, bam_fnp, bam_bai_fnp ->
        tuple(sid, bam_fnp, bam_bai_fnp, wgs_host1_filter_contigs_fnp, true)
    }
    //host1_bwa_filter_out =
    HOST1_BWA_BAM_FILTER_BY_CHROMS(host1_bwa_filt_in)

    // map minimap2 to host 2 filter

    // filter by host 2 filter on minimap2

    // map bwa to host 2 filter

    // filter by host 2 filter on bwa

    //combine all kept files

    // map to final bam


}
