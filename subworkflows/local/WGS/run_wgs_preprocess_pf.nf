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
include { COMBINE_KEPT_FILTERED_COUNTS } from '../../../modules/local/combine_kept_filtered_fastqs'

include { WGS_BUILD_FINAL_SAMPLE_SUMMARY } from "../../../modules/local/wgs_build_final_sample_summary"


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

    // channel
    //     .fromFilePairs("${input_fastq_dir}/*_{1,2}.fastq.gz", flat: true)
    //     .mix(
    //         channel.fromFilePairs("${input_fastq_dir}/*_R{1,2}.fastq.gz", flat: true)
    //     )
    //     .mix(
    //         channel.fromFilePairs("${input_fastq_dir}/*_R{1,2}.fq.gz", flat: true)
    //     )
    //     .mix(
    //         channel.fromFilePairs("${input_fastq_dir}/*_{1,2}.fq.gz", flat: true)
    //     )
    //     .map { sid, r1, r2 ->
    //         /*
    //         * sid is the common basename from fromFilePairs, e.g. ST131 or ST131_R
    //         * clean it a bit so:
    //         *   ST131_1.fastq.gz    -> ST131
    //         *   ST131_R1.fastq.gz   -> ST131
    //         */
    //         def base = sid
    //                 .replaceAll(/_R?1$/, '')
    //                 .replaceAll(/_R?2$/, '')
    //         def sample_id = base
    //         def final_id  = rename_map.get(sample_id, sample_id)
    //         tuple(sample_id, final_id, [r1, r2])
    //     }
    //     .unique()   // in case patterns overlap
    def fastp_trim_info_dir = file("${output_dir}/fastp_and_filter_info")
    fastp_trim_info_dir.mkdirs()
    channel
        .fromFilePairs(
            "${input_fastq_dir}/*_{R,}{1,2}*.{fq,fastq}.gz",
            flat: true
        )
        .map { sid, r1, r2 ->
            /*
            * sid examples produced by fromFilePairs:
            *   PAT-020
            *   PAT-025
            *   S3_WGS-batchD_8045353103-ST170_S19
            */
            def sample_id = sid
            def final_id  = rename_map.get(sample_id, sample_id)
            tuple(sample_id, final_id, [r1, r2])
        }
        .unique()
        .map { sample_id, final_id, reads ->
            // Skip if final BAM exists and skip_existing is true
            def final_bam = file("${output_dir}/${final_id}.sorted.bam")
            def skip = params.wgs_skip_existing && final_bam.exists()
            tuple(sample_id, final_id, reads, skip)
        }
        .filter { _sample_id, _final_id, _reads, skip -> !skip }
        .map { _sample_id, final_id, reads, _skip ->
            tuple(final_id, reads, fastp_trim_info_dir)
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

    // map minimap2 to host 1 filter
    host1_minimap2_out = HOST1_FILT_MAP_MINIMAP2(host1_minimap2_in)
    // filter by host 1 filter on minimap2
    host1_minimap2_filt_in = host1_minimap2_out.map{sid, bam_fnp, bam_bai_fnp ->
        tuple(sid, "host1_minimap2_filt", bam_fnp, bam_bai_fnp, wgs_host1_filter_contigs_fnp, true)
    }
    host1_minimap2_filter_out = HOST1_MINIMAP2_BAM_FILTER_BY_CHROMS(host1_minimap2_filt_in)

    // map bwa to host 1 filter
    host1_bwa_map_in = host1_minimap2_filter_out.map{sid, _kept_r1, _kept_r2, unmapped_r1, unmapped_r2, _chrom_tab, _total_tab ->
        tuple(sid, wgs_host1_filter_genome_fasta_fnp, file("${wgs_host1_filter_genome_fasta_fnp}.*"), unmapped_r1, unmapped_r2)
    }
    host1_bwa_filter_out = HOST1_FILT_MAP_BWA(host1_bwa_map_in)
    // filter by host 1 filter on bwa
    host1_bwa_filt_in = host1_bwa_filter_out.map{sid, bam_fnp, bam_bai_fnp ->
        tuple(sid, "host1_bwa_filt", bam_fnp, bam_bai_fnp, wgs_host1_filter_contigs_fnp, true)
    }
    host1_bwa_filter_out = HOST1_BWA_BAM_FILTER_BY_CHROMS(host1_bwa_filt_in)

    // map minimap2 to host 2 filter
    host2_minimap2_in = host1_bwa_filter_out.map { sid, _kept_r1, _kept_r2, unmapped_r1, unmapped_r2, _chrom_tab, _total_tab ->
        tuple(sid, host2_mmi_path, unmapped_r1, unmapped_r2)
    }
    host2_minimap2_out = HOST2_FILT_MAP_MINIMAP2(host2_minimap2_in)
    // filter by host 2 filter on minimap2
    host2_minimap2_filt_in = host2_minimap2_out.map{sid, bam_fnp, bam_bai_fnp ->
        tuple(sid, "host2_minimap2_filt", bam_fnp, bam_bai_fnp, wgs_host2_filter_contigs_fnp, true)
    }
    host2_minimap2_filter_out = HOST2_MINIMAP2_BAM_FILTER_BY_CHROMS(host2_minimap2_filt_in)

    // map bwa to host 2 filter
    host2_bwa_map_in = host2_minimap2_filter_out.map{sid, _kept_r1, _kept_r2, unmapped_r1, unmapped_r2, _chrom_tab, _total_tab ->
        tuple(sid, wgs_host2_filter_genome_fasta_fnp, file("${wgs_host2_filter_genome_fasta_fnp}.*"), unmapped_r1, unmapped_r2)
    }
    host2_bwa_filter_out = HOST2_FILT_MAP_BWA(host2_bwa_map_in)
    // filter by host 2 filter on bwa
    host2_bwa_filt_in = host2_bwa_filter_out.map{sid, bam_fnp, bam_bai_fnp ->
        tuple(sid, "host2_bwa_filt", bam_fnp, bam_bai_fnp, wgs_host2_filter_contigs_fnp, false)
    }
    host2_bwa_filter_out = HOST2_BWA_BAM_FILTER_BY_CHROMS(host2_bwa_filt_in)

    //combine all kept files
    // kept reads from each pass
    kept_host1_minimap2 = host1_minimap2_filter_out.map { sid, kept_r1, kept_r2, _unm1, _unm2, _chrom, _total ->
        tuple(sid, kept_r1, kept_r2)
    }

    kept_host1_bwa = host1_bwa_filter_out.map { sid, kept_r1, kept_r2, _unm1, _unm2, _chrom, _total ->
        tuple(sid, kept_r1, kept_r2)
    }

    kept_host2_minimap2 = host2_minimap2_filter_out.map { sid, kept_r1, kept_r2, _unm1, _unm2, _chrom, _total ->
        tuple(sid, kept_r1, kept_r2)
    }

    kept_host2_bwa = host2_bwa_filter_out.map { sid, kept_r1, kept_r2, _unm1, _unm2, _chrom, _total ->
        tuple(sid, kept_r1, kept_r2)
    }
    kept_host1_minimap2_tagged = kept_host1_minimap2.map { sid, r1, r2 -> tuple(sid, 'h1_mm',  r1, r2) }
    kept_host1_bwa_tagged      = kept_host1_bwa.map      { sid, r1, r2 -> tuple(sid, 'h1_bwa', r1, r2) }
    kept_host2_minimap2_tagged = kept_host2_minimap2.map { sid, r1, r2 -> tuple(sid, 'h2_mm',  r1, r2) }
    kept_host2_bwa_tagged      = kept_host2_bwa.map      { sid, r1, r2 -> tuple(sid, 'h2_bwa', r1, r2) }

    all_kept =
        kept_host1_minimap2_tagged
            .mix(kept_host1_bwa_tagged)
            .mix(kept_host2_minimap2_tagged)
            .mix(kept_host2_bwa_tagged)

    // group by sample_id; each record is [tag, r1, r2]
    combine_in =
        all_kept
            .map { sid, tag, r1, r2 -> tuple(sid, [tag, r1, r2]) }
            .groupTuple(size: 4)
            .map { sid, recs ->
                def byTag = recs.collectEntries { r -> [(r[0]): r] }

                tuple(
                    sid,
                    byTag['h1_mm'][1],  byTag['h1_mm'][2],
                    byTag['h1_bwa'][1], byTag['h1_bwa'][2],
                    byTag['h2_mm'][1],  byTag['h2_mm'][2],
                    byTag['h2_bwa'][1], byTag['h2_bwa'][2]
                )
            }
    combined_kept = COMBINE_KEPT_FILTERED_FASTQS(combine_in)
    // emits: tuple(sid, sid_kept_R1.fastq.gz, sid_kept_R2.fastq.gz)

    //combine the filtered counts
    filt_counts_host1_minimap2 = host1_minimap2_filter_out.map { sid, _kept_r1, _kept_r2, _unm1, _unm2, chrom, total ->
        tuple(sid, chrom, total)
    }

    filt_counts_host1_bwa = host1_bwa_filter_out.map { sid, _kept_r1, _kept_r2, _unm1, _unm2, chrom, total ->
        tuple(sid, chrom, total)
    }

    filt_counts_host2_minimap2 = host2_minimap2_filter_out.map { sid, _kept_r1, _kept_r2, _unm1, _unm2, chrom, total ->
        tuple(sid, chrom, total)
    }

    filt_counts_host2_bwa = host2_bwa_filter_out.map { sid, _kept_r1, _kept_r2, _unm1, _unm2, chrom, total ->
        tuple(sid, chrom, total)
    }
    filt_counts_host1_minimap2_tagged = filt_counts_host1_minimap2.map { sid, chrom, total -> tuple(sid, 'h1_mm',  chrom, total) }
    filt_counts_host1_bwa_tagged      = filt_counts_host1_bwa.map      { sid, chrom, total -> tuple(sid, 'h1_bwa', chrom, total) }
    filt_counts_host2_minimap2_tagged = filt_counts_host2_minimap2.map { sid, chrom, total -> tuple(sid, 'h2_mm',  chrom, total) }
    filt_counts_host2_bwa_tagged      = filt_counts_host2_bwa.map      { sid, chrom, total -> tuple(sid, 'h2_bwa', chrom, total) }

    all_filt_counts =
        filt_counts_host1_minimap2_tagged
            .mix(filt_counts_host1_bwa_tagged)
            .mix(filt_counts_host2_minimap2_tagged)
            .mix(filt_counts_host2_bwa_tagged)

    // group by sample_id; each record is [tag, chrom, total]
    combine_filt_counts_in =
        all_filt_counts
            .map { sid, tag, chrom, total -> tuple(sid, [tag, chrom, total]) }
            .groupTuple(size: 4)
            .map { sid, recs ->
                def byTag = recs.collectEntries { r -> [(r[0]): r] }
                tuple(
                    sid, fastp_trim_info_dir,
                    byTag['h1_mm'][1],  byTag['h1_mm'][2],
                    byTag['h1_bwa'][1], byTag['h1_bwa'][2],
                    byTag['h2_mm'][1],  byTag['h2_mm'][2],
                    byTag['h2_bwa'][1], byTag['h2_bwa'][2]
                )
            }
    combined_counts = COMBINE_KEPT_FILTERED_COUNTS(combine_filt_counts_in)

    // map to final bam
    final_map_in = combined_kept.map { sid, kept_r1, kept_r2 ->
        tuple(
            sid,
            wgs_final_genome_fasta_fnp,
            file("${wgs_final_genome_fasta_fnp}.*"),
            kept_r1,
            kept_r2
        )
    }

    final_bam_out = FINAL_MAP_BWA(final_map_in)

    //summarize
    fastp_json_tagged = trimmed_ch.map { sid, _r1, _r2, json ->
        tuple(sid, 'fastp', json)
    }

    total_counts_tagged = combined_counts.map { sid, _chrom, total ->
        tuple(sid, 'counts', total)
    }

    final_bam_tagged = final_bam_out.map { sid, bam, _bai ->
        tuple(sid, 'bam', bam)
    }

    all_summary_inputs =
        fastp_json_tagged
            .mix(total_counts_tagged)
            .mix(final_bam_tagged)

    summary_in =
        all_summary_inputs
            .map { sid, tag, obj -> tuple(sid, [tag, obj]) }
            .groupTuple(size: 3)
            .map { sid, recs ->
                def byTag = recs.collectEntries { r -> [(r[0]): r[1]] }
                tuple(sid, fastp_trim_info_dir, byTag['fastp'], byTag['counts'], byTag['bam'])
            }

    WGS_BUILD_FINAL_SAMPLE_SUMMARY(summary_in)

}
