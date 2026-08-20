#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Host read depletion only — for submitting otherwise-raw data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This is the host-filter core of run_wgs_preprocess_pf.nf with fastp trimming and the
    final Pf mapping removed. Reads that survive the filter are byte-identical to the
    input reads (BamFilterByChroms restores the original orientation off the bam), so the
    only thing that happens to the data is removal of host reads and of any singletons
    the filtering creates.
----------------------------------------------------------------------------------------
*/

include { MAP_BWA as HOST1_FILT_MAP_BWA } from '../../../modules/local/map_bwa'
include { MAP_BWA as HOST2_FILT_MAP_BWA } from '../../../modules/local/map_bwa'
include { MAP_MINIMAP2 as HOST1_FILT_MAP_MINIMAP2 } from '../../../modules/local/map_minimap2'
include { MAP_MINIMAP2 as HOST2_FILT_MAP_MINIMAP2 } from '../../../modules/local/map_minimap2'

include { BAM_FILTER_BY_CHROMS as HOST1_MINIMAP2_BAM_FILTER_BY_CHROMS } from '../../../modules/local/bam_filter_by_chroms'
include { BAM_FILTER_BY_CHROMS as HOST1_BWA_BAM_FILTER_BY_CHROMS } from '../../../modules/local/bam_filter_by_chroms'
include { BAM_FILTER_BY_CHROMS as HOST2_MINIMAP2_BAM_FILTER_BY_CHROMS } from '../../../modules/local/bam_filter_by_chroms'
include { BAM_FILTER_BY_CHROMS as HOST2_BWA_BAM_FILTER_BY_CHROMS } from '../../../modules/local/bam_filter_by_chroms'

include { COMBINE_KEPT_FILTERED_FASTQS } from '../../../modules/local/combine_kept_filtered_fastqs'
include { COMBINE_KEPT_FILTERED_COUNTS } from '../../../modules/local/combine_kept_filtered_fastqs'

include { SYNC_PAIRED_FASTQS } from '../../../modules/local/sync_paired_fastqs'

include { WGS_HOST_DEPLETION_SAMPLE_SUMMARY } from '../../../modules/local/wgs_host_depletion_summary'
include { WGS_HOST_DEPLETION_RUN_SUMMARY } from '../../../modules/local/wgs_host_depletion_summary'

include { WRITE_SAMPLE_NAME_KEY } from '../../../modules/local/write_sample_name_key/main.nf'


workflow RUN_WGS_HOST_DEPLETION_PF {
    take:
    input_fastq_dir // a directory with all fastqs to be analyzed
    output_dir      // the output directory
    wgs_host1_filter_genome_fasta_fnp
    wgs_host1_filter_contigs_fnp
    wgs_host2_filter_genome_fasta_fnp
    wgs_host2_filter_contigs_fnp

    main:

    /********************************
    * Output layout
    ********************************/
    def depleted_fastq_dir = file("${output_dir}/host_depleted_fastqs")
    depleted_fastq_dir.mkdirs()
    def filter_info_dir = file("${output_dir}/host_filter_info")
    filter_info_dir.mkdirs()
    def unpaired_dir = file("${output_dir}/discarded_singletons")
    if( params.wgs_depletion_save_unpaired ) {
        unpaired_dir.mkdirs()
    }

    /********************************
    * Helper: build rename map
    ********************************/
    def rename_map = SampleRename.loadRenameKey(params.wgs_rename_tsv)

    /********************************
    * Input channel: paired fastqs
    * supports *_1.fastq.gz/_2.fastq.gz and *_R1.fastq.gz/_R2.fastq.gz
    ********************************/
    def keep_set = null
    if( params.wgs_samples_keep_file ) {
        keep_set = new HashSet<String>()
        new File(params.wgs_samples_keep_file).eachLine { line ->
            def s = line.trim()
            if( !s || s.startsWith('#') ) return
            keep_set.add(s)
        }
    }

    all_samples_ch = channel
        .fromFilePairs(
            "${input_fastq_dir}/*_{R,}{1,2}*.{fq,fastq}.gz",
            flat: true
        )
        .filter { sid, _r1, _r2 ->
            // if no keep file, keep everything
            keep_set == null || keep_set.contains(sid)
        }
        .map { sid, r1, r2 ->
            def sample_id = sid
            def final_id  = rename_map.get(sample_id, sample_id)
            tuple(sample_id, final_id, [r1, r2])
        }
        .unique()

    // write the old -> new sample name key for the whole batch, including samples that
    // get skipped below because they were already depleted on an earlier run
    sample_key_lines_ch = all_samples_ch
        .map { sample_id, final_id, _reads -> "${sample_id}\t${final_id}" }
        .distinct()
        .collect()

    WRITE_SAMPLE_NAME_KEY(sample_key_lines_ch, filter_info_dir.toString())

    READS_CH = all_samples_ch
        .map { sample_id, final_id, reads ->
            // skip if the depleted R1 already exists and skip_existing is true
            def final_r1 = file("${depleted_fastq_dir}/${final_id}_R1.fastq.gz")
            def skip = params.wgs_skip_existing && final_r1.exists()
            tuple(sample_id, final_id, reads, skip)
        }
        .filter { _sample_id, _final_id, _reads, skip -> !skip }
        .map { _sample_id, final_id, reads, _skip ->
            tuple(final_id, reads[0], reads[1])
        }

    // set up paths for minimap2 index files
    def host1_mmi = wgs_host1_filter_genome_fasta_fnp
        .toString()
        .replaceFirst(/\.(fa|fasta|fna)(\.gz)?$/, '.mmi')
    def host1_mmi_path = file(host1_mmi)

    def host2_mmi = wgs_host2_filter_genome_fasta_fnp
        .toString()
        .replaceFirst(/\.(fa|fasta|fna)(\.gz)?$/, '.mmi')
    def host2_mmi_path = file(host2_mmi)

    // map minimap2 to host 1 filter — raw reads straight in, no trimming
    host1_minimap2_in = READS_CH.map { sid, r1, r2 ->
        tuple(sid, host1_mmi_path, r1, r2)
    }
    host1_minimap2_out = HOST1_FILT_MAP_MINIMAP2(host1_minimap2_in)
    // filter by host 1 filter on minimap2
    host1_minimap2_filt_in = host1_minimap2_out.map { sid, bam_fnp, bam_bai_fnp ->
        tuple(sid, "host1_minimap2_filt", bam_fnp, bam_bai_fnp, wgs_host1_filter_contigs_fnp, true)
    }
    host1_minimap2_filter_out = HOST1_MINIMAP2_BAM_FILTER_BY_CHROMS(host1_minimap2_filt_in)

    // map bwa to host 1 filter
    host1_bwa_map_in = host1_minimap2_filter_out.map { sid, _kept_r1, _kept_r2, unmapped_r1, unmapped_r2, _chrom_tab, _total_tab ->
        tuple(sid, wgs_host1_filter_genome_fasta_fnp, file("${wgs_host1_filter_genome_fasta_fnp}.*"), unmapped_r1, unmapped_r2)
    }
    host1_bwa_map_out = HOST1_FILT_MAP_BWA(host1_bwa_map_in)
    // filter by host 1 filter on bwa
    host1_bwa_filt_in = host1_bwa_map_out.map { sid, bam_fnp, bam_bai_fnp, _flagstat ->
        tuple(sid, "host1_bwa_filt", bam_fnp, bam_bai_fnp, wgs_host1_filter_contigs_fnp, true)
    }
    host1_bwa_filter_out = HOST1_BWA_BAM_FILTER_BY_CHROMS(host1_bwa_filt_in)

    // map minimap2 to host 2 filter
    host2_minimap2_in = host1_bwa_filter_out.map { sid, _kept_r1, _kept_r2, unmapped_r1, unmapped_r2, _chrom_tab, _total_tab ->
        tuple(sid, host2_mmi_path, unmapped_r1, unmapped_r2)
    }
    host2_minimap2_out = HOST2_FILT_MAP_MINIMAP2(host2_minimap2_in)
    // filter by host 2 filter on minimap2
    host2_minimap2_filt_in = host2_minimap2_out.map { sid, bam_fnp, bam_bai_fnp ->
        tuple(sid, "host2_minimap2_filt", bam_fnp, bam_bai_fnp, wgs_host2_filter_contigs_fnp, true)
    }
    host2_minimap2_filter_out = HOST2_MINIMAP2_BAM_FILTER_BY_CHROMS(host2_minimap2_filt_in)

    // map bwa to host 2 filter
    host2_bwa_map_in = host2_minimap2_filter_out.map { sid, _kept_r1, _kept_r2, unmapped_r1, unmapped_r2, _chrom_tab, _total_tab ->
        tuple(sid, wgs_host2_filter_genome_fasta_fnp, file("${wgs_host2_filter_genome_fasta_fnp}.*"), unmapped_r1, unmapped_r2)
    }
    host2_bwa_map_out = HOST2_FILT_MAP_BWA(host2_bwa_map_in)
    // filter by host 2 filter on bwa — last pass, so the still-unmapped reads stay in kept
    host2_bwa_filt_in = host2_bwa_map_out.map { sid, bam_fnp, bam_bai_fnp, _flagstat ->
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

    // drop any singleton left over from the filtering so the submitted data stays
    // strictly paired end
    sync_in = combined_kept.map { sid, kept_r1, kept_r2 ->
        tuple(sid, depleted_fastq_dir.toString(), unpaired_dir.toString(), kept_r1, kept_r2)
    }
    synced = SYNC_PAIRED_FASTQS(sync_in)

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
                    sid, filter_info_dir,
                    byTag['h1_mm'][1],  byTag['h1_mm'][2],
                    byTag['h1_bwa'][1], byTag['h1_bwa'][2],
                    byTag['h2_mm'][1],  byTag['h2_mm'][2],
                    byTag['h2_bwa'][1], byTag['h2_bwa'][2]
                )
            }
    combined_counts = COMBINE_KEPT_FILTERED_COUNTS(combine_filt_counts_in)

    // per sample summary of what was filtered off, then one table for the whole run
    summary_in = combined_counts
        .map { sid, _chrom, total -> tuple(sid, total) }
        .join(synced.stats)

    sample_summaries = WGS_HOST_DEPLETION_SAMPLE_SUMMARY(summary_in)

    WGS_HOST_DEPLETION_RUN_SUMMARY(sample_summaries.collect(), filter_info_dir.toString())

    emit:
    paired_fastqs = synced.paired
    sample_summary = sample_summaries
}
