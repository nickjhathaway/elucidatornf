#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { FASTP_TRIM } from '../../../modules/local/run_fastp_pf_rna_seq'
include { MAP_BWA as HOST1_FILT_MAP_BWA} from '../../../modules/local/map_bwa'
include { MAP_BWA as HOST2_FILT_MAP_BWA} from '../../../modules/local/map_bwa'
include { MAP_MINIMAP2 as HOST1_FILT_MAP_MINIMAP2 } from '../../../modules/local/map_minimap2'
include { MAP_MINIMAP2 as HOST2_FILT_MAP_MINIMAP2} from '../../../modules/local/map_minimap2'

include { BAM_FILTER_BY_CHROMS as HOST1_MINIMAP2_BAM_FILTER_BY_CHROMS } from '../../../modules/local/bam_filter_by_chroms'
include { BAM_FILTER_BY_CHROMS as HOST1_BWA_BAM_FILTER_BY_CHROMS} from '../../../modules/local/bam_filter_by_chroms'
include { BAM_FILTER_BY_CHROMS as HOST2_MINIMAP2_BAM_FILTER_BY_CHROMS} from '../../../modules/local/bam_filter_by_chroms'
include { BAM_FILTER_BY_CHROMS as HOST2_BWA_BAM_FILTER_BY_CHROMS} from '../../../modules/local/bam_filter_by_chroms'

include { COMBINE_KEPT_FILTERED_FASTQS } from '../../../modules/local/combine_kept_filtered_fastqs'
include { COMBINE_KEPT_FILTERED_COUNTS } from '../../../modules/local/combine_kept_filtered_fastqs'

include { SALMON_QUANT as FINAL_SALMON_QUANT } from '../../../modules/local/salmon_quant'
include { SALMON_QUANT as INITIAL_SALMON_QUANT } from '../../../modules/local/salmon_quant'
include { SALMON_QUANT_COMBINE as FINAL_SALMON_QUANT_COMBINE} from '../../../modules/local/salmon_quant_combine'
include { SALMON_QUANT_COMBINE as INITIAL_SALMON_QUANT_COMBINE} from '../../../modules/local/salmon_quant_combine'

include { RUN_MUSIC_DECOMP_PF } from '../../../modules/local/running_music_decomp_pf'
include { RUN_INSTAPRISM_DECONVOLUTION_PF } from '../../../modules/local/running_instaprism_deconvolution_pf'

include { KMER_SET_FILTER } from '../../../modules/local/extract_by_kmer_sets'
include { KMER_SET_SORT } from '../../../modules/local/extract_by_kmer_sets'
include { COMBINE_KMER_EXTRACTION_COUNTS as COMBINE_KMER_FILTER_COUNTS } from '../../../modules/local/extract_by_kmer_sets'
include { COMBINE_KMER_EXTRACTION_COUNTS as COMBINE_KMER_SORT_COUNTS } from '../../../modules/local/extract_by_kmer_sets'

include { RNASEQ_SAMPLE_SUMMARY } from '../../../modules/local/rnaseq_summary'
include { RNASEQ_RUN_SUMMARY } from '../../../modules/local/rnaseq_summary'
include { RNASEQ_RESULTS_README } from '../../../modules/local/rnaseq_summary'


// Published per-sample outputs that the cohort-level steps (quant combining,
// deconvolution, count combining, run summary) need from every sample, including
// samples processed by an earlier run into the same output directory.
def rnaseqPublishedOutputs(output_dir, sample_id) {
    def counts_dir = "${output_dir}/kmer_extraction/counts"
    def count_tabs = ['extractionCounts', 'multihitCounts', 'inversePairsCounts']
    return [
        initial_quant      : file("${output_dir}/initial_quants/${sample_id}_salmon_quants_quant.sf.gz"),
        filtered_quant     : file("${output_dir}/filtered_quants/${sample_id}_salmon_quants_quant.sf.gz"),
        kmer_filter_counts : count_tabs.collect { tab -> file("${counts_dir}/${sample_id}_kmer_filter_${tab}.tab.txt") },
        kmer_sort_counts   : count_tabs.collect { tab -> file("${counts_dir}/${sample_id}_kmer_sort_${tab}.tab.txt") },
        summary            : file("${output_dir}/summary/${sample_id}_rnaseq_summary.tab.txt"),
    ]
}

// A sample is skipped only when every output above exists; anything less reruns it in full.
def rnaseqSampleComplete(output_dir, sample_id) {
    def outs = rnaseqPublishedOutputs(output_dir, sample_id)
    def all_files = [outs.initial_quant, outs.filtered_quant, outs.summary] + outs.kmer_filter_counts + outs.kmer_sort_counts
    return all_files.every { f -> f.exists() }
}


workflow RUN_RNASEQ_PROCESS_PF{
    take:
    input_fastq_dir //a directory with all fastqs to be analyzed
    output_dir // the output directory
    rnaseq_host1_filter_genome_fasta_fnp
    rnaseq_host1_filter_contigs_fnp
    rnaseq_host2_filter_genome_fasta_fnp
    rnaseq_host2_filter_contigs_fnp
    salmon_index
    single_cell_seurat_r_object_fnp
    kmer_filter_table_fnp // k-mer sets whose reads are removed before sorting (e.g. mycoplasma, Pf rRNA)
    kmer_sort_table_fnp   // k-mer sets to sort the remaining reads into (e.g. vars, rifins, stevors)
    main:

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

    /*
    Create output directories
    Publish dirs are handed to processes as strings: a Path val is hashed like a
    file, by mtime, so publishing into the dir would invalidate -resume for every
    task that received it.
    */
    def fastp_trim_info_dir = "${output_dir}/fastp_and_filter_info"
    file(fastp_trim_info_dir).mkdirs()

    def human_removed_fastqs_dir = file("${output_dir}/human_removed_fastqs")
    human_removed_fastqs_dir.mkdirs()

    // CLI params arrive as strings, so a plain truthiness test would treat "false" as on
    def skip_existing = params.wgs_skip_existing.toString().toBoolean()

    samples_ch = channel
        .fromFilePairs(
            "${input_fastq_dir}/*_{R,}{1,2}*.{fq,fastq}.gz",
            flat: true
        )
        .filter { sid, _r1, _r2 ->
            // if no keep file, keep everything
            keep_set == null || keep_set.contains(sid)
        }
        .map { sid, r1, r2 ->
            tuple(rename_map.get(sid, sid), [r1, r2])
        }
        .unique()
        .branch { final_id, _reads ->
            done: skip_existing && rnaseqSampleComplete(output_dir, final_id)
            todo: true
        }

    READS_CH = samples_ch.todo.map { final_id, reads ->
        tuple(final_id, reads, fastp_trim_info_dir)
    }

    // processed by an earlier run into this output dir: only their published
    // per-sample outputs are fed into the cohort-level steps below
    samples_ch.done.subscribe { final_id, _reads ->
        log.info "skipping ${final_id}, already processed in ${output_dir}"
    }
    done_outputs_ch = samples_ch.done.map { final_id, _reads ->
        rnaseqPublishedOutputs(output_dir, final_id)
    }

    // set up paths for minimap2 index files
    def host1_mmi = rnaseq_host1_filter_genome_fasta_fnp
        .toString()
        .replaceFirst(/\.(fa|fasta|fna)(\.gz)?$/, '.mmi')
    def host1_mmi_path = file(host1_mmi)

    def host2_mmi = rnaseq_host2_filter_genome_fasta_fnp
        .toString()
        .replaceFirst(/\.(fa|fasta|fna)(\.gz)?$/, '.mmi')
    def host2_mmi_path = file(host2_mmi)

    // run fastp
    fastp_ch = FASTP_TRIM(READS_CH)
    trimmed_ch = fastp_ch.trimmed

    // initial quant sf
    def initial_quants_dir = "${output_dir}/initial_quants"
    file(initial_quants_dir).mkdirs()
    initial_salmon_in = trimmed_ch.map { sid, r1_trim, r2_trim, _json ->
        tuple(sid, salmon_index, r1_trim, r2_trim, initial_quants_dir)
    }
    //combine initial quants
    INITIAL_SALMON_QUANT(initial_salmon_in)
    initial_salmon_combine_in = INITIAL_SALMON_QUANT.out
        .map { _sample_id, quant_sf, _lib_format, _meta_info -> quant_sf }
        .mix(done_outputs_ch.map { outs -> outs.initial_quant })
        .collect()
        .map {all_quant_fnps -> tuple(all_quant_fnps, initial_quants_dir) }

    INITIAL_SALMON_QUANT_COMBINE(initial_salmon_combine_in)

    // map minimap2 to host 1 filter
    host1_minimap2_in = trimmed_ch.map { sid, r1_trim, r2_trim, _json ->
        tuple(sid, host1_mmi_path, r1_trim, r2_trim)
    }
    host1_minimap2_out = HOST1_FILT_MAP_MINIMAP2(host1_minimap2_in)
    // filter by host 1 filter on minimap2
    host1_minimap2_filt_in = host1_minimap2_out.map{sid, bam_fnp, bam_bai_fnp ->
        tuple(sid, "host1_minimap2_filt", bam_fnp, bam_bai_fnp, rnaseq_host1_filter_contigs_fnp, true)
    }
    host1_minimap2_filter_out = HOST1_MINIMAP2_BAM_FILTER_BY_CHROMS(host1_minimap2_filt_in)

    // map bwa to host 1 filter
    host1_bwa_map_in = host1_minimap2_filter_out.map{sid, _kept_r1, _kept_r2, unmapped_r1, unmapped_r2, _chrom_tab, _total_tab ->
        tuple(sid, rnaseq_host1_filter_genome_fasta_fnp, files("${rnaseq_host1_filter_genome_fasta_fnp}.*"), unmapped_r1, unmapped_r2)
    }
    host1_bwa_map_out = HOST1_FILT_MAP_BWA(host1_bwa_map_in)
    // filter by host 1 filter on bwa
    host1_bwa_filt_in = host1_bwa_map_out.map{sid, bam_fnp, bam_bai_fnp, _flagstat ->
        tuple(sid, "host1_bwa_filt", bam_fnp, bam_bai_fnp, rnaseq_host1_filter_contigs_fnp, true)
    }
    host1_bwa_filter_out = HOST1_BWA_BAM_FILTER_BY_CHROMS(host1_bwa_filt_in)

    // map minimap2 to host 2 filter
    host2_minimap2_in = host1_bwa_filter_out.map { sid, _kept_r1, _kept_r2, unmapped_r1, unmapped_r2, _chrom_tab, _total_tab ->
        tuple(sid, host2_mmi_path, unmapped_r1, unmapped_r2)
    }
    host2_minimap2_out = HOST2_FILT_MAP_MINIMAP2(host2_minimap2_in)
    // filter by host 2 filter on minimap2
    host2_minimap2_filt_in = host2_minimap2_out.map{sid, bam_fnp, bam_bai_fnp ->
        tuple(sid, "host2_minimap2_filt", bam_fnp, bam_bai_fnp, rnaseq_host2_filter_contigs_fnp, true)
    }
    host2_minimap2_filter_out = HOST2_MINIMAP2_BAM_FILTER_BY_CHROMS(host2_minimap2_filt_in)

    // map bwa to host 2 filter
    host2_bwa_map_in = host2_minimap2_filter_out.map{sid, _kept_r1, _kept_r2, unmapped_r1, unmapped_r2, _chrom_tab, _total_tab ->
        tuple(sid, rnaseq_host2_filter_genome_fasta_fnp, files("${rnaseq_host2_filter_genome_fasta_fnp}.*"), unmapped_r1, unmapped_r2)
    }
    host2_bwa_map_out = HOST2_FILT_MAP_BWA(host2_bwa_map_in)
    // filter by host 2 filter on bwa
    host2_bwa_filt_in = host2_bwa_map_out.map{sid, bam_fnp, bam_bai_fnp, _flagstat ->
        tuple(sid, "host2_bwa_filt", bam_fnp, bam_bai_fnp, rnaseq_host2_filter_contigs_fnp, false)
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
    //combined_counts =
    COMBINE_KEPT_FILTERED_COUNTS(combine_filt_counts_in)


    // quant filtered results
    def filtered_quants_dir = "${output_dir}/filtered_quants"
    file(filtered_quants_dir).mkdirs()
    final_salmon_in = combined_kept.map { sid, kept_r1, kept_r2 ->
        tuple(
            sid,
            salmon_index,
            kept_r1,
            kept_r2,
            filtered_quants_dir
        )
    }

    //combine the initial quants so summary can be run on them

    // run quant on filtered fastqs
    FINAL_SALMON_QUANT(final_salmon_in)

    // combining the quants and running music on it
    final_salmon_combine_in = FINAL_SALMON_QUANT.out
        .map { _sample_id, quant_sf, _lib_format, _meta_info -> quant_sf }
        .mix(done_outputs_ch.map { outs -> outs.filtered_quant })
        .collect()
        .map {all_quant_fnps -> tuple(all_quant_fnps, filtered_quants_dir) }

    FINAL_SALMON_QUANT_COMBINE(final_salmon_combine_in)
    RUN_MUSIC_DECOMP_PF(FINAL_SALMON_QUANT_COMBINE.out.combined_quants_fnp.map { combined_quants_fnp -> tuple(combined_quants_fnp, single_cell_seurat_r_object_fnp, filtered_quants_dir)})
    RUN_INSTAPRISM_DECONVOLUTION_PF(FINAL_SALMON_QUANT_COMBINE.out.combined_quants_fnp.map { combined_quants_fnp -> tuple(combined_quants_fnp, single_cell_seurat_r_object_fnp, filtered_quants_dir)})

    // k-mer set extraction on the host-depleted reads: first remove mycoplasma and
    // Pf rRNA reads, then sort what remains by gene/gene family
    def kmer_extraction_dir = "${output_dir}/kmer_extraction"
    file(kmer_extraction_dir).mkdirs()

    KMER_SET_FILTER(combined_kept.map { sid, kept_r1, kept_r2 ->
        tuple(sid, kept_r1, kept_r2, kmer_filter_table_fnp, kmer_extraction_dir)
    })
    KMER_SET_SORT(KMER_SET_FILTER.out.kept.map { sid, kept_r1, kept_r2 ->
        tuple(sid, kept_r1, kept_r2, kmer_sort_table_fnp, kmer_extraction_dir)
    })

    COMBINE_KMER_FILTER_COUNTS(KMER_SET_FILTER.out.counts
        .map { _sid, extraction, multihit, inverse_pairs -> [extraction, multihit, inverse_pairs] }
        .mix(done_outputs_ch.map { outs -> outs.kmer_filter_counts })
        .collect()
        .map { all_tabs -> tuple('kmer_filter', all_tabs, kmer_extraction_dir) })
    COMBINE_KMER_SORT_COUNTS(KMER_SET_SORT.out.counts
        .map { _sid, extraction, multihit, inverse_pairs -> [extraction, multihit, inverse_pairs] }
        .mix(done_outputs_ch.map { outs -> outs.kmer_sort_counts })
        .collect()
        .map { all_tabs -> tuple('kmer_sort', all_tabs, kmer_extraction_dir) })

    // read accounting and QC per sample: input -> fastp -> host / non-host ->
    // Pf, mycoplasma, undetermined -> sorted genes, plus salmon mapping rates
    def summary_dir = "${output_dir}/summary"
    file(summary_dir).mkdirs()

    sample_summary_in = trimmed_ch
        .map { sid, _r1_trim, _r2_trim, fastp_json -> tuple(sid, fastp_json) }
        .join(COMBINE_KEPT_FILTERED_COUNTS.out.map { sid, _by_chrom, total_counts -> tuple(sid, total_counts) }, failOnMismatch: true)
        .join(KMER_SET_FILTER.out.counts.map { sid, extraction, _multihit, _inverse_pairs -> tuple(sid, extraction) }, failOnMismatch: true)
        .join(KMER_SET_SORT.out.counts.map { sid, extraction, _multihit, _inverse_pairs -> tuple(sid, extraction) }, failOnMismatch: true)
        .join(INITIAL_SALMON_QUANT.out.map { sid, _quant_sf, lib_format, meta_info -> tuple(sid, meta_info, lib_format) }, failOnMismatch: true)
        .join(FINAL_SALMON_QUANT.out.map { sid, _quant_sf, _lib_format, meta_info -> tuple(sid, meta_info) }, failOnMismatch: true)
        .map { sid, fastp_json, total_counts, filter_counts, sort_counts, initial_meta, initial_lib_format, filtered_meta ->
            tuple(sid, fastp_json, total_counts, filter_counts, sort_counts, initial_meta, initial_lib_format, filtered_meta, summary_dir)
        }
    RNASEQ_SAMPLE_SUMMARY(sample_summary_in)

    RNASEQ_RUN_SUMMARY(RNASEQ_SAMPLE_SUMMARY.out
        .map { _sid, sample_summary -> sample_summary }
        .mix(done_outputs_ch.map { outs -> outs.summary })
        .collect()
        .map { summaries -> tuple(summaries, summary_dir) })

    // a README describing this results directory, written after the run summary
    RNASEQ_RESULTS_README(RNASEQ_RUN_SUMMARY.out.map { run_summary ->
        tuple(run_summary, file("${projectDir}/assets/rnaseq_results_README.md"), "${output_dir}")
    })

    // running assembly on the vars

    // running the domain determination on the assembled vars

}
