# Bash completion for nickjhathaway/elucidatornf workflows  (GENERATED -- see
# completions/generate_completion.sh; edits here are overwritten on regen)
#
# Completes `nextflow run <...>/<workflow>.nf` with the flags that belong to the
# detected workflow (by basename, so full or partial paths both work), plus
# file/dir completion for path-valued flags and the workflow script itself.
# For any command that is NOT one of these workflows, it delegates to nextflow's
# own completion if one was registered before this file was sourced.
#
# Install (any one of):
#   * echo 'source /path/to/elucidatornf/completions/elucidatornf.bash' >> ~/.bashrc
#   * cp completions/elucidatornf.bash ~/.local/share/bash-completion/completions/nextflow
#   * sudo cp completions/elucidatornf.bash /etc/bash_completion.d/elucidatornf

_ELU_NXF='-profile -resume -params-file -c -config -w -work-dir -with-report -with-trace -with-timeline -with-dag -ansi-log -bg -name -stub-run -dump-channels -preview'
_ELU_COMMON='--outdir --publish_dir_mode --max_task_cpus --max_task_memory --max_task_time --max_executor_cpus --max_executor_memory --max_retry --max_cpus --max_time'
_ELU_VARIANT='--do_variant_calling --variant_calling_ncpus --known_amino_acid_changes_fnp --vc_variant_frequency_cut_off --vc_variant_occurrence_cut_off --vc_no_pairwise_comps --vc_primary_genome --vc_correct_small_homopolymer_errors --meta_fields_to_calc_pop_diffs --vc_extra_args --meta_fnp'
_ELU_AMPLICON='--primers_fnp --genome_dir --gff_dir --nanopore_primers_errors_allowed --nanopore_rename_key_fnp --nanopore_render_clustering_report --amplicon_clustering_pop_clustering_extra_args --amplicon_clustering_test_number --amplicon_clustering_number_of_files_to_investigate'
_ELU_WGS_COMMON='--wgs_samples_keep_file --wgs_rename_tsv --wgs_skip_existing --wgs_host_filter_any_mate --wgs_host_filter_filter_with_unmapped_mate --wgs_host_filter_min_mapq'
_ELU_ILLUMINA='--illumina_clustering_paired_end_length --illumina_clustering_error_profile --illumina_clustering_min_sample_read_count --illumina_clustering_population_ncpus --illumina_clustering_trim_front --illumina_clustering_trim_back'
_ELU_NANOPORE='--nanopore_clustering_min_len --nanopore_clustering_min_sample_read_count --nanopore_clustering_ncpus --nanopore_clustering_trim_front --nanopore_clustering_trim_back --nanopore_clustering_max_reads_use --nanopore_extractor_use_inner_primers --nanopore_clustering_lower_base --nanopore_extractor_primer_within_start'
_ELU_PATHWEAVER='--pw_ncpus --pw_pop_clustering_ncpus --pw_seqs_table --pw_seqs_table_seqs_col --pw_seqs_table_name_col --pw_seqs_table_target_col --pw_gene_ids --pw_samples_file --pw_bed_fnp --pw_extract_region_assemblies_extra_args --genome_fnp --bams_dir --bams_file_ending --render_pw_report --pw_pop_clustering_extra_args --pw_overwrite_dir --run_sub_segments_determination --sub_segments_correction_occurence_cut_off --sub_segments_low_freq_cut_off --sub_segments_uniqueHapCountCutOff --sub_segments_extra_args --primers_fnp'
_ELU_RNASEQ='--rnaseq_host1_filter_genome_fasta_fnp --rnaseq_host1_filter_contigs_fnp --rnaseq_host2_filter_genome_fasta_fnp --rnaseq_host2_filter_contigs_fnp --rnaseq_salmon_index --single_cell_seurat_r_object_fnp'
_ELU_WGS='--wgs_host1_filter_genome_fasta_fnp --wgs_host1_filter_contigs_fnp --wgs_host2_filter_genome_fasta_fnp --wgs_host2_filter_contigs_fnp --wgs_final_genome_fasta_fnp'

# file/dir completion that works with or without the bash-completion package
_elu_filedir() {
    if declare -F _filedir >/dev/null 2>&1; then
        _filedir "$@"
    elif [[ "${1:-}" == nf ]]; then
        local IFS=$'\n'
        COMPREPLY=( $(compgen -d -- "$cur") $(compgen -f -X '!*.nf' -- "$cur") )
    else
        local IFS=$'\n'
        COMPREPLY=( $(compgen -f -- "$cur") )
    fi
}

# delegate to a pre-existing `nextflow` completion (captured at load time, see
# bottom of file); otherwise fall back to plain file/dir completion
_elu_fallback() {
    if [[ -n "${_ELU_PREV_FUNC:-}" ]] && declare -F "$_ELU_PREV_FUNC" >/dev/null 2>&1; then
        "$_ELU_PREV_FUNC" "$@"
    else
        _elu_filedir
    fi
}

_elucidatornf_complete() {
    local cur prev w base wf=""
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"

    # detect an elucidatornf workflow on the line (by basename, path-agnostic)
    for w in "${COMP_WORDS[@]:1}"; do
        base="${w##*/}"
        case "$base" in
            illumina_clustering.nf)        wf=illumina ;;
            nanopore_clustering.nf)        wf=nanopore ;;
            pathweaver_extract_regions.nf) wf=pathweaver ;;
            wgs_preprocess_pf.nf)          wf=wgs ;;
            rnaseq_process_pf.nf)          wf=rnaseq ;;
            index_genomes.nf)              wf=index ;;
        esac
    done

    # not one of ours -> hand off to nextflow's own completion if present
    if [[ -z "$wf" ]]; then
        _elu_fallback "$@"
        return
    fi

    # value completion for the previous flag
    case "$prev" in
        -profile)
            COMPREPLY=( $(compgen -W "standard docker singularity apptainer test slurm sge" -- "$cur") ); return ;;
        --*_fnp|--*_dir|--outdir|--*_tsv|--*_file|--*_index|--pw_seqs_table|--pw_gene_ids|--genome_fnp|\
        -params-file|-c|-config|-w|-work-dir|-with-report|-with-trace|-with-timeline|-with-dag)
            _elu_filedir; return ;;
    esac

    # candidate flags for the detected workflow
    local flags
    case "$wf" in
        illumina)   flags="$_ELU_COMMON --input_fastq_dir $_ELU_VARIANT $_ELU_AMPLICON $_ELU_ILLUMINA $_ELU_NXF" ;;
        nanopore)   flags="$_ELU_COMMON --input_fastq_dir $_ELU_VARIANT $_ELU_AMPLICON $_ELU_NANOPORE $_ELU_NXF" ;;
        pathweaver) flags="$_ELU_COMMON $_ELU_VARIANT $_ELU_PATHWEAVER $_ELU_NXF" ;;
        wgs)        flags="$_ELU_COMMON --input_fastq_dir $_ELU_WGS_COMMON $_ELU_WGS $_ELU_NXF" ;;
        rnaseq)     flags="$_ELU_COMMON --input_fastq_dir $_ELU_WGS_COMMON $_ELU_RNASEQ $_ELU_NXF" ;;
        index)      flags="$_ELU_COMMON --genome_dir $_ELU_NXF" ;;
        *)          flags="$_ELU_NXF" ;;
    esac

    if [[ "$cur" == -* ]]; then
        COMPREPLY=( $(compgen -W "$flags" -- "$cur") )
    else
        _elu_filedir nf
    fi
}

# Capture any pre-existing `nextflow` completion BEFORE we register ours, so the
# handler above can delegate to it for non-elucidatornf invocations. For this to
# find it, source this file AFTER nextflow's own completion is registered (only
# `complete -F <func>` style is proxied).
_elu_prev_spec="$(complete -p nextflow 2>/dev/null || true)"
if [[ "$_elu_prev_spec" =~ -F[[:space:]]+([^[:space:]]+) && "${BASH_REMATCH[1]}" != _elucidatornf_complete ]]; then
    _ELU_PREV_FUNC="${BASH_REMATCH[1]}"
fi
: "${_ELU_PREV_FUNC:=}"
unset _elu_prev_spec
complete -F _elucidatornf_complete nextflow
