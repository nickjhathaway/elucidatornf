#!/usr/bin/env bash
# Regenerate completions/elucidatornf.bash from conf/params/*.config so the
# completion flag lists always mirror the actual params.
#
#   ./completions/generate_completion.sh        # writes completions/elucidatornf.bash
set -euo pipefail
cd "$(dirname "$0")/.."   # repo root
OUT="completions/elucidatornf.bash"

# extract "--name --name ..." from a conf/params file, dropping internal params
extract() {
    grep -oE '^[[:space:]]+[a-zA-Z_][a-zA-Z0-9_]*[[:space:]]*=' "conf/params/$1.config" \
      | sed -E 's/[[:space:]]//g; s/=//' \
      | grep -avE '^(help|help_full|show_hidden|version|validate_params|monochrome_logs|trace_timestamp|machine_cpus|machine_memory|empty_file_fnp|config_profile_name|config_profile_description|config_profile_contact|config_profile_url|custom_config_version)$' \
      | sed 's/^/--/' | tr '\n' ' ' | sed -E 's/ +$//'
}

GLOBAL_ALL=$(extract global)
COMMON=$(echo "$GLOBAL_ALL" | tr ' ' '\n' | grep -avx -- '--input_fastq_dir' | tr '\n' ' ' | sed -E 's/ +$//')
VARIANT=$(extract variant_calling)
AMPLICON=$(extract amplicon)
WGS_COMMON=$(extract wgs_common)
ILLUMINA=$(extract illumina)
NANOPORE=$(extract nanopore)
PATHWEAVER="$(extract pathweaver) --primers_fnp"
RNASEQ=$(extract rnaseq)
WGS=$(extract wgs)
WGS_DEPLETION=$(extract wgs_host_depletion)

{
cat <<'HEADER'
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

HEADER

echo "_ELU_NXF='-profile -resume -params-file -c -config -w -work-dir -with-report -with-trace -with-timeline -with-dag -ansi-log -bg -name -stub-run -dump-channels -preview'"
echo "_ELU_COMMON='${COMMON}'"
echo "_ELU_VARIANT='${VARIANT}'"
echo "_ELU_AMPLICON='${AMPLICON}'"
echo "_ELU_WGS_COMMON='${WGS_COMMON}'"
echo "_ELU_ILLUMINA='${ILLUMINA}'"
echo "_ELU_NANOPORE='${NANOPORE}'"
echo "_ELU_PATHWEAVER='${PATHWEAVER}'"
echo "_ELU_RNASEQ='${RNASEQ}'"
echo "_ELU_WGS='${WGS}'"
echo "_ELU_WGS_DEPLETION='${WGS_DEPLETION}'"

cat <<'BODY'

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
            wgs_host_depletion_pf.nf)      wf=wgs_depletion ;;
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
        wgs_depletion) flags="$_ELU_COMMON --input_fastq_dir $_ELU_WGS_COMMON $_ELU_WGS $_ELU_WGS_DEPLETION $_ELU_NXF" ;;
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
BODY
} > "$OUT"

echo "wrote $OUT"
