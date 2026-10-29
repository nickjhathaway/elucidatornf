process SCAN_RAW_INPUT_SAMPLE_NAMES {
    tag "scan_raw_input_sample_names"
    label 'process_single'
    // a name collision fails the same way every time, so don't retry it
    errorStrategy 'terminate'

    publishDir { "${pub_dir}" }, mode: 'copy', overwrite: true

    input:
    path input_dir // directory of fastq(.gz) and/or unaligned bam files
    val technology
    val recursive // search input_dir recursively
    val prefix_with_top_dir // prefix names with the top directory under input_dir, e.g. run1/sub/barcode01.bam -> run1_barcode01
    val use_bam_sm // name bams by the SM of their read group rather than their file name
    val pub_dir

    output:
    // columns: file (relative to input_dir), sample, baseName, nameSource, topDir, status (ok, ignored or duplicate)
    path "rawInputSampleNames.tsv", emit: sample_names

    script:
    def extra_args = task.ext.args ? task.ext.args : ''
    def recursive_arg = recursive ? "--recursive" : ""
    def prefix_arg = prefix_with_top_dir ? "--prefixWithTopDir" : ""
    def bam_sm_arg = use_bam_sm ? "" : "--doNotUseBamSM"

    // only bam headers are read, so this finds sample name collisions before the long extraction step
    """
    SeekDeep getPossibleSampleNamesFromRawInput \
            --inputDir ${input_dir} \
            --technology ${technology} \
            ${recursive_arg} \
            ${prefix_arg} \
            ${bam_sm_arg} \
            --failOnDuplicates \
            --out rawInputSampleNames.tsv \
            ${extra_args}
    """
}
