
process WRITE_SAMPLE_NAME_KEY {
  tag "sample_name_key"

  publishDir { meta_dir }, mode: 'copy', overwrite: true

  input:
    val lines
    val meta_dir

  output:
    path "sample_name_key.tsv"

  script:
    def header = "old_sample_name\tnew_sample_name\n"

    // sort by old_sample_name (first column)
    def sorted = (lines ?: []).sort { a, b ->
      def aOld = a.split('\t', 2)[0]
      def bOld = b.split('\t', 2)[0]
      aOld <=> bOld
    }

    def body = sorted.join('\n')
    def text = header + (body ? body : "")

    """
    cat > sample_name_key.tsv <<'EOF'
${text}
EOF
    """
}
