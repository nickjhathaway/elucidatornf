process WRITE_AMPLICON_SAMPLE_KEY {
  tag "sample_key"

  publishDir { meta_dir }, mode: 'copy', overwrite: true

  input:
    val lines // SampleKey.keyLine() lines
    val meta_dir

  output:
    path "sample_name_key.tsv"

  script:
    // sort by input then replicate
    def sorted = (lines ?: []).sort { a, b ->
      def aToks = a.split('\t', -1)
      def bToks = b.split('\t', -1)
      aToks[0] <=> bToks[0] ?: aToks[2] <=> bToks[2]
    }

    def body = sorted.join('\n')
    def text = SampleKey.keyHeader + "\n" + (body ? body : "")

    """
    cat > sample_name_key.tsv <<'EOF'
${text}
EOF
    """
}
