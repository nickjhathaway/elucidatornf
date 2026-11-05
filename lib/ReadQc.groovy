class ReadQc {
/**
* What's wrong with the params for AMPLICON_READ_QC, shared by illumina_read_qc.nf and --illumina_clustering_read_qc
*
* @return a message for each problem found
*/
  static List<String> validate(def params, def genomeDir, Map<String, List<Map>> sampleKey) {
    def errors = []

    def primaryGenome = params.illumina_read_qc_primary_genome?.toString()?.trim()
    if( !primaryGenome ) {
      errors << "--illumina_read_qc_primary_genome must be supplied for the read QC"
    } else if( !new File("${genomeDir}/${primaryGenome}.fasta").exists() ) {
      errors << "Primary genome not found for --illumina_read_qc_primary_genome '${primaryGenome}': '${genomeDir}/${primaryGenome}.fasta' (expected in --genome_dir)"
    }

    def hostGenome = params.illumina_read_qc_host_genome_fnp
    if( hostGenome && !new File(hostGenome.toString()).isFile() ) {
      errors << "File does not exist for --illumina_read_qc_host_genome_fnp: '${hostGenome}'"
    }

    // the barcodes would stop the primers from being found so every pair would just be reported as failing primers
    def barcoded = sampleKey.findAll { _input, rows -> rows[0].fbarcode }.keySet().sort()
    if( barcoded ) {
      errors << "The read QC doesn't handle MIDs yet and these inputs have barcodes in the sample key: ${barcoded.join(", ")}"
    }
    return errors
  }

/**
* Put the fastq pairs from fromFilePairs, sorted by name so the task caches, into the input for AMPLICON_READ_QC
*
* @return a list of the names, a list of the first mates and a list of the second mates
*/
  static List samplesTuple(List fastqPairs) {
    return [
      fastqPairs.collect { pair -> pair[0] },
      fastqPairs.collect { pair -> pair[1][0] },
      fastqPairs.collect { pair -> pair[1][1] }
    ]
  }
}
