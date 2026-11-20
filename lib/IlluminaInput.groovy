import java.nio.file.FileVisitOption
import java.nio.file.Files
import java.nio.file.Path

class IlluminaInput {
  // the same as SeekDeep's illumina input pattern, the first group is the sample name, e.g. 9 for
  // 9_S390_R1_001.fastq.gz and Lib1 for Lib1_AAGAAGGC-AGCTACCA_S2_L001_R1_001.fastq.gz
  static final String illuminaNamePat = /^(.*?)(_[AGCT]{8,}-[AGCT]{8,})?(_S[0-9]+)?(_L[0-9]+)?_R[12](_[0-9]+)?\.fastq(\.gz)?$/

  // the first mates the pairGlobs find, the first group is the fromFilePairs prefix
  static final String firstMatePat = /^(.+)_R?1(_001)?\.fastq(\.gz)?$/

/**
* The fromFilePairs patterns for the paired fastqs in the input directory, ** also searches its subdirectories
*
* @return the patterns for the _1/_2 and _R1/_R2 namings
*/
  static List<String> pairGlobs(def inputDir, boolean recursive) {
    def star = recursive ? "**" : "*"
    return [
      "${inputDir}/${star}_{1,2}{_001,}.fastq{.gz,}".toString(),
      "${inputDir}/${star}_{R1,R2}{_001,}.fastq{.gz,}".toString()
    ]
  }

/**
* The sample name for a fastq pair, with strip the illumina suffix (index barcodes, _S, _L, _R1/_R2, _001) is
* removed like SeekDeep does, unless the sample key has the full prefix as an input
*
* @param pairPrefix the name fromFilePairs gave the pair
* @param firstMate the first mate's fastq
* @return the sample name
*/
  static String sampleName(String pairPrefix, def firstMate, boolean strip, Set<String> keyInputs) {
    if( !strip || keyInputs.contains(pairPrefix) ) {
      return pairPrefix
    }
    def m = firstMate.getFileName().toString() =~ illuminaNamePat
    return (m.matches() && m.group(1)) ? m.group(1) : pairPrefix
  }

/**
* What's wrong with the fastq pairs in the input directory, more than one pair getting the same sample name
* (the same name in several subdirectories, or several lanes once the suffix is stripped)
*
* @return a message for each sample name given to more than one pair
*/
  static List<String> validate(Path inputDir, boolean recursive, boolean strip, Set<String> keyInputs) {
    if( !Files.isDirectory(inputDir) ) {
      return []
    }
    def pairsBySample = [:].withDefault { [] }
    def stream = Files.walk(inputDir, recursive ? Integer.MAX_VALUE : 1, FileVisitOption.FOLLOW_LINKS)
    try {
      stream.filter { Path p -> Files.isRegularFile(p) }.each { Path p ->
        def m = p.fileName.toString() =~ firstMatePat
        if( m.matches() ) {
          pairsBySample[sampleName(m.group(1), p, strip, keyInputs)] << inputDir.relativize(p).toString()
        }
      }
    } finally {
      stream.close()
    }
    def hint = strip ? ", lanes aren't combined and --illumina_strip_sample_suffix false keeps the full names" : ""
    return pairsBySample.findAll { _sample, fastqs -> fastqs.size() > 1 }.collect { sample, fastqs ->
      "Sample name '${sample}' would be given to more than one fastq pair in --input_fastq_dir, rename them so they're unique${hint}: ${fastqs.sort().join(", ")}".toString()
    }
  }
}
