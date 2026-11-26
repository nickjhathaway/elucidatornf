import java.nio.file.FileVisitOption
import java.nio.file.Files
import java.nio.file.Path

class IlluminaInput {
  // the same as SeekDeep's illumina input pattern, the first group is the sample name, e.g. 9 for
  // 9_S390_R1_001.fastq.gz and Lib1 for Lib1_AAGAAGGC-AGCTACCA_S2_L001_R1_001.fastq.gz
  static final String illuminaNamePat = /^(.*?)(_[AGCT]{8,}-[AGCT]{8,})?(_S[0-9]+)?(_L[0-9]+)?_R[12](_[0-9]+)?\.fastq(\.gz)?$/

  // the first mates the pairGlobs find, the first group is the fromFilePairs prefix
  static final String firstMatePat = /^(.+)_R?1(_001)?\.fastq(\.gz)?$/

  // single end fastqs (--illumina_single_end), the first group is the file name minus the extension
  static final String singleEndPat = /^(.+?)\.fastq(\.gz)?$/

  // the read length genTargetInfoFromGenomes gets for single end input, single end reads cover the whole target so
  // it's longer than any amplicon, which makes every target overlap and keeps its refSeqs whole
  static final int singleEndReadLength = 10000

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
* The sample names of the fastq pairs in the input directory, like the pairGlobs channel a first mate without its
* second mate is left out
*/
  static Set<String> sampleNames(Path inputDir, boolean recursive, boolean strip, Set<String> keyInputs) {
    def names = new TreeSet<String>()
    eachFirstMate(inputDir, recursive) { Path p, String pairPrefix ->
      def mateName = p.fileName.toString().replaceFirst(/1((_001)?\.fastq(\.gz)?)$/, '2$1')
      if( Files.exists(p.resolveSibling(mateName)) ) names << sampleName(pairPrefix, p, strip, keyInputs)
    }
    return names
  }

/**
* The fromPath pattern for the single end fastqs in the input directory (--illumina_single_end), ** also searches its
* subdirectories
*/
  static String singleEndGlob(def inputDir, boolean recursive) {
    def star = recursive ? "**" : "*"
    return "${inputDir}/${star}.fastq{.gz,}".toString()
  }

/**
* The sample name for a single end fastq, the file name minus its extension, with strip an illumina named file
* (9_S390_L001_R1_001.fastq.gz) loses its suffix like a pair would, unless the sample key has the full name as an input
*
* @param fastq the fastq
* @return the sample name
*/
  static String singleEndSampleName(def fastq, boolean strip, Set<String> keyInputs) {
    def fileName = fastq.getFileName().toString()
    def m = fileName =~ singleEndPat
    def base = m.matches() ? m.group(1) : fileName
    if( !strip || keyInputs.contains(base) ) {
      return base
    }
    def illuminaMatch = fileName =~ illuminaNamePat
    return (illuminaMatch.matches() && illuminaMatch.group(1)) ? illuminaMatch.group(1) : base
  }

/**
* The sample names of the single end fastqs in the input directory
*/
  static Set<String> singleEndSampleNames(Path inputDir, boolean recursive, boolean strip, Set<String> keyInputs) {
    def names = new TreeSet<String>()
    eachSingleEnd(inputDir, recursive) { Path p -> names << singleEndSampleName(p, strip, keyInputs) }
    return names
  }

/**
* Like validate for the single end fastqs, more than one fastq getting the same sample name (the same name in several
* subdirectories, or the R1 and R2 of a pair, which aren't single end)
*/
  static List<String> validateSingleEnd(Path inputDir, boolean recursive, boolean strip, Set<String> keyInputs, Closure<Boolean> keepSample = { _name -> true }) {
    def fastqsBySample = [:].withDefault { [] }
    eachSingleEnd(inputDir, recursive) { Path p ->
      def sample = singleEndSampleName(p, strip, keyInputs)
      if( keepSample(sample) ) {
        fastqsBySample[sample] << inputDir.relativize(p).toString()
      }
    }
    def hint = strip ? ", --illumina_strip_sample_suffix false keeps the full names" : ""
    return fastqsBySample.findAll { _sample, fastqs -> fastqs.size() > 1 }.collect { sample, fastqs ->
      "Sample name '${sample}' would be given to more than one fastq in --input_fastq_dir with --illumina_single_end, rename them so they're unique (paired fastqs need the default paired end mode)${hint}: ${fastqs.sort().join(", ")}".toString()
    }
  }

  private static void eachSingleEnd(Path inputDir, boolean recursive, Closure action) {
    if( !Files.isDirectory(inputDir) ) return
    def stream = Files.walk(inputDir, recursive ? Integer.MAX_VALUE : 1, FileVisitOption.FOLLOW_LINKS)
    try {
      stream.filter { Path p -> Files.isRegularFile(p) && p.fileName.toString() ==~ singleEndPat }.each { Path p -> action(p) }
    } finally {
      stream.close()
    }
  }

  private static void eachFirstMate(Path inputDir, boolean recursive, Closure action) {
    if( !Files.isDirectory(inputDir) ) return
    def stream = Files.walk(inputDir, recursive ? Integer.MAX_VALUE : 1, FileVisitOption.FOLLOW_LINKS)
    try {
      stream.filter { Path p -> Files.isRegularFile(p) }.each { Path p ->
        def m = p.fileName.toString() =~ firstMatePat
        if( m.matches() ) action(p, m.group(1))
      }
    } finally {
      stream.close()
    }
  }

/**
* What's wrong with the fastq pairs in the input directory, more than one pair getting the same sample name
* (the same name in several subdirectories, or several lanes once the suffix is stripped)
*
* @param keepSample whether a sample name is processed, collisions between samples that aren't are ignored
* @return a message for each sample name given to more than one pair
*/
  static List<String> validate(Path inputDir, boolean recursive, boolean strip, Set<String> keyInputs, Closure<Boolean> keepSample = { _name -> true }) {
    def pairsBySample = [:].withDefault { [] }
    eachFirstMate(inputDir, recursive) { Path p, String pairPrefix ->
      def sample = sampleName(pairPrefix, p, strip, keyInputs)
      if( keepSample(sample) ) {
        pairsBySample[sample] << inputDir.relativize(p).toString()
      }
    }
    def hint = strip ? ", lanes aren't combined and --illumina_strip_sample_suffix false keeps the full names" : ""
    return pairsBySample.findAll { _sample, fastqs -> fastqs.size() > 1 }.collect { sample, fastqs ->
      "Sample name '${sample}' would be given to more than one fastq pair in --input_fastq_dir, rename them so they're unique${hint}: ${fastqs.sort().join(", ")}".toString()
    }
  }
}
