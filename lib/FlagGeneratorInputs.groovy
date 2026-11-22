import java.nio.file.Path

/**
* Manifest lines for AMPLICON_CLUSTER_AUTO_SEEKDEEP_FLAG_GENERATOR, which stages only the input directory and links the
* listed files into its own directory, so the number of staged inputs doesn't grow with the number of samples
*/
class FlagGeneratorInputs {

/**
* @param inputDir the input directory given to the process
* @param fnp a file under inputDir
* @param linkName the name to link it as, has to be unique among the lines
* @return the file's path relative to inputDir and the link name, tab separated
*/
  static String line(Path inputDir, Path fnp, String linkName) {
    def rel = inputDir.relativize(fnp).toString()
    if( rel.startsWith('..') )
      throw new IllegalArgumentException("${fnp} isn't under the input directory ${inputDir}")
    return "${rel}\t${linkName}".toString()
  }

/**
* The link name for a nanopore input, its sample name plus the file's extension (e.g. run1_barcode01.bam), as SeekDeep
* names nanopore samples by file name rather than by a bam's SM and files from different runs can share a name
*/
  static String nanoporeLinkName(String relPath, String sample) {
    def m = relPath =~ /(\.[^.\/]+\.gz|\.[^.\/]+)$/
    return sample + (m.find() ? m.group(1) : "")
  }
}
