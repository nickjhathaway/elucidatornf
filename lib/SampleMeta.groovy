import java.nio.file.Files
import java.nio.file.Path

class SampleMeta {
  static final int MAX_LISTED = 20

/**
* The samples in the meta file (--meta_fnp), population clustering fails on a sample that isn't in it
*
* @return the names, or null when there's no meta (the empty file the param defaults to)
*/
  static Set<String> load(Path fnp) {
    if( !Files.isRegularFile(fnp) || Files.size(fnp) == 0 ) return null
    def rows = readRows(fnp)
    def pos = rows ? rows[0].indexOf('sample') : -1
    return pos < 0 ? new LinkedHashSet<String>() : rows.drop(1).findAll { row -> row.size() > pos }.collect { row -> row[pos] } as LinkedHashSet<String>
  }

/**
* What's wrong with the meta file itself: it has to be tab delimited with a header, a sample column and at
* least one more column, with a sample on only one row
*/
  static List<String> validate(Path fnp) {
    if( !Files.exists(fnp) ) return ["Path does not exist for --meta_fnp: '${fnp}'".toString()]
    if( !Files.isRegularFile(fnp) ) return ["Expected a file for --meta_fnp but found a directory: '${fnp}'".toString()]
    if( Files.size(fnp) == 0 ) return []
    def rows = readRows(fnp)
    def header = rows ? rows[0] : []
    def pos = header.indexOf('sample')
    if( pos < 0 ) return ["--meta_fnp needs a tab delimited header with a 'sample' column, found: ${header.join(', ')} ('${fnp}')".toString()]
    def problems = []
    if( header.size() < 2 ) problems << "--meta_fnp needs at least one column besides 'sample': '${fnp}'".toString()
    def dups = rows.drop(1).findAll { row -> row.size() > pos }.countBy { row -> row[pos] }.findAll { _name, count -> count > 1 }.keySet()
    if( dups ) problems << "--meta_fnp has samples on more than one row: ${listNames(dups)} ('${fnp}')".toString()
    return problems
  }

/**
* A message for the samples that aren't in the meta, null when they all are or there's no meta
*
* @param samples the final sample names, after the sample key
*/
  static String missingMessage(Set<String> metaSamples, Collection<String> samples) {
    if( metaSamples == null ) return null
    def missing = samples.findAll { name -> !metaSamples.contains(name) } as LinkedHashSet<String>
    if( !missing ) return null
    return "${missing.size()} of ${(samples as Set).size()} samples aren't in --meta_fnp: ${listNames(missing)}. " +
      "Add them to it, give them the names it uses with --amplicon_sample_key_fnp, or leave them out with --amplicon_samples_keep_file"
  }

  private static List<List<String>> readRows(Path fnp) {
    return Files.readAllLines(fnp).findAll { line -> line.trim() }.collect { line -> line.replaceAll(/\r$/, "").split('\t', -1) as List<String> }
  }

  private static String listNames(Collection<String> names) {
    def sorted = names.sort(false)
    return sorted.take(MAX_LISTED).join(', ') + (sorted.size() > MAX_LISTED ? ", and ${sorted.size() - MAX_LISTED} more" : "")
  }
}
