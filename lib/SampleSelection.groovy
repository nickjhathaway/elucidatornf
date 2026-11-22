class SampleSelection {
/**
* Load the samples to process, one name per line (only the first tab/comma/space delimited column is used);
* ignores blank lines and lines starting with '#'
*
* @return the names, or null when no file is given, meaning every sample is processed
*/
  static Set<String> load(def fnp) {
    if( fnp == null || fnp.toString().trim() == "" ) return null
    def f = new File(fnp.toString())
    if( !f.isFile() )
      throw new IllegalArgumentException("samples keep file does not exist: ${fnp}")
    def names = new LinkedHashSet<String>()
    f.eachLine { raw ->
      def line = raw.trim()
      if( !line || line.startsWith('#') ) return
      names << line.split(/[\t, ]/)[0].trim()
    }
    if( names.isEmpty() )
      throw new IllegalArgumentException("samples keep file has no sample names: ${fnp}")
    return names
  }

/**
* A sample key row is selected when its input, sample or replicate name is in keep
*/
  static boolean selected(Set<String> keep, Map row) {
    return keep == null || [row.input, row.sample, row.replicate].any { name -> keep.contains(name) }
  }

/**
* An input is processed when any of its rows in the sample key (or the input itself when it isn't in the key) is selected
*/
  static boolean selectsInput(Set<String> keep, Map<String, List<Map>> key, String input) {
    return keep == null || SampleKey.rowsForInput(key, input).any { row -> selected(keep, row) }
  }

/**
* The sample key with only the selected rows, inputs with none left are dropped
*/
  static Map<String, List<Map>> filterKey(Set<String> keep, Map<String, List<Map>> key) {
    if( keep == null ) return key
    return key.collectEntries { input, rows -> [(input): rows.findAll { row -> selected(keep, row) }] }
      .findAll { _input, rows -> !rows.isEmpty() }
  }

/**
* The names in keep that don't select any of the inputs found
*/
  static List<String> unmatched(Set<String> keep, Map<String, List<Map>> key, Collection<String> inputNames) {
    if( keep == null ) return []
    def found = inputNames.collectMany { input ->
      SampleKey.rowsForInput(key, input).collectMany { row -> [row.input, row.sample, row.replicate] }
    } as Set
    return keep.findAll { name -> !found.contains(name) }.sort()
  }
}
