class SampleRename {
/**
* Load and validate a 2-column rename key: old_name -> new_name
*
* Rules:
*  - old names must be unique (no duplicate keys)
*  - new names must be unique (no collisions)
*  - new names must not contain '/', '\' or '?'
*
* Accepts TSV/CSV/whitespace-separated lines; ignores blank lines and lines starting with '#'.
*/
  static Map loadRenameKey(def renameKeyFnp) {
    if( renameKeyFnp == null ) return [:]

    def f = new File(renameKeyFnp.toString())
    if( !f.exists() )
      throw new IllegalArgumentException("rename key does not exist: ${renameKeyFnp}")

    def badCharsRe = /[\/\\\?]/
    def map = [:]
    def seenNew = new LinkedHashSet<String>()
    def dupOld = new LinkedHashSet<String>()
    def dupNew = new LinkedHashSet<String>()
    def badNew = new LinkedHashSet<String>()

    int lineNo = 0
    f.eachLine { raw ->
      lineNo++
      def line = raw?.trim()
      if( !line || line.startsWith('#') ) return

      def parts = line.split(/[,\t ]+/).findAll { it }
      if( parts.size() < 2 )
        throw new IllegalArgumentException("Invalid rename key line ${lineNo} in ${f}: '${raw}'")

      def oldName = parts[0].trim()
      def newName = parts[1].trim()

      if( map.containsKey(oldName) ) dupOld << oldName
      if( seenNew.contains(newName) ) dupNew << newName
      if( (newName =~ badCharsRe).find() ) badNew << newName

      map[oldName] = newName
      seenNew << newName
    }

    if( dupOld ) throw new IllegalArgumentException("Duplicate old sample names: ${dupOld.join(', ')}")
    if( dupNew ) throw new IllegalArgumentException("Duplicate new sample names: ${dupNew.join(', ')}")
    if( badNew ) throw new IllegalArgumentException("Invalid new sample names (/, \\, ?): ${badNew.join(', ')}")

    return map
  }
}
