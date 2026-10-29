class SampleKey {
/**
* Load and validate an amplicon sample key, one row per input (and barcode) with a header:
*
*   input     sample   replicate   fbarcode   rbarcode
*
*  - input: the input name as scanned (fastq pair prefix for illumina, e.g. Lib1_S1 for Lib1_S1_R1_001.fastq.gz)
*  - sample: the final sample name, a sample on more than one row has replicates that are
*    combined during population clustering
*  - replicate (optional): name for the replicate, when blank it's the sample name if the sample
*    only has one row, otherwise the input name (plus the barcodes for barcoded rows)
*  - fbarcode/rbarcode (optional): the MID(s) of this sample within the input, rbarcode only for dual barcoding
*
* A headerless 2 column file (old_name new_name) is still read as a plain rename key.
*
* Rules:
*  - an input + barcodes combination can only appear once
*  - replicate names must be unique (they name the per replicate files)
*  - within an input either all rows have barcodes or none do, with only one row when none do
*  - names must not contain '/', '\', '?' or whitespace
*
* Tab or comma delimited, blank optional columns are kept; ignores blank lines and lines starting with '#'.
*
* @return Map of input name -> list of rows, each a map with keys input, sample, replicate, fbarcode, rbarcode
*/
  static Map<String, List<Map>> load(def keyFnp) {
    if( keyFnp == null ) return [:]

    def f = new File(keyFnp.toString())
    if( !f.exists() )
      throw new IllegalArgumentException("sample key does not exist: ${keyFnp}")

    def lines = []
    int lineNo = 0
    f.eachLine { raw ->
      lineNo++
      def line = raw?.trim()
      if( !line || line.startsWith('#') ) return
      lines << [lineNo, raw]
    }
    if( lines.isEmpty() ) return [:]

    def delim = lines[0][1].contains('\t') ? '\t' : ','
    def header = lines[0][1].split(delim, -1).collect { it -> it.trim().toLowerCase() }
    if( !header.contains('input') || !header.contains('sample') ) {
      // headerless 2 column rename key, the old format
      def renameKey = SampleRename.loadRenameKey(keyFnp)
      return renameKey.collectEntries { oldName, newName ->
        [(oldName): [[input: oldName, sample: newName, replicate: newName, fbarcode: "", rbarcode: ""]]]
      }
    }
    def unknown = header.findAll { col -> !(col in ['input', 'sample', 'replicate', 'fbarcode', 'rbarcode']) }
    if( unknown )
      throw new IllegalArgumentException("Unrecognized columns in sample key ${f}: ${unknown.join(', ')}, allowed columns are input, sample, replicate, fbarcode, rbarcode")
    if( header.contains('rbarcode') && !header.contains('fbarcode') )
      throw new IllegalArgumentException("Sample key ${f} has an rbarcode column without an fbarcode column")

    def rows = []
    lines.drop(1).each { entry ->
      def (no, raw) = entry
      def toks = raw.split(delim, -1).collect { it -> it.trim() }
      if( toks.size() > header.size() && toks.drop(header.size()).any { it -> it } )
        throw new IllegalArgumentException("Sample key line ${no} in ${f} has more columns than the header: '${raw}'")
      def get = { col -> def pos = header.indexOf(col); (pos >= 0 && pos < toks.size()) ? toks[pos] : "" }
      def row = [input: get('input'), sample: get('sample'), replicate: get('replicate'),
                 fbarcode: get('fbarcode').toUpperCase(), rbarcode: get('rbarcode').toUpperCase(), line: no]
      if( !row.input || !row.sample )
        throw new IllegalArgumentException("Sample key line ${no} in ${f} is missing an input or sample: '${raw}'")
      if( row.rbarcode && !row.fbarcode )
        throw new IllegalArgumentException("Sample key line ${no} in ${f} has an rbarcode without an fbarcode: '${raw}'")
      rows << row
    }

    def errors = []
    def badCharsRe = /[\/\\\?\s]/
    def barRe = /^[ACGTN]*$/

    def rowsBySample = rows.groupBy { it -> it.sample }
    rows.each { row ->
      if( !row.replicate ) {
        if( rowsBySample[row.sample].size() == 1 ) {
          row.replicate = row.sample
        } else {
          row.replicate = [row.input, row.fbarcode, row.rbarcode].findAll { it -> it }.join('-')
        }
      }
      ['input', 'sample', 'replicate'].each { col ->
        if( (row[col] =~ badCharsRe).find() ) errors << "line ${row.line}: ${col} '${row[col]}' can't contain /, \\, ? or whitespace"
      }
      ['fbarcode', 'rbarcode'].each { col ->
        if( !(row[col] ==~ barRe) ) errors << "line ${row.line}: ${col} '${row[col]}' should only contain A, C, G, T or N"
      }
    }

    rows.groupBy { it -> [it.input, it.fbarcode, it.rbarcode] }.each { key, dups ->
      if( dups.size() > 1 ) errors << "input ${key[0]} with barcodes '${key[1]}' '${key[2]}' appears more than once (lines ${dups*.line.join(', ')})"
    }
    rows.groupBy { it -> it.replicate }.each { rep, dups ->
      if( dups.size() > 1 ) errors << "replicate name ${rep} is used more than once (lines ${dups*.line.join(', ')})"
    }
    // a replicate named after another sample would read as that sample in per replicate outputs
    def sampleNames = rowsBySample.keySet()
    rows.each { row ->
      if( row.replicate != row.sample && row.replicate in sampleNames )
        errors << "line ${row.line}: replicate name ${row.replicate} is also the name of a sample"
    }
    def rowsByInput = rows.groupBy { it -> it.input }
    rowsByInput.each { input, inputRows ->
      def barcoded = inputRows.findAll { it -> it.fbarcode }
      if( barcoded && barcoded.size() != inputRows.size() ) {
        errors << "input ${input} has rows with and without barcodes, either all rows for an input have barcodes or none do"
      } else if( !barcoded && inputRows.size() > 1 ) {
        errors << "input ${input} has ${inputRows.size()} rows but no barcodes to separate them"
      }
      def dual = barcoded.findAll { it -> it.rbarcode }
      if( dual && dual.size() != barcoded.size() ) {
        errors << "input ${input} mixes single and dual barcoded rows"
      }
    }

    if( errors ) throw new IllegalArgumentException("Errors in sample key ${f}:\n  ${errors.join('\n  ')}")

    return rowsByInput.collectEntries { input, inputRows ->
      [(input): inputRows.collect { it -> it.findAll { k, _v -> k != 'line' } }]
    }
  }

/**
* The sample key file from the params, --amplicon_sample_key_fnp or the older --nanopore_rename_key_fnp
*/
  static def keyFnp(def params) {
    if( params.amplicon_sample_key_fnp && params.nanopore_rename_key_fnp )
      throw new IllegalArgumentException("Only one of --amplicon_sample_key_fnp and --nanopore_rename_key_fnp can be given")
    return params.amplicon_sample_key_fnp ?: params.nanopore_rename_key_fnp
  }

/**
* The rows for an input, an input that isn't in the key is a single sample named after the input;
* each row gets extraction_name, the name the extractor uses in its output files for that row
*/
  static List<Map> rowsForInput(Map<String, List<Map>> key, String input) {
    def rows = key.containsKey(input) ? key[input] : [[input: input, sample: input, replicate: input, fbarcode: "", rbarcode: ""]]
    return rows.collect { row -> row + [extraction_name: row.fbarcode ? row.replicate : input] }
  }

/**
* Names that inputs left out of the key can't have, as they would be read as a keyed sample or replicate
*/
  static Set<String> keyedNames(Map<String, List<Map>> key) {
    return key.values().flatten().collectMany { row -> [row.sample, row.replicate] } as Set
  }

  static final String keyHeader = "input\tsample\treplicate\tfbarcode\trbarcode\textraction_name"

  static String keyLine(Map row) {
    return [row.input, row.sample, row.replicate, row.fbarcode, row.rbarcode, row.extraction_name].join('\t')
  }

/**
* Write the primers plus the MIDs for one barcoded input as a SeekDeep id file, the MIDs are named by
* their replicate so the extractor writes <target><replicate>.fastq.gz
*/
  static void writeIdFile(def primersFnp, List<Map> inputRows, def outFnp) {
    def out = new File(outFnp.toString())
    out.parentFile?.mkdirs()
    def primerLines = new File(primersFnp.toString()).readLines().findAll { it -> it.trim() }
    def dual = inputRows.any { it -> it.rbarcode }
    out.withWriter { w ->
      primerLines.each { it -> w.println(it) }
      w.println(dual ? "id\tbarcode\trevBarcode" : "id\tbarcode")
      inputRows.each { row ->
        w.println(dual ? "${row.replicate}\t${row.fbarcode}\t${row.rbarcode}" : "${row.replicate}\t${row.fbarcode}")
      }
    }
  }

/**
* Write the primers plus every distinct barcode combination in the key as a SeekDeep id file, used to
* guess the extractor flags across all inputs; returns null when the key has no barcodes
*/
  static def writeCombinedIdFile(def primersFnp, Map<String, List<Map>> key, def outFnp) {
    def barcodes = key.values().flatten().findAll { it -> it.fbarcode }.collect { it -> [it.fbarcode, it.rbarcode] }.unique()
    if( barcodes.isEmpty() ) return null
    def dual = barcodes.any { it -> it[1] }
    def width = barcodes.size().toString().length()
    def rows = barcodes.withIndex().collect { bars, idx ->
      [replicate: "MID${(idx + 1).toString().padLeft(width, '0')}".toString(), fbarcode: bars[0], rbarcode: dual ? bars[1] : ""]
    }
    writeIdFile(primersFnp, rows, outFnp)
    return outFnp
  }
}
