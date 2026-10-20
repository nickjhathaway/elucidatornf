import java.nio.file.Files
import java.nio.file.Path
import java.util.zip.GZIPInputStream

class RnaseqRefFiles {
/**
* Reference inputs of rnaseq_process_pf.nf, resolved from --rnaseq_ref_files_dir and/or
* their individual flags, then sanity checked before anything runs.
*
* Each input comes from its own flag when that is set, otherwise from the expected name
* under --rnaseq_ref_files_dir, so one file can be swapped without changing the directory.
*/

  // param name -> expected name under --rnaseq_ref_files_dir
  static final Map<String, String> EXPECTED = [
    rnaseq_host1_filter_genome_fasta_fnp : 'host1_genome.fasta',
    rnaseq_host1_filter_contigs_fnp      : 'host1_remove_contigs.txt',
    rnaseq_host2_filter_genome_fasta_fnp : 'host2_genome.fasta',
    rnaseq_host2_filter_contigs_fnp      : 'host2_remove_contigs.txt',
    rnaseq_salmon_index                  : 'salmon_index',
    single_cell_seurat_r_object_fnp      : 'single_cell_seurat.rds',
    rnaseq_kmer_filter_table_fnp         : 'kmer_filter_table.tsv.gz',
    rnaseq_kmer_sort_table_fnp           : 'kmer_sort_table.tsv.gz',
  ]

  static final List<String> BWA_INDEX_EXTENSIONS = ['amb', 'ann', 'bwt', 'pac', 'sa']

  /**
  * elucidator extractByCountingUniqKmersFromSets uses nonUniqueKmers_<table name> from the
  * table's directory when it exists; returns that file, or null when there is none.
  */
  static Path nonUniqueKmerTable(Path kmerTable) {
    if( kmerTable == null ) return null
    def candidate = kmerTable.resolveSibling("nonUniqueKmers_${kmerTable.fileName}")
    return Files.isRegularFile(candidate) ? candidate : null
  }

  /**
  * Resolve every reference input. Returns [files: param -> Path, sources: param -> description,
  * errors: [...]]; `toPath` turns a string into a Path (the workflow passes Nextflow's file()).
  */
  static Map resolve(Map params, Closure toPath) {
    def files = [:]
    def sources = [:]
    def errors = []
    def refDir = params.rnaseq_ref_files_dir ? toPath(params.rnaseq_ref_files_dir.toString()) : null
    if( refDir != null && !Files.isDirectory(refDir) )
      errors << "--rnaseq_ref_files_dir is not a directory: ${refDir}"

    EXPECTED.each { param, expectedName ->
      def flagValue = params[param]
      if( flagValue ) {
        files[param] = toPath(flagValue.toString())
        sources[param] = "--${param}"
      }
      else if( refDir != null ) {
        files[param] = refDir.resolve(expectedName)
        sources[param] = "--rnaseq_ref_files_dir (${expectedName})"
      }
      else {
        errors << "no ${expectedName} given: set --rnaseq_ref_files_dir or --${param}"
      }
    }
    return [files: files, sources: sources, errors: errors]
  }

  /** One line per input with where it came from, plus the non-unique k-mer tables picked up. */
  static String describe(Map files, Map sources) {
    def lines = EXPECTED.keySet().collect { param -> "  ${param}: ${files[param]}  [${sources[param]}]".toString() }
    ['rnaseq_kmer_filter_table_fnp', 'rnaseq_kmer_sort_table_fnp'].each { param ->
      def nonUnique = nonUniqueKmerTable(files[param] as Path)
      lines << "  non-unique k-mers for ${param}: ${nonUnique ?: 'none found'}".toString()
    }
    return lines.join('\n')
  }

  /** Check every resolved file exists and looks like what it should be; returns problems found. */
  static List<String> validate(Map files, Map sources) {
    def errors = []
    def check = { String param, Closure test ->
      def f = files[param] as Path
      if( f == null ) return
      try {
        test(f).each { problem -> errors << "${problem}\n      file: ${f}\n      from: ${sources[param]}".toString() }
      }
      catch( Exception e ) {
        errors << "could not read ${f} (from ${sources[param]}): ${e.message}".toString()
      }
    }

    ['host1', 'host2'].each { host ->
      check("rnaseq_${host}_filter_genome_fasta_fnp") { Path f -> checkGenomeFasta(f, host) }
      check("rnaseq_${host}_filter_contigs_fnp") { Path f ->
        def genome = files["rnaseq_${host}_filter_genome_fasta_fnp"] as Path
        checkContigsInGenome(f, genome, host)
      }
    }
    check('rnaseq_salmon_index') { Path f -> checkSalmonIndex(f) }
    check('rnaseq_kmer_filter_table_fnp') { Path f ->
      def nonUnique = nonUniqueKmerTable(f)
      checkKmerTable(f, 'k-mer filter table') + (nonUnique ? checkKmerTable(nonUnique, "its non-unique k-mer table ${nonUnique.fileName}") : [])
    }
    check('rnaseq_kmer_sort_table_fnp') { Path f ->
      def nonUnique = nonUniqueKmerTable(f)
      checkKmerTable(f, 'k-mer sort table') + (nonUnique ? checkKmerTable(nonUnique, "its non-unique k-mer table ${nonUnique.fileName}") : [])
    }
    check('single_cell_seurat_r_object_fnp') { Path f -> checkRds(f) }
    return errors
  }

  static List<String> checkGenomeFasta(Path f, String host) {
    if( !Files.isRegularFile(f) ) return ["${host} genome FASTA does not exist"]
    def problems = []
    if( Files.size(f) == 0 ) problems << "${host} genome FASTA is empty"
    else if( firstNonSpaceChar(openMaybeGz(f)) != '>' ) problems << "${host} genome is not a FASTA file (does not start with '>')"
    def missingBwa = BWA_INDEX_EXTENSIONS.findAll { ext ->
      def idx = f.resolveSibling("${f.fileName}.${ext}")
      !Files.isRegularFile(idx) || Files.size(idx) == 0
    }
    if( missingBwa ) problems << "${host} genome is missing its bwa index: ${missingBwa.collect { ext -> "${f.fileName}.${ext}" }.join(', ')} (make with: bwa index ${f.fileName})"
    def mmi = mmiPath(f)
    if( !Files.isRegularFile(mmi) ) problems << "${host} genome is missing its minimap2 index ${mmi.fileName} (make with: minimap2 -d ${mmi.fileName} ${f.fileName})"
    else if( !startsWith(mmi, [0x4d, 0x4d, 0x49, 0x02] as byte[]) ) problems << "${mmi.fileName} is not a minimap2 index (no MMI header)"
    return problems
  }

  // matches how the subworkflow derives the .mmi path from the genome FASTA
  static Path mmiPath(Path fasta) {
    return fasta.resolveSibling(fasta.fileName.toString().replaceFirst(/\.(fa|fasta|fna)(\.gz)?$/, '.mmi'))
  }

  /** Every contig to remove must be a contig of that host's genome, read from the bwa .ann. */
  static List<String> checkContigsInGenome(Path f, Path genome, String host) {
    if( !Files.isRegularFile(f) ) return ["${host} contigs-to-remove file does not exist"]
    def wanted = new LinkedHashSet<String>()
    f.eachLine { line ->
      def name = line.trim()
      if( name && !name.startsWith('#') ) wanted << name.split(/\s+/)[0]
    }
    if( wanted.isEmpty() ) return ["${host} contigs-to-remove file lists no contigs"]
    def ann = genome?.resolveSibling("${genome.fileName}.ann")
    if( ann == null || !Files.isRegularFile(ann) ) return []  // reported with the genome
    // .ann: a header line, then two lines per contig: "<gi> <name> <description>" and "<offset> <len> <n_ambs>"
    def contigs = new HashSet<String>()
    int lineNo = 0
    ann.eachLine { line ->
      lineNo++
      if( lineNo > 1 && lineNo % 2 == 0 ) {
        def parts = line.split(' ', 3)
        if( parts.size() >= 2 ) contigs << parts[1]
      }
    }
    def missing = wanted.findAll { name -> !contigs.contains(name) }
    if( missing ) {
      def shown = missing.take(5).join(', ') + (missing.size() > 5 ? ", ... (${missing.size()} in all)" : '')
      return ["${missing.size()} of ${wanted.size()} ${host} contigs to remove are not in the ${host} genome (${genome.fileName}): ${shown}".toString()]
    }
    return []
  }

  static List<String> checkSalmonIndex(Path d) {
    if( !Files.isDirectory(d) ) return ['salmon index directory does not exist']
    def info = d.resolve('info.json')
    if( Files.isRegularFile(info) ) {
      def version = (info.text =~ /"salmon_version"\s*:\s*"([^"]+)"/)
      if( version.find() && !version.group(1).startsWith('2') )
        return ["salmon index was built with salmon ${version.group(1)}; the pipeline image runs salmon 2.x".toString()]
      return []
    }
    def oldInfo = d.resolve('versionInfo.json')
    if( Files.isRegularFile(oldInfo) ) {
      def version = (oldInfo.text =~ /"salmonVersion"\s*:\s*"([^"]+)"/)
      def built = version.find() ? version.group(1) : '1.x'
      return ["salmon index was built with salmon ${built}, which salmon 2.x cannot read; rebuild it with salmon 2.x (salmon index -t <transcripts> -i <index>)".toString()]
    }
    return ['not a salmon index (no info.json)']
  }

  /** k-mer set table: gzipped, lines "<set><TAB><kmer>". */
  static List<String> checkKmerTable(Path f, String label) {
    if( !Files.isRegularFile(f) ) return ["${label} does not exist"]
    if( !startsWith(f, [0x1f, 0x8b] as byte[]) ) return ["${label} is not gzip compressed (expected <name>.tsv.gz)"]
    def first = null
    new BufferedReader(new InputStreamReader(new GZIPInputStream(Files.newInputStream(f)))).withCloseable { r -> first = r.readLine() }
    def parts = first?.split('\t')
    if( parts == null || parts.size() != 2 || !(parts[1] ==~ /[ACGTacgt]+/) )
      return ["${label} does not look like a k-mer set table (expected '<set><TAB><kmer>' lines, first line: '${first?.take(80)}')".toString()]
    return []
  }

  /** R serialization: gzip/bzip2/xz compressed or plain, starting with an X/A/B format byte and a newline. */
  static List<String> checkRds(Path f) {
    if( !Files.isRegularFile(f) ) return ['single-cell RDS file does not exist']
    byte[] head = new byte[6]
    def n = 0
    Files.newInputStream(f).withCloseable { s -> n = s.read(head) }
    def isGzip = n >= 2 && (head[0] & 0xff) == 0x1f && (head[1] & 0xff) == 0x8b
    def isBzip2 = n >= 3 && new String(head, 0, 3) == 'BZh'
    def isXz = n >= 6 && (head[0] & 0xff) == 0xfd && new String(head, 1, 4) == '7zXZ'
    if( isBzip2 || isXz ) return []
    byte[] start = head
    if( isGzip ) {
      start = new byte[2]
      new GZIPInputStream(Files.newInputStream(f)).withCloseable { s -> s.read(start) }
    }
    // compare as Strings: a Groovy Character never equals the String 'X'
    def format = new String(start, 0, 2, 'ISO-8859-1')
    if( format in ['X\n', 'A\n', 'B\n'] ) return []
    return ['single-cell file is not an R .rds file (no R serialization header)']
  }

  static boolean startsWith(Path f, byte[] magic) {
    byte[] head = new byte[magic.length]
    def n = 0
    Files.newInputStream(f).withCloseable { s -> n = s.read(head) }
    return n == magic.length && Arrays.equals(head, magic)
  }

  static InputStream openMaybeGz(Path f) {
    return startsWith(f, [0x1f, 0x8b] as byte[]) ? new GZIPInputStream(Files.newInputStream(f)) : Files.newInputStream(f)
  }

  static String firstNonSpaceChar(InputStream s) {
    return s.withCloseable { stream ->
      int c
      while( (c = stream.read()) != -1 ) {
        if( !Character.isWhitespace(c) ) return String.valueOf((char) c)
      }
      return null
    }
  }
}
