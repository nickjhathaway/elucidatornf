import java.net.URI
import java.net.http.HttpClient
import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.StandardCopyOption
import java.security.MessageDigest
import java.time.Duration
import java.util.zip.GZIPInputStream
import java.util.zip.ZipInputStream

/**
* Download and check the reference genomes from https://github.com/nickjhathaway/plasmodium_genomes for
* --download_plasmodium_genomes. The repo layout is kept, so the genomes are in <root>/genomes, the gffs in
* <root>/info/gff and the drug resistance positions in <root>/info, which is where the workflows look for them.
*
* The download runs on the head job (not in a process) so the files exist when the workflows check their
* inputs; the indexing needs the container so that's INDEX_PLASMODIUM_GENOMES
*/
class PlasmodiumGenomes {
  static final String REPO = "nickjhathaway/plasmodium_genomes"
  static final String DEFAULT_DIR_NAME = "plasmodium_genomes"
  static final String DEFAULT_GENOME = "Pf3D7"
  static final String LFS_POINTER_START = "version https://git-lfs.github.com/spec/v1"

  // what elucidator bioIndexGenomes adds to each <genome> name (2bit, bowtie2, bwa, samtools faidx, blast)
  static final List<String> INDEX_SUFFIXES = [".2bit", ".1.bt2", ".fasta.bwt", ".fasta.fai", "_blastdb.ndb"]

  static boolean enabled(def params) {
    return params.download_plasmodium_genomes != null && params.download_plasmodium_genomes.toString().toBoolean()
  }

/**
* --plasmodium_genomes_dir, or <outdir>/plasmodium_genomes when it isn't given
*/
  static Path rootDir(def params) {
    def dir = params.plasmodium_genomes_dir?.toString()?.trim()
    if( !dir ) {
      def outdir = params.outdir?.toString()?.trim()
      if( !outdir )
        throw new IllegalArgumentException("--outdir or --plasmodium_genomes_dir must be given with --download_plasmodium_genomes")
      dir = "${outdir}/${DEFAULT_DIR_NAME}"
    }
    return Paths.get(dir).toAbsolutePath().normalize()
  }

  static Path genomeDir(Path root) { root.resolve("genomes") }

  static Path gffDir(Path root) { root.resolve("info").resolve("gff") }

  static List<String> genomeNames(Path root) {
    def dir = genomeDir(root).toFile()
    if( !dir.isDirectory() ) return []
    return dir.listFiles().findAll { f -> f.isFile() && f.name.endsWith(".fasta") }
      .collect { f -> f.name.substring(0, f.name.length() - ".fasta".length()) }
      .sort()
  }

/**
* Every genome is uncompressed and has a gff, and every .gz has its gunzipped file next to it
*/
  static List<String> extractionProblems(Path root) {
    def names = genomeNames(root)
    if( names.isEmpty() ) return ["no .fasta files in ${genomeDir(root)}".toString()]
    def problems = []
    names.each { g ->
      if( !Files.exists(gffDir(root).resolve("${g}.gff")) ) problems << "missing ${gffDir(root).resolve("${g}.gff")}".toString()
    }
    if( Files.exists(root.resolve(".download_tmp")) ) problems << "unfinished download in ${root.resolve(".download_tmp")}".toString()
    root.toFile().eachFileRecurse { f ->
      if( f.isFile() && f.name.endsWith(".gz") && !new File(f.parentFile, f.name.replaceAll(/\.gz$/, "")).exists() )
        problems << "not gunzipped ${f}".toString()
    }
    return problems
  }

  static List<String> missingIndexFiles(Path root) {
    return genomeNames(root).collectMany { g ->
      INDEX_SUFFIXES.collect { suffix -> genomeDir(root).resolve("${g}${suffix}".toString()) }.findAll { p -> !Files.exists(p) }
    }.collect { p -> p.toString() }
  }

  static boolean isExtracted(Path root) { extractionProblems(root).isEmpty() }

  static boolean isIndexed(Path root) { isExtracted(root) && missingIndexFiles(root).isEmpty() }

/**
* Download <ref> of the repo into root and gunzip it unless root is already extracted. Files already there are
* left alone (every file is written to a temp file then moved) so an interrupted download picks up where it stopped
*/
  static void download(Path root, String ref, def log) {
    if( isExtracted(root) ) return
    Files.createDirectories(root)
    def client = HttpClient.newBuilder()
      .followRedirects(HttpClient.Redirect.NORMAL)
      .connectTimeout(Duration.ofSeconds(60))
      .build()
    def tmpDir = root.resolve(".download_tmp")
    Files.createDirectories(tmpDir)
    def archiveUrl = "https://github.com/${REPO}/archive/${ref}.zip"
    log.info "Downloading ${REPO} ${ref} into ${root}"
    def zipFnp = tmpDir.resolve("archive.zip")
    fetch(client, archiveUrl, zipFnp)

    new ZipInputStream(Files.newInputStream(zipFnp)).withCloseable { zip ->
      def entry
      while( (entry = zip.getNextEntry()) != null ) {
        // drop the <repo>-<ref>/ top directory the archive adds
        def relative = entry.name.contains("/") ? entry.name.substring(entry.name.indexOf("/") + 1) : ""
        if( entry.isDirectory() || !relative ) continue
        def target = root.resolve(relative).normalize()
        if( !target.startsWith(root) ) throw new IllegalStateException("unexpected path in ${archiveUrl}: ${entry.name}")
        def finalTarget = relative.endsWith(".gz") ? Paths.get(target.toString().replaceAll(/\.gz$/, "")) : target
        if( Files.exists(finalTarget) ) continue
        Files.createDirectories(target.getParent())

        def content = zip.readAllBytes()
        def tmp = tmpDir.resolve(relative.replace("/", "__"))
        def text = new String(content, 0, Math.min(content.length, LFS_POINTER_START.length()), "UTF-8")
        if( text == LFS_POINTER_START ) {
          def pointer = parseLfsPointer(new String(content, "UTF-8"))
          log.info "Downloading ${relative} (${String.format('%.1f', pointer.bytes / 1e6)} MB)"
          fetch(client, "https://media.githubusercontent.com/media/${REPO}/${ref}/${relative}", tmp)
          checkLfsObject(tmp, pointer, relative)
        } else {
          Files.write(tmp, content)
        }
        if( relative.endsWith(".gz") ) {
          def unzipped = tmpDir.resolve(tmp.fileName.toString().replaceAll(/\.gz$/, ""))
          new GZIPInputStream(Files.newInputStream(tmp)).withCloseable { input ->
            Files.copy(input, unzipped, StandardCopyOption.REPLACE_EXISTING)
          }
          Files.delete(tmp)
          tmp = unzipped
        }
        Files.move(tmp, finalTarget, StandardCopyOption.REPLACE_EXISTING)
      }
    }
    tmpDir.toFile().deleteDir()

    def problems = extractionProblems(root)
    if( problems )
      throw new IllegalStateException("${REPO} ${ref} didn't download as expected into ${root}:\n  ${problems.join('\n  ')}")
  }

  private static void fetch(HttpClient client, String url, Path out) {
    def request = HttpRequest.newBuilder(URI.create(url)).timeout(Duration.ofMinutes(30)).GET().build()
    def response = client.send(request, HttpResponse.BodyHandlers.ofFile(out))
    if( response.statusCode() != 200 ) {
      Files.deleteIfExists(out)
      throw new IllegalStateException("failed to download ${url}: HTTP ${response.statusCode()}")
    }
  }

  private static Map parseLfsPointer(String text) {
    def fields = text.readLines().findAll { line -> line.contains(" ") }.collectEntries { line ->
      def idx = line.indexOf(" ")
      [(line.substring(0, idx)): line.substring(idx + 1).trim()]
    }
    return [sha256: fields['oid'].toString().replaceFirst(/^sha256:/, ""), bytes: fields['size'] as long]
  }

  private static void checkLfsObject(Path fnp, Map pointer, String relative) {
    if( Files.size(fnp) != pointer.bytes )
      throw new IllegalStateException("downloaded ${relative} is ${Files.size(fnp)} bytes, expected ${pointer.bytes}")
    def digest = MessageDigest.getInstance("SHA-256")
    Files.newInputStream(fnp).withCloseable { input ->
      def buf = new byte[1 << 16]
      int n
      while( (n = input.read(buf)) > 0 ) digest.update(buf, 0, n)
    }
    def sha = digest.digest().encodeHex().toString()
    if( sha != pointer.sha256 )
      throw new IllegalStateException("downloaded ${relative} has sha256 ${sha}, expected ${pointer.sha256}")
  }
}
