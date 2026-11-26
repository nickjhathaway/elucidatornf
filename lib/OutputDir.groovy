import java.nio.file.FileVisitResult
import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.SimpleFileVisitor
import java.nio.file.attribute.BasicFileAttributes

class OutputDir {
/**
* Delete dir, except for keep (and the directories leading to it) when keep is inside dir, so the
* downloaded genomes in <outdir>/plasmodium_genomes survive the output directory being reset
*
* @param keep a path to leave alone, or null to delete everything
* @param runTimestamp params.trace_timestamp, files with it in their name are this run's and are left alone
*   (nextflow already has the trace file open)
* @param protect paths that mustn't be deleted, e.g. [workflow.launchDir, workflow.workDir]; it's an error for dir to be
*   one of them or to contain one, rather than wiping the run's own .nextflow cache, log and work directory
*/
  static void clear(def dir, def keep, def runTimestamp, List protect = []) {
    def dirPath = Paths.get(dir.toString()).toAbsolutePath().normalize()
    if( !Files.exists(dirPath, LinkOption.NOFOLLOW_LINKS) ) return
    def protectedInside = protect.findAll { p -> p != null }
      .collect { p -> Paths.get(p.toString()).toAbsolutePath().normalize() }
      .findAll { p -> p.startsWith(dirPath) }
    if( protectedInside ) {
      throw new IllegalArgumentException("Won't reset the output directory ${dirPath}, it is or contains ${protectedInside.join(", ")} (the launch or work directory), which would delete this run's own .nextflow cache, log or work files; give an --outdir that is separate from them")
    }
    def keepPath = keep == null ? null : Paths.get(keep.toString()).toAbsolutePath().normalize()
    if( keepPath == null || !keepPath.startsWith(dirPath) ) {
      deleteTree(dirPath, runTimestamp)
      return
    }
    if( keepPath == dirPath ) return
    Files.list(dirPath).withCloseable { children ->
      children.toList().each { Path child ->
        if( keepPath.startsWith(child) ) {
          clear(child, keepPath, runTimestamp)
        } else {
          deleteTree(child, runTimestamp)
        }
      }
    }
  }

  // symlinks are removed, not followed; directories are left when all that's in them was left
  private static void deleteTree(Path path, def runTimestamp) {
    def left = new HashSet<Path>()
    Files.walkFileTree(path, new SimpleFileVisitor<Path>() {
      @Override
      FileVisitResult visitFile(Path file, BasicFileAttributes attrs) {
        if( isLeft(file, runTimestamp) ) {
          left << file
        } else {
          // deleteIfExists, a .nfs placeholder can go away once its last handle closes mid walk
          Files.deleteIfExists(file)
        }
        return FileVisitResult.CONTINUE
      }

      // a file that went away between being listed and being visited (e.g. a .nfs placeholder) is fine
      @Override
      FileVisitResult visitFileFailed(Path file, IOException exc) {
        if( exc instanceof NoSuchFileException ) return FileVisitResult.CONTINUE
        throw exc
      }

      @Override
      FileVisitResult postVisitDirectory(Path d, IOException exc) {
        if( exc != null ) throw exc
        def remaining = Files.list(d).withCloseable { children -> children.toList() }
        if( remaining && remaining.every { Path child -> left.contains(child) || isNfsPlaceholder(child) } ) {
          left << d
        } else {
          Files.deleteIfExists(d)
        }
        return FileVisitResult.CONTINUE
      }
    })
  }

  static boolean isLeft(Path file, def runTimestamp) {
    return isNfsPlaceholder(file) || (runTimestamp && file.fileName.toString().contains(runTimestamp.toString()))
  }

  // NFS renames a deleted file that's still open to .nfs<id>, it goes away once the file is closed
  static boolean isNfsPlaceholder(Path file) {
    return file.fileName.toString().startsWith(".nfs")
  }
}
