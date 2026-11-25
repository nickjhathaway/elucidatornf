import java.nio.file.FileVisitResult
import java.nio.file.Files
import java.nio.file.LinkOption
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
*/
  static void clear(def dir, def keep) {
    def dirPath = Paths.get(dir.toString()).toAbsolutePath().normalize()
    if( !Files.exists(dirPath, LinkOption.NOFOLLOW_LINKS) ) return
    def keepPath = keep == null ? null : Paths.get(keep.toString()).toAbsolutePath().normalize()
    if( keepPath == null || !keepPath.startsWith(dirPath) ) {
      deleteTree(dirPath)
      return
    }
    if( keepPath == dirPath ) return
    Files.list(dirPath).withCloseable { children ->
      children.toList().each { Path child ->
        if( keepPath.startsWith(child) ) {
          clear(child, keepPath)
        } else {
          deleteTree(child)
        }
      }
    }
  }

  // symlinks are removed, not followed
  private static void deleteTree(Path path) {
    Files.walkFileTree(path, new SimpleFileVisitor<Path>() {
      @Override
      FileVisitResult visitFile(Path file, BasicFileAttributes attrs) {
        Files.delete(file)
        return FileVisitResult.CONTINUE
      }

      @Override
      FileVisitResult postVisitDirectory(Path d, IOException exc) {
        if( exc != null ) throw exc
        Files.delete(d)
        return FileVisitResult.CONTINUE
      }
    })
  }
}
