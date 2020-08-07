package microglia;

import de.embl.cba.morphometry.commands.MicrogliaMorphometryCommand;
import de.embl.cba.tables.FileUtils;
import de.embl.cba.tables.command.ExploreObjectsTableCommand;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import java.io.File;
import java.util.List;

public class MeasureMicrogliaMorphometryOneImageManyCells
{
	public static < T extends RealType< T > & NativeType< T > > void main( String[] args )
	{
		final net.imagej.ImageJ ij = new net.imagej.ImageJ();

		final File inputDir = new File( "/Users/tischer/Documents/fiji-plugin-morphometry/src/test/resources/microglia/oneTimePointManyCells" );

		final List< File > files = FileUtils.getFileList( inputDir, ".*-labelMasks.tif", true );

		final MicrogliaMorphometryCommand< T > measure = new MicrogliaMorphometryCommand<>();
		measure.labelMaskFile = new File( inputDir + File.separator + "im-labelMasks.tif" );
		measure.outputDirectory = inputDir;
		measure.opService = ij.op();
		measure.showIntermediateResults = false;
		measure.run();

		final ExploreObjectsTableCommand explore = new ExploreObjectsTableCommand();
		explore.tableFile = measure.getTableOutputFile();
		explore.is2D = true;
		explore.isRelativeImagePath = true;
		explore.imageRootFolder = measure.outputDirectory;
		explore.isOneBasedTimePoint = true;
		explore.imagePathColumnsId = "Path_";

		explore.run();
	}
}
