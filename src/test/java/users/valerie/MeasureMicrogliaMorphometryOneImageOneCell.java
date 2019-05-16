package users.valerie;

import de.embl.cba.morphometry.commands.MicrogliaMorphometryCommand;
import de.embl.cba.tables.command.ExploreObjectsTableCommand;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import java.io.File;

public class MeasureMicrogliaMorphometryOneImageOneCell
{

	public static < T extends RealType< T > & NativeType< T > > void main( String[] args )
	{
		final net.imagej.ImageJ ij = new net.imagej.ImageJ();

		final MicrogliaMorphometryCommand< T > measure = new MicrogliaMorphometryCommand<>();
		measure.labelMaskFile = new File( "/Users/tischer/Documents/fiji-plugin-morphometry/src/test/resources/microglia/MAX_5C-crop-t1-3-labelMasks.tif" );
		measure.outputDirectory = measure.labelMaskFile.getParentFile();
		measure.opService = ij.op();
		measure.showIntermediateResults = false;

		measure.run();


		final ExploreObjectsTableCommand explore = new ExploreObjectsTableCommand();
		explore.tableFile = measure.getTableOutputFile();
		explore.is2D = true;
		explore.isPathMapping = false;
		explore.isRelativeImagePath = true;
		explore.imageRootFolder = measure.outputDirectory;
		explore.isOneBasedTimePoint = true;
		explore.imagePathColumnsId = "Path_";

		explore.run();
	}
}
