package users.valerie;

import de.embl.cba.morphometry.commands.MicrogliaMorphometryCommand;
import de.embl.cba.tables.FileUtils;
import de.embl.cba.tables.command.ConcatTablesCommand;
import de.embl.cba.tables.command.ExploreObjectsTableCommand;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import java.io.File;
import java.util.List;

public class MeasureMicrogliaMorphometryMultipleLargeData
{
	public static < T extends RealType< T > & NativeType< T > > void main( String[] args )
	{
		final net.imagej.ImageJ ij = new net.imagej.ImageJ();

		final File inputDir = new File( "/Users/tischer/Documents/valerie-blanche-petegnief-CORBEL-microglia-quantification--data/2019-February/several" );

		final List< File > files = FileUtils.getFileList( inputDir, ".*_20-labelMasks.tif", true );

		for ( File file : files )
		{
			final MicrogliaMorphometryCommand< T > measure = new MicrogliaMorphometryCommand<>();
			measure.labelMaskFile = file;
			measure.outputDirectory = inputDir;
			measure.opService = ij.op();
			measure.showIntermediateResults = false;
			measure.run();
		}


		final ConcatTablesCommand concatTablesCommand = new ConcatTablesCommand();
		concatTablesCommand.directory = inputDir;
		concatTablesCommand.regExp = ".*.tif.csv";
		concatTablesCommand.outputTable  = new File( inputDir + File.separator + "concat.csv" );
		concatTablesCommand.run();

		final ExploreObjectsTableCommand explore = new ExploreObjectsTableCommand();
		explore.tableFile = concatTablesCommand.outputTable;
		explore.is2D = true;
		explore.isPathMapping = false;
		explore.imageRootFolder = inputDir;
		explore.isOneBasedTimePoint = true;
		explore.imagePathColumnsId = "Path_";



	}
}
