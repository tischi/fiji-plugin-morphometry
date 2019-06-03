package users.tobias;

import de.embl.cba.morphometry.commands.SpindleMorphometryCommand;
import net.imagej.ImageJ;
import net.imglib2.type.numeric.RealType;

import java.io.File;

public class RunSpindleMorphometry
{
	public static < R extends RealType< R > > void main( String[] args )
	{
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		final SpindleMorphometryCommand< R > command = new SpindleMorphometryCommand<>();
		command.opService = ij.op();

		command.inputImageFile = new File("/Users/tischer/Desktop/kletter/HighZoom--W0000--P0001-T0004--0001.tif");

		//		command.inputImageFile = new File( "/Users/tischer/Desktop/kletter/HighZoom--W0000--P0001-T0007.tif" );



		command.outputDirectory = new File( "/Users/tischer/Desktop/kletter" );
		command.inputImageFilesParentDirectory = new File( "/Users/tischer/Desktop/kletter" );

		command.spindleChannelIndexOneBased = 1;
		command.dnaChannelIndexOneBased = 2;

		command.voxelSpacingDuringAnalysis = 0.25;

		command.showIntermediateResults = true;

		command.run();

	}
}