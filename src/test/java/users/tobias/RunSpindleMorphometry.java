package users.tobias;

import de.embl.cba.morphometry.commands.Spindle3DCommand;
import net.imagej.ImageJ;
import net.imglib2.type.numeric.RealType;

import java.io.File;

public class RunSpindleMorphometry
{
	public static < R extends RealType< R > > void main( String[] args )
	{
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		final Spindle3DCommand< R > command = new Spindle3DCommand<>();
		command.opService = ij.op();

		command.inputImageFile = new File("/Users/tischer/Downloads/20190827_T0248_A-8degrees.tif");
//		command.inputImageFile = new File("/Users/tischer/Downloads/20190827_T0145_C-28degrees.tif");
//		command.inputImageFile = new File( "/Users/tischer/Downloads/20190827_T0075_B-8degrees.tif" );

		command.outputDirectory = new File( "/Users/tischer/Desktop/kletter" );
		command.inputImageFilesParentDirectory = new File( "/Users/tischer/Desktop/kletter" );

		command.spindleChannelIndexOneBased = 1;
		command.dnaChannelIndexOneBased = 2;

		command.voxelSpacingDuringAnalysis = 0.25; // normally 0.25

		command.showIntermediateResults = false;

		command.run();
	}
}
