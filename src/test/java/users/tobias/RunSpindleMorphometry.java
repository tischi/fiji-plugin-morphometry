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

		command.inputImageFile = new File("/Users/tischer/Downloads/Composite_calibrated.tif");

		command.outputDirectory = new File( "/Users/tischer/Desktop/kletter" );
		command.inputImageFilesParentDirectory = new File( "/Users/tischer/Desktop/kletter" );

		command.spindleChannelIndexOneBased = 1;
		command.dnaChannelIndexOneBased = 2;

		command.voxelSpacingDuringAnalysis = 0.24; // normally 0.25

		command.showIntermediateResults = false;

		command.run();
	}
}
