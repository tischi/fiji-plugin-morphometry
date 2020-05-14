package users.tobias;

import de.embl.cba.morphometry.spindle.SpindleMorphometryCommand;
import loci.common.DebugTools;
import net.imagej.ImageJ;
import net.imglib2.type.numeric.RealType;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

public class TestSpindleAxialExtendUnderestimation
{
	public static < R extends RealType< R > > void main( String[] args )
	{
		DebugTools.setRootLevel("OFF"); // Bio-Formats

		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		final SpindleMorphometryCommand< R > command = new SpindleMorphometryCommand<>();
		command.opService = ij.op();

//		command.inputImageFile = new File( "/Users/tischer/Downloads/Incorrect-t82-crop.tif" );
//		command.spindleChannelIndexOneBased = 1;
//		command.dnaChannelIndexOneBased = 2;

		command.inputImageFile = new File( "/Users/tischer/Downloads/R1EWT_Undiffd0_aTub_568_CDK5RAP2_647_001-1.tif" );
		command.spindleChannelIndexOneBased = 2;
		command.dnaChannelIndexOneBased = 3;

		command.showIntermediateResults = true;
		command.saveResults = false;
		command.settings.showOutputImage = true;
		command.outputDirectory = new File("/Users/tischer/Downloads" );
		command.run();

		final HashMap< Integer, Map< String, Object > > measurements =
				command.getObjectMeasurements();
	}
}


