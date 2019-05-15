package headless;

import de.embl.cba.morphometry.commands.PlatynereisRegistrationCommand;
import net.imagej.ImageJ;

import java.io.File;

public class PlatynereisRegistrationCommandKimberly
{
	public static void main(final String... args)
	{

		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		final PlatynereisRegistrationCommand command = new PlatynereisRegistrationCommand();

		// input
		command.inputResolutionMicrometer = 0.325;
		command.inputDirectory = new File( "/Volumes/cba/exchange/Kimberly/High_res/tomo_w3" );
		command.outputDirectory = new File( "/Volumes/cba/exchange/Kimberly/High_res/tomo_w3_aligned" );
		command.invertImage = true;
		command.fileNameEndsWith = "tiff";

//		command.inputResolutionMicrometer = 1.0;
//		command.inputDirectory = new File( "/Volumes/cba/exchange/Kimberly/Low_res/020_3_2__w7" );
//		command.outputDirectory = new File( "/Volumes/cba/exchange/Kimberly/Low_res/020_3_2__w7_aligned" );
//		command.invertImage = false;
//		command.fileNameEndsWith = "png";

		// output
		command.outputResolution = command.inputResolutionMicrometer;

		// other
		command.registrationResolution = 8.0;
		command.opService = ij.op();
		command.showIntermediateResults = false;

		command.run();
	}
}
