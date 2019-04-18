package command;

import de.embl.cba.morphometry.commands.PlatynereisRegistrationCommand;
import net.imagej.ImageJ;

import java.io.File;

public class RunHeadlessPlatynereisRegistrationCommand
{
	public static void main(final String... args)
	{

		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		final PlatynereisRegistrationCommand command = new PlatynereisRegistrationCommand();

		command.inputResolution = 1.0;
		command.inputDirectory = new File( "/Volumes/cba/exchange/Kimberly/Low_res/020_3_2__w7" );
		command.fileNameEndsWith = "png";
		command.outputDirectory = new File( "/Volumes/cba/exchange/Kimberly/Low_res/020_3_2__w7_aligned" );
		command.opService = ij.op();
		command.showIntermediateResults = false;
		command.outputResolution = 1.0;

		command.run();
	}
}
