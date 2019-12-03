package headless;

import de.embl.cba.morphometry.commands.BDImageProcessingCommand;
import de.embl.cba.morphometry.fccf.FCCF;
import net.imagej.ImageJ;

import java.io.File;

public class HeadlessBDImageProcessingCommand
{
	public static void main(final String... args)
	{
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		final BDImageProcessingCommand command = new BDImageProcessingCommand();
		command.logService = ij.log();
		command.viewingModality = FCCF.VIEW_PROCESSED_MONTAGE;
		command.inputDirectory = new File("/Users/tischer/Documents/BD-image-processing/sample_data/five_images");
		command.outputImageDirectory = new File( command.inputDirectory + File.separator + "output" );
		command.minBF = 0.0;
		command.maxBF = 1.0;
		command.minGFP = 0.0;
		command.maxGFP = 1.0;
		command.run();
	}
}
