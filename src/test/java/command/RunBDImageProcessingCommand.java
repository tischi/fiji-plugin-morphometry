package command;

import de.embl.cba.morphometry.commands.BDImageProcessingCommand;
import de.embl.cba.morphometry.commands.DrosophilaRegistrationCommand;
import net.imagej.ImageJ;

public class RunBDImageProcessingCommand
{
	public static void main(final String... args)
	{
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		// invoke the plugin
		ij.command().run( BDImageProcessingCommand.class, true );
	}
}
