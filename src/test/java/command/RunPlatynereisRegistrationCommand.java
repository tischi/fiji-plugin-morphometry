package command;

import de.embl.cba.morphometry.commands.PlatynereisRegistrationCommand;
import net.imagej.ImageJ;

public class RunPlatynereisRegistrationCommand
{
	public static void main(final String... args)
	{
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		// invoke the plugin
		ij.command().run( PlatynereisRegistrationCommand.class, true );
	}
}
