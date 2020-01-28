package command;

import de.embl.cba.morphometry.commands.BDOpenFolderCommand;
import de.embl.cba.morphometry.commands.BDOpenTableCommand;
import net.imagej.ImageJ;

public class RunBDOpenFolderCommand
{
	public static void main(final String... args)
	{
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		// invoke the plugin
		ij.command().run( BDOpenFolderCommand.class, true );
	}
}
