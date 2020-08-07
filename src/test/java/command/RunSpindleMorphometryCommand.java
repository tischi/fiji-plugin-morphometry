package command;

import de.embl.cba.morphometry.commands.Spindle3DCommand;
import net.imagej.ImageJ;

public class RunSpindleMorphometryCommand
{
	public static void main( String[] args )
	{
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		// invoke the plugin
		ij.command().run( Spindle3DCommand.class, true );
	}
}
