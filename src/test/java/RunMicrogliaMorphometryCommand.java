import de.embl.cba.morphometry.commands.MicrogliaMorphometryCommand;
import net.imagej.ImageJ;

public class RunMicrogliaMorphometryCommand
{
	public static void main(final String... args)
	{
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		// invoke the plugin
		ij.command().run( MicrogliaMorphometryCommand.class, true );
	}
}
