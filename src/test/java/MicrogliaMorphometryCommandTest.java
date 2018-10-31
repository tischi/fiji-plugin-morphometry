import de.embl.cba.morphometry.microglia.MicrogliaMorphometryCommand;
import de.embl.cba.morphometry.microglia.MicrogliaTrackingCommand;
import net.imagej.ImageJ;

public class MicrogliaMorphometryCommandTest
{

	public static void main(final String... args) throws Exception
	{
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		// invoke the plugin
		ij.command().run( MicrogliaMorphometryCommand.class, true );
	}

}
