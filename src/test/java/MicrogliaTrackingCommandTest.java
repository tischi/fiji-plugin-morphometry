import de.embl.cba.morphometry.drosophila.shavenbaby.ShavenBabyRegistrationCommand;
import de.embl.cba.morphometry.microglia.MicrogliaTrackingCommand;
import net.imagej.ImageJ;

public class MicrogliaTrackingCommandTest
{

	public static void main(final String... args) throws Exception
	{
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		// invoke the plugin
		ij.command().run( MicrogliaTrackingCommand.class, true );
	}

}
