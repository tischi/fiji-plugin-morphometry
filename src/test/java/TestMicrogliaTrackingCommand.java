import de.embl.cba.morphometry.microglia.MicrogliaTrackingCommand;
import net.imagej.ImageJ;

public class TestMicrogliaTrackingCommand
{

	public static void main(final String... args) throws Exception
	{
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		ij.command().run( MicrogliaTrackingCommand.class, true );
	}

}
