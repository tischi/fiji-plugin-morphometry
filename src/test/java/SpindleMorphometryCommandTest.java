import de.embl.cba.morphometry.microglia.MicrogliaTrackingCommand;
import de.embl.cba.morphometry.spindle.SpindleMorphometryCommand;
import net.imagej.ImageJ;

public class SpindleMorphometryCommandTest
{
	public static void main( String[] args )
	{
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		// invoke the plugin
		ij.command().run( SpindleMorphometryCommand.class, true );
	}
}
