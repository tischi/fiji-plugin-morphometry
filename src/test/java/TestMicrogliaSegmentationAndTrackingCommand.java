import de.embl.cba.morphometry.commands.MicrogliaSegmentationAndTrackingCommand;
import net.imagej.ImageJ;

public class TestMicrogliaSegmentationAndTrackingCommand
{

	public static void main(final String... args) throws Exception
	{
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		ij.command().run( MicrogliaSegmentationAndTrackingCommand.class, true );
	}

}
