import de.embl.cba.morphometry.commands.TranslocationCommand;
import ij.IJ;
import ij.plugin.frame.RoiManager;
import net.imagej.ImageJ;

public class TranslocationCommandRun
{
	public static void main( String[] args )
	{
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		IJ.open( TranslocationCommandTest.class.getResource(
						"translocation/test01.zip" ).getFile() );

		final RoiManager rm = new RoiManager();
		rm.runCommand( "open", TranslocationCommandTest.class.getResource( "translocation/test01-rois.zip" ).getFile() );

		// invoke the plugin
		ij.command().run( TranslocationCommand.class, true );
	}
}