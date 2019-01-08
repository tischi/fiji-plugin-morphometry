import de.embl.cba.morphometry.commands.TranslocationCommand;
import ij.IJ;
import ij.plugin.frame.RoiManager;
import net.imagej.ImageJ;

public class RunTranslocationCommand
{
	public static void main( String[] args )
	{
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		// stage image
		//
		IJ.open( TranslocationExample.class.getResource(
						"translocation/test01.zip" ).getFile() );

		// stage roimanager
		//
		final RoiManager rm = new RoiManager();
		rm.runCommand( "open", TranslocationExample.class.getResource( "translocation/test01-rois.zip" ).getFile() );

		// invoke the plugin
		ij.command().run( TranslocationCommand.class, true );
	}
}
