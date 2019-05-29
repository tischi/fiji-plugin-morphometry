package command;

import de.embl.cba.morphometry.commands.SpindleMorphometryCommand;
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

		IJ.open( "/Users/tischer/Documents/fiji-plugin-morphometry/src/test/resources/translocation/test01.zip" );

		final RoiManager rm = new RoiManager();
		rm.runCommand( "open", "/Users/tischer/Documents/fiji-plugin-morphometry/src/test/resources/translocation/test01-rois.zip" );

		// invoke the plugin
		ij.command().run( TranslocationCommand.class, true );
	}
}
