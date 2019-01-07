import de.embl.cba.morphometry.Rois;
import de.embl.cba.morphometry.spindle.SpindleMorphometryCommand;
import de.embl.cba.morphometry.translocation.TranslocationCommand;
import ij.IJ;
import ij.ImagePlus;
import ij.plugin.frame.RoiManager;
import net.imagej.ImageJ;
import net.imglib2.FinalInterval;

import java.util.ArrayList;

public class RunTranslocationCommand
{
	public static void main( String[] args )
	{
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		// stage image
		//
		IJ.open( TranslocationExample01.class.getResource(
						"translocation/test01.zip" ).getFile() );

		// stage roimanager
		//
		final RoiManager rm = new RoiManager();
		rm.runCommand( "open", TranslocationExample01.class.getResource( "translocation/test01-rois.zip" ).getFile() );

		// invoke the plugin
		ij.command().run( TranslocationCommand.class, true );
	}
}
