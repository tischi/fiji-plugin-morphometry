package tests;

import de.embl.cba.morphometry.commands.TranslocationCommand;
import ij.IJ;
import ij.plugin.frame.RoiManager;
import net.imagej.ImageJ;

import java.io.File;

public class TestTranslocationCommand
{
	public static void main( String[] args )
	{
		final ImageJ ij = new ImageJ();
		ij.ui().showUI();

		IJ.open( "translocation/test01.zip" );

		final RoiManager rm = new RoiManager();
		rm.runCommand( "open", "translocation/test01-rois.zip" );

		final TranslocationCommand command = new TranslocationCommand();

		command.reviewMembraneSegmentations = true;
		command.intensitiesImp = IJ.getImage();
		command.showTranslocationPlots = true;
		command.opService = ij.op();
		command.outputDirectory = new File( "translocation/output" );
		command.showIntermediateResults = false;

		command.run();




		// TODO: Add a test to check whether it works
	}
}
