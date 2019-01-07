import de.embl.cba.morphometry.Rois;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.translocation.TranslocationComputer;
import de.embl.cba.morphometry.translocation.TranslocationResult;
import ij.IJ;
import ij.ImagePlus;
import ij.plugin.frame.RoiManager;
import net.imagej.ImageJ;
import net.imglib2.FinalInterval;

import java.util.ArrayList;

public class ExampleTranslocation01
{
	public static void main( String[] args )
	{
		final ImageJ ij = new ImageJ();

		final ImagePlus imagePlus = IJ.openImage( ExampleTranslocation01.class.getResource( "translocation/test01.zip" ).getFile() );

		final RoiManager rm = new RoiManager();
		rm.runCommand( "open", ExampleTranslocation01.class.getResource( "translocation/test01-rois.zip" ).getFile() );
		final ArrayList< FinalInterval > intervals = Rois.asIntervals( rm.getRoisAsArray() );

		final TranslocationComputer computer = new TranslocationComputer(
				Utils.get2DImagePlusMovieAsFrameList( imagePlus, 1 ),
				intervals,
				ij.op() );

		final ArrayList< TranslocationResult > results = computer.getResults();

		for ( int i = 0; i < results.size(); i++ )
		{
			Utils.listOf2DImagesAsImagePlusMovie(
					results.get( 0 ).cellMasks,
					"cell masks 0" ).show();
		}

	}

}
