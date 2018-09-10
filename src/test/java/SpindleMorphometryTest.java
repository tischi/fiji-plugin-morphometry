import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.spindle.SpindleMorphometrySettings;
import ij.IJ;
import ij.ImagePlus;
import net.imagej.ImageJ;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import de.embl.cba.morphometry.spindle.*;

import java.io.File;

public class SpindleMorphometryTest <T extends RealType<T> & NativeType< T > >
{

	public void run()
	{
		ImageJ imagej = new ImageJ();
		imagej.ui().showUI();

		final String dapiPath = SpindleMorphometryTest.class.getResource( "spindle/test002-C01.tif" ).getFile();
		final String tubulinPath = SpindleMorphometryTest.class.getResource( "spindle/test002-C00.tif" ).getFile();

		final ImagePlus dapiImp = IJ.openImage( dapiPath );
		final ImagePlus tubulinImp = IJ.openImage( tubulinPath );

		final double[] calibration = Utils.getCalibration( dapiImp );

		final Img< T > dapi = ImageJFunctions.wrapReal( dapiImp );
		final Img< T > tubulin = ImageJFunctions.wrapReal( tubulinImp );

		SpindleMorphometrySettings settings = new SpindleMorphometrySettings();
		settings.showIntermediateResults = true;
		settings.inputCalibration = Utils.getCalibration( dapiImp );
		settings.dapi = dapi;
		settings.tubulin = tubulin;
		settings.workingVoxelSize = 0.25;
		settings.maxPossibleValueInDataSet = Math.pow( 2, dapiImp.getBitDepth() ) - 1.0;
		settings.maxShortAxisDist = 6;
		settings.thresholdInUnitsOfBackgroundPeakHalfWidth = 5.0;
		settings.watershedSeedsLocalMaximaDistanceThreshold = 1.0;
		settings.watershedSeedsGlobalDistanceThreshold = 2.0;
		settings.interestPointsRadius = 0.5;
		settings.outputDirectory = new File( dapiPath ).getParentFile();
		settings.inputDataSetName = "test";

		SpindleMorphometry morphometry = new SpindleMorphometry( settings, imagej.op() );
		morphometry.run();

	}

	public static void main( String... args )
	{
		new SpindleMorphometryTest().run();

	}
}
