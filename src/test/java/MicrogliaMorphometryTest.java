import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.microglia.MicrogliaMorphometry;
import de.embl.cba.morphometry.microglia.MicrogliaMorphometrySettings;
import de.embl.cba.morphometry.spindle.SpindleMorphometry;
import ij.IJ;
import ij.ImagePlus;
import net.imagej.ImageJ;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;

import java.io.File;

public class MicrogliaMorphometryTest<T extends RealType<T> & NativeType< T > >
{

	public void run()
	{
		ImageJ imagej = new ImageJ();
		imagej.ui().showUI();

		final String path = MicrogliaMorphometryTest.class.getResource( "microglia/MAX_18C__t1-2.tif" ).getFile();

		final ImagePlus imagePlus = IJ.openImage( path );
		final Img< T > img = ImageJFunctions.wrapReal( imagePlus );

		MicrogliaMorphometrySettings settings = new MicrogliaMorphometrySettings();
		settings.showIntermediateResults = true;
		settings.inputCalibration = Utils.get2dCalibration( imagePlus ) ;
		settings.image = Views.hyperSlice( img, 2, 0); // extract time point
		settings.workingVoxelSize = 0.5;
		settings.maxPossibleValueInDataSet = Math.pow( 2, imagePlus.getBitDepth() ) - 1.0;
		settings.maxShortAxisDist = 6;
		settings.thresholdInUnitsOfBackgroundPeakHalfWidth = 5.0;
		settings.watershedSeedsLocalMaximaDistanceThreshold = Double.MAX_VALUE;
		settings.watershedSeedsGlobalDistanceThreshold = 5.0;
		settings.interestPointsRadius = 0.5;
		settings.outputDirectory = new File( path ).getParentFile();
		settings.inputDataSetName = "test";
		settings.returnEarly = true;
		settings.minimalObjectSize = 75;

		MicrogliaMorphometry morphometry = new MicrogliaMorphometry( settings, imagej.op() );
		morphometry.run();

	}

	public static void main( String... args )
	{
		new MicrogliaMorphometryTest().run();

	}
}
