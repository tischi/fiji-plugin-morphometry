import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.microglia.MicrogliaMorphometry;
import de.embl.cba.morphometry.microglia.MicrogliaMorphometrySettings;
import de.embl.cba.morphometry.spindle.SpindleMorphometry;
import ij.IJ;
import ij.ImagePlus;
import net.imagej.ImageJ;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

import java.io.File;
import java.util.ArrayList;

public class MicrogliaMorphometryTest <T extends RealType< T > & NativeType< T > >
{

	public void run()
	{
		ImageJ imagej = new ImageJ();
		imagej.ui().showUI();

//		final String path = MicrogliaMorphometryTest.class.getResource( "microglia/MAX_18C__t1-2.tif" ).getFile();
		String path = "/Users/tischer/Documents/valerie-blanche-petegnief-CORBEL-microglia-quantification/data/MAX_18C__t1-25.tif";

		final ImagePlus imagePlus = IJ.openImage( path );
		final Img< T > img = ImageJFunctions.wrapReal( imagePlus );

		MicrogliaMorphometrySettings settings = new MicrogliaMorphometrySettings();
		settings.inputCalibration = Utils.get2dCalibration( imagePlus ) ;
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
		settings.minimalObjectSize = 30;

		settings.showIntermediateResults = false;

		ArrayList< RandomAccessibleInterval< T > > timepoints = new ArrayList<>();
		for ( long t = img.min( 2 ); t <= img.max( 2 ); ++t )
		{
			Utils.log( "# Processing timepoint " + t );
			settings.image = Views.hyperSlice( img, 2, t); // extract time point
			MicrogliaMorphometry morphometry = new MicrogliaMorphometry( settings, imagej.op() );
			timepoints.add( morphometry.run() );
		}

		RandomAccessibleInterval< T > movie = Views.addDimension( Views.stack( timepoints ), 0, 0 );
		movie = Views.permute( movie, 3, 4 );
		ImagePlus show = ImageJFunctions.show( movie, settings.inputDataSetName );
		IJ.saveAsTiff( show, path + "--output.tif" );


	}

	public static void main( String... args )
	{
		new MicrogliaMorphometryTest().run();

	}
}
