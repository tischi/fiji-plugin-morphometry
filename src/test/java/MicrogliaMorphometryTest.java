import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.microglia.MicrogliaMorphometry;
import de.embl.cba.morphometry.microglia.MicrogliaMorphometrySettings;
import de.embl.cba.morphometry.objects.Measurements;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Overlay;
import ij.gui.TextRoi;
import net.imagej.ImageJ;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;

import java.awt.*;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class MicrogliaMorphometryTest <T extends RealType< T > & NativeType< T > >
{

	public void run()
	{
		ImageJ imagej = new ImageJ();
		imagej.ui().showUI();

//		final String path = MicrogliaMorphometryTest.class.getResource( "microglia/MAX_18C__t1-2.tif" ).getFile();
		String path = "/Users/tischer/Documents/valerie-blanche-petegnief-CORBEL-microglia-quantification/data/MAX_18C.tif";

		final ImagePlus imagePlus = IJ.openImage( path );
		final Img< T > img = ImageJFunctions.wrapReal( imagePlus );

		MicrogliaMorphometrySettings settings = new MicrogliaMorphometrySettings();
		settings.inputCalibration = Utils.get2dCalibration( imagePlus ) ;
		settings.workingVoxelSize = settings.inputCalibration[ 0 ];
		settings.maxPossibleValueInDataSet = Math.pow( 2, imagePlus.getBitDepth() ) - 1.0;
		settings.maxShortAxisDist = 6;
		settings.thresholdInUnitsOfBackgroundPeakHalfWidth = 5.0;
		settings.watershedSeedsLocalMaximaDistanceThreshold = Double.MAX_VALUE;
		settings.watershedSeedsGlobalDistanceThreshold = 2.5;
		settings.interestPointsRadius = 0.5;
		settings.outputDirectory = new File( path ).getParentFile();
		settings.inputDataSetName = "test";
		settings.returnEarly = true;
		settings.skeletonMaxLength = 600 * settings.workingVoxelSize;
		settings.minimalObjectSize = 500;
		settings.minimalObjectCenterDistance = 6;
		settings.maximalWatershedLength = 10;

		settings.closingRadius = 3;

		settings.showIntermediateResults = false;

		ArrayList< RandomAccessibleInterval< T > > timepoints = new ArrayList<>();
		ArrayList< HashMap > measurements = new ArrayList<>();

		long tMax = 0; //img.max( 2 );

		for ( long t = img.min( 2 ); t <= tMax; ++t )
		{
			Utils.log( "# Processing timepoint " + t );
			settings.image = Views.hyperSlice( img, 2, t); // extract time point
			MicrogliaMorphometry morphometry = new MicrogliaMorphometry( settings, imagej.op() );
			morphometry.run();
			timepoints.add( morphometry.getResultImageStack() );
			measurements.add( morphometry.getObjectMeasurements() );
		}

		RandomAccessibleInterval< T > movie = Views.addDimension( Views.stack( timepoints ), 0, 0 );
		movie = Views.permute( movie, 3, 4 );
		ImagePlus show = ImageJFunctions.show( movie, settings.inputDataSetName );
		IJ.run( show, "Grays", "");

		final Overlay overlay = new Overlay();

		Font font = new Font("SansSerif", Font.PLAIN, 10);

		for ( long t = img.min( 2 ); t <= tMax; ++t )
		{
			final HashMap hashMap = measurements.get( ( int ) t );
			for ( Object label : hashMap.keySet() )
			{
				final Map< String, Object > objectMeasurements = ( Map< String, Object > ) hashMap.get( label );

				long correctedSumIntensity = getCorrectedSumIntensity( objectMeasurements );

				String text = "";
				text += label;
				text += ", size: " + objectMeasurements.get( Measurements.SIZE_PIXEL_UNITS );
				text += ", intens: " + correctedSumIntensity;
				text += ", skel: " + Math.round( settings.workingVoxelSize * (long) objectMeasurements.get( Measurements.SUM_INTENSITY + "_skeleton" ) );

				double[] pixelPosition = getPixelPosition( settings, objectMeasurements );

				final TextRoi textRoi = new TextRoi( pixelPosition[ 0 ], pixelPosition[ 1 ], text );
				textRoi.setPosition( 1,1, (int) ( t+1 ) );
				textRoi.setCurrentFont( font );
				overlay.add( textRoi );

			}
		}

		show.setOverlay( overlay );
		final Overlay overlay1 = show.getOverlay();
		int a = 1;

	}

	public double[] getPixelPosition( MicrogliaMorphometrySettings settings, Map< String, Object > objectMeasurements )
	{
		final double[] position = ( double[] ) objectMeasurements.get( Measurements.CALIBRATED_POSITION );

		double[] pixelPosition = new double[ position.length ];
		for ( int d = 0; d < pixelPosition.length; ++d )
		{
			pixelPosition[ d ] = position[ d ] / settings.workingVoxelSize;
		}
		return pixelPosition;
	}

	public long getCorrectedSumIntensity( Map< String, Object > objectMeasurements )
	{
		final long sumIntensity = ( long ) objectMeasurements.get( Measurements.SUM_INTENSITY + "_channel01" );
		final long size = ( long ) objectMeasurements.get( Measurements.SIZE_PIXEL_UNITS );
		final double bgIntensity = ( double ) objectMeasurements.get( Measurements.GOBAL_BACKGROUND_INTENSITY );
		return (long) ( sumIntensity - ( bgIntensity * size ) );
	}

	public static void main( String... args )
	{
		new MicrogliaMorphometryTest().run();

	}
}
