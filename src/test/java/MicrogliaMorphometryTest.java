import de.embl.cba.morphometry.Algorithms;
import de.embl.cba.morphometry.ImageIO;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.microglia.MicrogliaTrackingSettings;
import de.embl.cba.morphometry.measurements.ObjectMeasurements;
import de.embl.cba.morphometry.segmentation.SimpleSegmenter;
import de.embl.cba.morphometry.splitting.TrackingSplitter;
import de.embl.cba.morphometry.tracking.MaximalOverlapTracker;
import ij.IJ;
import ij.ImagePlus;
import net.imagej.ImageJ;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;

import java.io.File;
import java.util.ArrayList;
import java.util.Map;

public class MicrogliaMorphometryTest <T extends RealType< T > & NativeType< T > >
{

	public static final String SIMPLE_SEGMENTATION_TRACKING_SPLITTING_SIMPLE_TRACKING = "Simple segmentation - Tracking splitting - Simple tracking";
	public static final String INTENSITIES = "intensities";

	public void run()
	{
		ImageJ imagej = new ImageJ();
		imagej.ui().showUI();

//		final String path = MicrogliaMorphometryTest.class.getResource( "microglia/MAX_18C__t1-2.tif" ).getFile();
		//String path = "/Users/tischer/Documents/valerie-blanche-petegnief-CORBEL-microglia-quantification/data/MAX_18C.tif";

		String path = "/Users/tischer/Documents/valerie-blanche-petegnief-CORBEL-microglia-quantification--data/10C.lif";

		final ImagePlus imagePlus = ImageIO.openWithBioFormats( path );

		imagePlus.show();

		final Img< T > inputImages = ImageJFunctions.wrapReal( imagePlus );

		MicrogliaTrackingSettings settings = new MicrogliaTrackingSettings();
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
		settings.minimalObjectSize = 200;  // um2
		settings.minimalObjectCenterDistance = 6;
		settings.maximalWatershedLength = 10;
		settings.closingRadius = 3;
		settings.showIntermediateResults = false;
		settings.opService = imagej.op();
		settings.microgliaChannelIndexOneBased = 2;

		long tMin = inputImages.min( 4 );
		long tMax = 2; //inputImages.max( 2 );

		ArrayList< RandomAccessibleInterval< T > > intensities =
				Algorithms.createMaximumProjectedIntensitiesAssumingImagePlusDimensionOrder(
					inputImages,
					settings.microgliaChannelIndexOneBased - 1,
					tMin, tMax );

		ArrayList<  RandomAccessibleInterval< T > > masks = new ArrayList<>();

		for ( long t = tMin; t <= tMax; ++t )
		{
			final SimpleSegmenter simpleSegmenter = new SimpleSegmenter( intensities.get( ( int ) t ), settings );
			simpleSegmenter.run();
			masks.add( simpleSegmenter.getMask() );
		}

		ImageJFunctions.show( Views.stack( masks ), "Simple segmentation" );

		final TrackingSplitter splitter = new TrackingSplitter( masks, intensities, settings );
		splitter.run();
		masks = splitter.getSplitMasks();

		final MaximalOverlapTracker maximalOverlapTracker = new MaximalOverlapTracker( masks );
		maximalOverlapTracker.run();
		final ArrayList< RandomAccessibleInterval< T > > labelings = maximalOverlapTracker.getLabelings();

		Utils.createIJ1Movie( labelings, SIMPLE_SEGMENTATION_TRACKING_SPLITTING_SIMPLE_TRACKING );
		IJ.run("3-3-2 RGB", "");

		Utils.createIJ1Movie( intensities, INTENSITIES );
		IJ.wait( 500 );
		IJ.run("16-bit", "");

		IJ.wait( 500 );
		IJ.run("Merge Channels...",
				"c1=intensities c2=[Simple segmentation - Tracking splitting - Simple tracking] create");


//		ArrayList< RandomAccessibleInterval< T > > results = new ArrayList<>();
//		ArrayList< ImgLabeling< Integer, IntType > > imgLabelings = new ArrayList<>();




//
//		MicrogliaMorphometry morphometry = new MicrogliaMorphometry( settings, imagej.op() );
//		morphometry.run();
//		results.add( morphometry.getResultImageStack() );
//		measurements.add( morphometry.getObjectMeasurements() );
//		imgLabelings.add( morphometry.getImgLabeling() );
//		intensities.add( Views.hyperSlice( inputImages, 2, t) );
//
//
//
//		ArrayList< HashMap > measurements = new ArrayList<>();
//
//
//		RandomAccessibleInterval< T > movie = Views.addDimension( Views.stack( results ), 0, 0 );
//		movie = Views.permute( movie, 3, 4 );
//		ImagePlus show = ImageJFunctions.show( movie, settings.inputDataSetName );
//		IJ.run( show, "Grays", "");
//
//		final Overlay overlay = new Overlay();
//
//		Font font = new Font("SansSerif", Font.PLAIN, 10);
//
//		for ( long t = inputImages.min( 2 ); t <= tMax; ++t )
//		{
//			final HashMap hashMap = measurements.get( ( int ) t );
//			for ( Object label : hashMap.keySet() )
//			{
//				final Map< String, Object > objectMeasurements = ( Map< String, Object > ) hashMap.get( label );
//
//				long correctedSumIntensity = getCorrectedSumIntensity( objectMeasurements );
//
//				String text = "";
//				text += label;
//				text += ", size: " + objectMeasurements.get( ObjectMeasurements.SIZE_PIXEL_UNITS );
//				text += ", intens: " + correctedSumIntensity;
//				text += ", skel: " + Math.round( settings.workingVoxelSize * (long) objectMeasurements.get( ObjectMeasurements.SUM_INTENSITY + "_skeleton" ) );
//
//				double[] pixelPosition = getPixelPosition( settings, objectMeasurements );
//
//				final TextRoi textRoi = new TextRoi( pixelPosition[ 0 ], pixelPosition[ 1 ], text );
//				textRoi.setPosition( 1,1, (int) ( t+1 ) );
//				textRoi.setCurrentFont( font );
//				overlay.add( textRoi );
//
//			}
//		}
//
//		show.setOverlay( overlay );
//		final RandomAccessibleInterval< IntType > labelings = ( RandomAccessibleInterval ) Views.hyperSlice( movie, 2, 1 );
//		ImageJFunctions.show( labelings, "labelings" );
//
//		final RandomAccessibleInterval< LabelingType< Integer > > stack = Views.stack( imgLabelings );
//
//		final Tracker tracker = new Tracker( imgLabelings, intensities );
//		tracker.run();
//
//		ImagePlus relabelled =  ImageJFunctions.show(
//			Views.permute(
//					Views.addDimension(
//							Views.stack(  tracker.getLabelings() ), 0, 0 ), 2, 3)
//		);
//
//		IJ.saveAsTiff( relabelled, "/Users/tischer/Desktop/relabelled.tif" );


	}

	private ArrayList< RandomAccessibleInterval< T > > getIntensities( Img< T > inputImages, long channel, long tMin, long tMax )
	{
		ArrayList<  RandomAccessibleInterval< T > > intensities = new ArrayList<>();

		for ( long t = tMin; t <= tMax; ++t )
		{
			intensities.add(
					Views.hyperSlice( Views.hyperSlice( inputImages, 2, channel), 2, t ) );
		}

		return intensities;
	}


	public double[] getPixelPosition( MicrogliaTrackingSettings settings, Map< String, Object > objectMeasurements )
	{
		final double[] position = ( double[] ) objectMeasurements.get( ObjectMeasurements.CALIBRATED_POSITION );

		double[] pixelPosition = new double[ position.length ];
		for ( int d = 0; d < pixelPosition.length; ++d )
		{
			pixelPosition[ d ] = position[ d ] / settings.workingVoxelSize;
		}
		return pixelPosition;
	}

	public long getCorrectedSumIntensity( Map< String, Object > objectMeasurements )
	{
		final long sumIntensity = ( long ) objectMeasurements.get( ObjectMeasurements.SUM_INTENSITY + "_channel01" );
		final long size = ( long ) objectMeasurements.get( ObjectMeasurements.SIZE_PIXEL_UNITS );
		final double bgIntensity = ( double ) objectMeasurements.get( ObjectMeasurements.GOBAL_BACKGROUND_INTENSITY );
		return (long) ( sumIntensity - ( bgIntensity * size ) );
	}

	public static void main( String... args )
	{
		new MicrogliaMorphometryTest().run();

	}
}
