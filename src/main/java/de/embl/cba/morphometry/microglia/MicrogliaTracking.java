package de.embl.cba.morphometry.microglia;

import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.segmentation.SimpleSegmenter;
import de.embl.cba.morphometry.splitting.TrackingSplitter;
import de.embl.cba.morphometry.tracking.MaximalOverlapTracker;
import ij.IJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

import java.io.File;
import java.util.ArrayList;

import static de.embl.cba.morphometry.microglia.Constants.INTENSITIES;
import static de.embl.cba.morphometry.microglia.Constants.SEGMENTATION;

public class MicrogliaTracking < T extends RealType< T > & NativeType< T > >
{
	final MicrogliaTrackingSettings settings;
	final ArrayList< RandomAccessibleInterval< T  > > intensities;
	private ArrayList< RandomAccessibleInterval< T > > labelings;

	public MicrogliaTracking ( ImagePlus imagePlus,
							   boolean showIntermediateResults,
							   OpService opService,
							   long microgliaChannelIndexOneBased,
							   long tMin,
							   long tMax )
	{

		settings = new MicrogliaTrackingSettings();

		settings.showIntermediateResults = showIntermediateResults;
		settings.opService = opService;
		settings.microgliaChannelIndexOneBased = microgliaChannelIndexOneBased;
		settings.tMin = tMin;
		settings.tMax = Math.min( tMax, imagePlus.getNFrames() );
		settings.inputCalibration = Utils.get2dCalibration( imagePlus ) ;
		settings.workingVoxelSize = settings.inputCalibration[ 0 ];
		settings.maxPossibleValueInDataSet = Math.pow( 2, imagePlus.getBitDepth() ) - 1.0;
		settings.maxShortAxisDist = 6;
		settings.thresholdInUnitsOfBackgroundPeakHalfWidth = 5.0;
		settings.watershedSeedsLocalMaximaDistanceThreshold = Double.MAX_VALUE;
		settings.watershedSeedsGlobalDistanceThreshold = 2.5;
		settings.interestPointsRadius = 0.5;
		settings.outputDirectory = null; //new File( path ).getParentFile();
		settings.inputDataSetName = "test";
		settings.returnEarly = true;
		settings.skeletonMaxLength = 600 * settings.workingVoxelSize;
		settings.minimalObjectSize = 200;  // um2
		settings.minimalTrackingSplittingObjectArea = 20; // this can be very small, um2
		settings.minimalObjectCenterDistance = 6;
		settings.maximalWatershedLength = 10;
		settings.closingRadius = 3;

		this.intensities = extractIntensities( imagePlus );

	}

	public void run()
	{

		ArrayList< RandomAccessibleInterval< T > > masks = createBinaryMasks( intensities );

		masks = splitTouchingObjects( intensities, masks );

		labelings = createTrackingBasedLabels( masks );

	}

	private ArrayList< RandomAccessibleInterval< T > > extractIntensities( ImagePlus imagePlus )
	{
		final Img< T > wrap = ( Img< T > ) ImageJFunctions.wrap( imagePlus );

		ArrayList<  RandomAccessibleInterval< T > > intensities = new ArrayList<>();
		for ( long t = 0; t < imagePlus.getNFrames() ; ++t )
		{
			// Channel
			final RandomAccessibleInterval< T > channel = Views.hyperSlice( wrap, 2, settings.microgliaChannelIndexOneBased - 1 );
			final RandomAccessibleInterval< T > timepoint = Views.hyperSlice( channel, 2, t );
			intensities.add( Utils.copyAsArrayImg( timepoint ) );
		}

		return intensities;
	}


	private ArrayList< RandomAccessibleInterval< T > > createTrackingBasedLabels( ArrayList< RandomAccessibleInterval< T > > masks )
	{
		final MaximalOverlapTracker maximalOverlapTracker = new MaximalOverlapTracker( masks );
		maximalOverlapTracker.run();
		return (ArrayList< RandomAccessibleInterval< T > > ) maximalOverlapTracker.getLabelings();
	}

	private ArrayList< RandomAccessibleInterval< T > > splitTouchingObjects( ArrayList< RandomAccessibleInterval< T > > intensities, ArrayList< RandomAccessibleInterval< T > > masks )
	{
		final TrackingSplitter splitter = new TrackingSplitter( masks, intensities, settings );
		splitter.run();
		masks = splitter.getSplitMasks();
		return masks;
	}

	private ArrayList< RandomAccessibleInterval< T > > createBinaryMasks( ArrayList< RandomAccessibleInterval< T > > intensities )
	{
		ArrayList<  RandomAccessibleInterval< T > > masks = new ArrayList<>();
		for ( long t = 0; t < intensities.size() ; ++t )
		{
			Utils.log("Creating mask for frame " + ( t + 1 ) );
			final SimpleSegmenter simpleSegmenter = new SimpleSegmenter( intensities.get( ( int ) t ), settings );
			simpleSegmenter.run();
			masks.add( simpleSegmenter.getMask() );
		}
		return masks;
	}


	public ArrayList< RandomAccessibleInterval< T > > getLabelings( )
	{
		return labelings;
	}


}
