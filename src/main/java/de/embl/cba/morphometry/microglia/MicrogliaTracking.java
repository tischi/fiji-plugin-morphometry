package de.embl.cba.morphometry.microglia;

import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.segmentation.SimpleSegmenter;
import de.embl.cba.morphometry.splitting.TrackingSplitter;
import ij.ImagePlus;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import java.util.ArrayList;

public class MicrogliaTracking < T extends RealType< T > & NativeType< T > >
{
	final MicrogliaTrackingSettings settings;
	final ArrayList< RandomAccessibleInterval< T  > > intensities;
	private ArrayList< RandomAccessibleInterval< T > > labelings;

	public MicrogliaTracking ( ImagePlus imagePlus,
							   boolean showIntermediateResults,
							   OpService opService,
							   long microgliaChannelIndexOneBased,
							   long tMinOneBased,
							   long tMaxOneBased )
	{
		this.settings = configureMicrogliaTrackingSettings( imagePlus, showIntermediateResults, opService, microgliaChannelIndexOneBased, tMinOneBased, tMaxOneBased );

		this.intensities = Utils.get2DImagePlusMovieAsFrameList( imagePlus, settings.microgliaChannelIndexOneBased, tMinOneBased, tMaxOneBased );
	}

	public static MicrogliaTrackingSettings configureMicrogliaTrackingSettings( ImagePlus imagePlus, boolean showIntermediateResults, OpService opService, long microgliaChannelIndexOneBased, long tMin, long tMax )
	{
		MicrogliaTrackingSettings settings;
		settings = new MicrogliaTrackingSettings();

		settings.microgliaChannelIndexOneBased = 1;
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

		return settings;
	}

	public void run()
	{
		ArrayList< RandomAccessibleInterval< T > > masks = createBinaryMasks( intensities );

		labelings = splitTouchingObjectsAndTrack( intensities, masks );
	}

	private ArrayList< RandomAccessibleInterval< T > > createBinaryMasks( ArrayList< RandomAccessibleInterval< T > > intensities )
	{
		ArrayList<  RandomAccessibleInterval< T > > masks = new ArrayList<>();
		for ( long t = 0; t < intensities.size() ; ++t )
		{
			Logger.log("Creating mask for frame " + ( t + 1 ) );
			final SimpleSegmenter simpleSegmenter = new SimpleSegmenter( intensities.get( ( int ) t ), settings );
			simpleSegmenter.run();
			masks.add( simpleSegmenter.getMask() );
		}
		return masks;
	}

	private ArrayList< RandomAccessibleInterval< T > > splitTouchingObjectsAndTrack( ArrayList< RandomAccessibleInterval< T > > intensities, ArrayList< RandomAccessibleInterval< T > > masks )
	{
		final TrackingSplitter splitter = new TrackingSplitter( masks, intensities, settings );
		splitter.run();
		return splitter.getLabelings();
	}

	public ArrayList< RandomAccessibleInterval< T > > getLabelings( )
	{
		return labelings;
	}


}
