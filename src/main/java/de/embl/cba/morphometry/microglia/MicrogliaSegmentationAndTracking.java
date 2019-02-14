package de.embl.cba.morphometry.microglia;

import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.segmentation.SimpleSegmenterMicroglia;
import de.embl.cba.morphometry.tracking.SemiAutomatedTrackingSplitter;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import java.util.ArrayList;

public class MicrogliaSegmentationAndTracking< T extends RealType< T > & NativeType< T > >
{
	final MicrogliaSegmentationAndTrackingSettings settings;
	final ArrayList< RandomAccessibleInterval< T  > > intensities;
	private ArrayList< RandomAccessibleInterval< T > > labelings;

	public MicrogliaSegmentationAndTracking( ArrayList< RandomAccessibleInterval< T  > > intensities,
											 double[] calibration,
											 String outputLabelingsPath,
											 boolean showIntermediateResults,
											 OpService opService )
	{
		this.intensities = intensities;
		this.settings = configureSettings( calibration, outputLabelingsPath, showIntermediateResults, opService );
	}

	public static MicrogliaSegmentationAndTrackingSettings configureSettings(
			double[] calibration,
			String outputLabelingsPath,
			boolean showIntermediateResults,
			OpService opService )
	{
		MicrogliaSegmentationAndTrackingSettings settings;
		settings = new MicrogliaSegmentationAndTrackingSettings();

		settings.showIntermediateResults = showIntermediateResults;
		settings.opService = opService;
		settings.inputCalibration = calibration;
		settings.workingVoxelSize = settings.inputCalibration[ 0 ];
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
		settings.outputLabelingsPath = outputLabelingsPath;

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
			final SimpleSegmenterMicroglia simpleSegmenterMicroglia = new SimpleSegmenterMicroglia( intensities.get( ( int ) t ), settings );
			simpleSegmenterMicroglia.run();
			masks.add( simpleSegmenterMicroglia.getMask() );
		}
		return masks;
	}

	private ArrayList< RandomAccessibleInterval< T > > splitTouchingObjectsAndTrack( ArrayList< RandomAccessibleInterval< T > > intensities, ArrayList< RandomAccessibleInterval< T > > masks )
	{
		final SemiAutomatedTrackingSplitter splitter = new SemiAutomatedTrackingSplitter( masks, intensities, settings );
		splitter.run();
		return splitter.getLabelings();
	}

	public ArrayList< RandomAccessibleInterval< T > > getLabelings( )
	{
		return labelings;
	}


}
