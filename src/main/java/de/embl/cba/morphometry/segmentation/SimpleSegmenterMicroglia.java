package de.embl.cba.morphometry.segmentation;

import de.embl.cba.morphometry.*;
import de.embl.cba.morphometry.microglia.MicrogliaSegmentationAndTrackingSettings;
import de.embl.cba.morphometry.regions.Regions;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import static de.embl.cba.morphometry.viewing.BdvViewer.show;
import static de.embl.cba.transforms.utils.Scalings.createRescaledArrayImg;
import static de.embl.cba.transforms.utils.Transforms.getScalingFactors;


public class SimpleSegmenterMicroglia< T extends RealType< T > & NativeType< T > >
{

	final MicrogliaSegmentationAndTrackingSettings settings;
	private RandomAccessibleInterval< BitType > mask;
	final private RandomAccessibleInterval< T > intensity;
	final private boolean showIntermediateResults;

	public SimpleSegmenterMicroglia(
			RandomAccessibleInterval< T > intensity,
			MicrogliaSegmentationAndTrackingSettings settings )
	{
		this.intensity = intensity;
		this.settings = settings;
		this.showIntermediateResults = false; //settings.showIntermediateResults;

	}

	public void run()
	{

		/**
		 *  Create working image
		 */

		final double[] workingCalibration = Utils.as2dDoubleArray( settings.workingVoxelSize );

		final RandomAccessibleInterval< T > image =
				createRescaledArrayImg( intensity,
				getScalingFactors( settings.inputCalibration, settings.workingVoxelSize ) );

		if ( showIntermediateResults ) show( image, "image isotropic resolution", null, workingCalibration, false );


		/**
		 *  Smooth
		 */

		// TODO


		/**
		 *  Compute offset and threshold
		 */

		final IntensityHistogram intensityHistogram = new IntensityHistogram( image, 65535, 2 );

		CoordinateAndValue mode = intensityHistogram.getMode();

		final CoordinateAndValue rightHandHalfMaximum = intensityHistogram.getRightHandHalfMaximum();

		double threshold = ( rightHandHalfMaximum.coordinate - mode.coordinate ) * settings.thresholdInUnitsOfBackgroundPeakHalfWidth;
		double offset = mode.coordinate;
		Logger.debug( "Intensity offset: " + offset );
		Logger.debug( "Threshold: " + ( threshold + offset ) );

		/**
		 * Create mask
		 */

		mask = Algorithms.createMask( image, threshold );

		if ( showIntermediateResults ) show( mask, "mask", null, workingCalibration, false );


		/**
		 * Remove small objects from mask
		 */

		Regions.removeSmallRegionsInMask( mask, settings.minimalObjectSize, settings.workingVoxelSize );

		if ( showIntermediateResults ) show( mask, "size filtered mask", null, workingCalibration, false );


	}

	public RandomAccessibleInterval< BitType > getMask()
	{
		return mask;
	}


}