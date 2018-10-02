package de.embl.cba.morphometry.segmentation;

import de.embl.cba.morphometry.Algorithms;
import de.embl.cba.morphometry.CoordinateAndValue;
import de.embl.cba.morphometry.IntensityHistogram;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.microglia.MicrogliaSettings;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;

import static de.embl.cba.morphometry.Transforms.getScalingFactors;
import static de.embl.cba.morphometry.viewing.BdvImageViewer.show;


public class SimpleSegmenter< T extends RealType< T > & NativeType< T > >
{

	final MicrogliaSettings settings;
	private RandomAccessibleInterval< BitType > mask;
	final private RandomAccessibleInterval< T > intensity;

	public SimpleSegmenter( RandomAccessibleInterval< T > intensity, MicrogliaSettings settings )
	{
		this.intensity = intensity;
		this.settings = settings;
	}

	public void run()
	{

		/**
		 *  Create working image
		 */

		Utils.log( "Creating working resolution image..." );

		final double[] workingCalibration = Utils.as2dDoubleArray( settings.workingVoxelSize );

		final RandomAccessibleInterval< T > image = Algorithms.createIsotropicArrayImg( intensity, getScalingFactors( settings.inputCalibration, settings.workingVoxelSize ) );

		if ( settings.showIntermediateResults ) show( image, "image isotropic resolution", null, workingCalibration, false );


		/**
		 *  Smooth
		 */

		// TODO


		/**
		 *  Compute offset and threshold
		 */

		Utils.log( "Computing offset and threshold..." );

		final IntensityHistogram intensityHistogram = new IntensityHistogram( image, settings.maxPossibleValueInDataSet, 2 );

		CoordinateAndValue mode = intensityHistogram.getMode();

		final CoordinateAndValue rightHandHalfMaximum = intensityHistogram.getRightHandHalfMaximum();

		double threshold = ( rightHandHalfMaximum.coordinate - mode.coordinate ) * settings.thresholdInUnitsOfBackgroundPeakHalfWidth;
		double offset = mode.coordinate;
		Utils.log( "Intensity offset: " + offset );
		Utils.log( "Threshold: " + ( threshold + offset ) );

		/**
		 * Create mask
		 */

		mask = Algorithms.createMask( image, threshold );

		if ( settings.showIntermediateResults ) show( mask, "mask", null, workingCalibration, false );


		/**
		 * Remove small objects from mask
		 */

		mask = Algorithms.removeSmallObjectsAndReturnMask( mask, settings.minimalObjectSize, settings.workingVoxelSize );

		if ( settings.showIntermediateResults ) show( mask, "size filtered mask", null, workingCalibration, false );


	}

	public RandomAccessibleInterval< BitType > getMask()
	{
		return mask;
	}


}
