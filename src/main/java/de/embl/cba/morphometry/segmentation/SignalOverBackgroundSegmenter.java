package de.embl.cba.morphometry.segmentation;

import de.embl.cba.morphometry.*;
import de.embl.cba.morphometry.regions.Regions;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;

import static de.embl.cba.morphometry.viewing.BdvViewer.show;


public class SignalOverBackgroundSegmenter< T extends RealType< T > & NativeType< T > >
{
	final private RandomAccessibleInterval< T > intensity;
	final private double signalToNoise;
	final private long minimalRegionSize;

	private RandomAccessibleInterval< BitType > mask;

	public SignalOverBackgroundSegmenter(
			RandomAccessibleInterval< T > intensity,
			double signalToNoise, long minimalRegionSize )
	{
		this.intensity = intensity;
		this.signalToNoise = signalToNoise;
		this.minimalRegionSize = minimalRegionSize;
	}

	public RandomAccessibleInterval< BitType > createMask()
	{

		/**
		 *  Compute offset and threshold
		 */

		final IntensityHistogram intensityHistogram =
				new IntensityHistogram( intensity, 65535, 2 );

		CoordinateAndValue mode = intensityHistogram.getMode();

		final CoordinateAndValue rightHandHalfMaximum = intensityHistogram.getRightHandHalfMode();

		double threshold = mode.coordinate + ( rightHandHalfMaximum.coordinate - mode.coordinate ) * signalToNoise;


		/**
		 * Create mask
		 */

		mask = Algorithms.createMask( intensity, threshold );

		/**
		 * Remove small objects from mask
		 */

		Regions.removeSmallRegionsInMask( mask, minimalRegionSize );

		return mask;

	}

}
