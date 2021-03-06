package de.embl.cba.morphometry.segmentation;

import anisotropic_diffusion.Anisotropic_Diffusion_2D;
import de.embl.cba.morphometry.*;
import de.embl.cba.morphometry.microglia.MicrogliaSettings;
import de.embl.cba.morphometry.regions.Regions;
import ij.IJ;
import ij.ImageJ;
import ij.ImagePlus;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import static de.embl.cba.morphometry.viewing.BdvViewer.show;
import static de.embl.cba.transforms.utils.Scalings.createRescaledArrayImg;
import static de.embl.cba.transforms.utils.Transforms.getScalingFactors;


public class MicrogliaSegmenter< T extends RealType< T > & NativeType< T > >
{
	final MicrogliaSettings settings;
	private RandomAccessibleInterval< BitType > mask;
	final private RandomAccessibleInterval< T > intensity;
	final private boolean showIntermediateResults;

	public MicrogliaSegmenter(
			RandomAccessibleInterval< T > intensity,
			MicrogliaSettings settings )
	{
		this.intensity = intensity;
		this.settings = settings;
		this.showIntermediateResults = settings.showIntermediateResults;

	}

	public void run()
	{
		/**
		 *  Create working image
		 */
		final double[] workingCalibration = Utils.as2dDoubleArray( settings.workingVoxelSize );

		RandomAccessibleInterval< T > image = createRescaledArrayImg( intensity, getScalingFactors( settings.calibration2D, settings.workingVoxelSize ) );

		if ( showIntermediateResults ) show( image, "rescaled image", null, workingCalibration, false );


		/**
		 *  Smooth
		 */

		final ImagePlus wrap = ImageJFunctions.wrap( image, "" );
		final Anisotropic_Diffusion_2D diffusion2D = new Anisotropic_Diffusion_2D();
		diffusion2D.setup( "", wrap );
		final ImagePlus imagePlus = diffusion2D.runTD( wrap.getProcessor() );
		image = ImageJFunctions.wrapReal( imagePlus );

		if ( showIntermediateResults ) ImageJFunctions.show( image, "smoothed image" );



		/**
		 *  Compute offset and threshold
		 */

		final IntensityHistogram intensityHistogram = new IntensityHistogram( image, 65535, 2 );

		CoordinateAndValue mode = intensityHistogram.getMode();

		final CoordinateAndValue rightHandHalfMaximum = intensityHistogram.getRightHandHalfMode();

		double offset = mode.coordinate;
		double threshold = offset + ( rightHandHalfMaximum.coordinate - mode.coordinate ) * settings.thresholdInUnitsOfBackgroundPeakHalfWidth;

		Logger.debug( "Intensity offset: " + offset );
		Logger.debug( "Threshold: " + threshold );

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
