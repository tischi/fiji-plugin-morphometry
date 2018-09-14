package de.embl.cba.morphometry.microglia;

import de.embl.cba.morphometry.*;
import de.embl.cba.morphometry.objects.Measurements;
import net.imagej.ops.OpService;
import net.imglib2.*;
import net.imglib2.algorithm.morphology.distance.DistanceTransform;
import net.imglib2.algorithm.neighborhood.HyperSphereShape;
import net.imglib2.converter.Converters;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;


import java.util.*;

import static de.embl.cba.morphometry.Transforms.getScalingFactors;
import static de.embl.cba.morphometry.viewing.BdvImageViewer.show;


public class MicrogliaMorphometry< T extends RealType< T > & NativeType< T > >
{

	final MicrogliaMorphometrySettings settings;
	final OpService opService;

	private ArrayList< RandomAccessibleInterval< T > > resultImages;


	private HashMap< Integer, Map< String, Object > > objectMeasurements;

	public MicrogliaMorphometry( MicrogliaMorphometrySettings settings, OpService opService )
	{
		this.settings = settings;
		this.opService = opService;
	}

	public RandomAccessibleInterval< T > getResultImageStack()
	{
		return Views.stack( resultImages );
	}

	public HashMap< Integer, Map< String, Object > > getObjectMeasurements()
	{
		return objectMeasurements;
	}


	public void run()
	{

		/**
		 *  Create working image
		 */

		Utils.log( "Creating working resolution image..." );

		final double[] workingCalibration = Utils.get2dDoubleArray( settings.workingVoxelSize );

		final RandomAccessibleInterval< T > image = Algorithms.createIsotropicArrayImg( settings.image, getScalingFactors( settings.inputCalibration, settings.workingVoxelSize ) );

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

		double threshold = ( rightHandHalfMaximum.position - mode.position ) * settings.thresholdInUnitsOfBackgroundPeakHalfWidth;
		double offset = mode.position;
		Utils.log( "Offset: " + offset );
		Utils.log( "Threshold: " + ( threshold + offset ) );

		/**
		 * Create mask
		 */

		RandomAccessibleInterval< BitType > mask = Algorithms.createMask( image, threshold );

		if ( settings.showIntermediateResults ) show( mask, "mask", null, workingCalibration, false );

		/**
		 * Remove small objects from mask
		 */

		mask = Algorithms.removeSmallObjectsAndReturnMask( mask, settings.minimalObjectSize, settings.workingVoxelSize );

		if ( settings.showIntermediateResults ) show( mask, "size filtered mask", null, workingCalibration, false );

		/**
		 * Get objects
		 */

		final ImgLabeling< Integer, IntType > imgLabeling = Utils.asImgLabeling( mask );

		/**
		 * Compute skeleton
		 */

		// TODO: smooth boundaries to have less of a skeleton

		RandomAccessibleInterval< BitType > skeleton = opService.morphology().thinGuoHall(  mask );
		
		if ( settings.showIntermediateResults ) show( skeleton, "skeleton", null, workingCalibration, false );


		// final ArrayList< PositionAndValue > localMaxima = getLocalMaxima( image, new HyperSphereShape( 7 ), 0 );

		HashMap< Integer, Integer > numObjectsPerRegion = getNumObjectsFromSkeleton( imgLabeling, skeleton, settings );


		Algorithms.splitTouchingObjects(
				imgLabeling,
				numObjectsPerRegion,
				image,
				( int ) ( settings.minimalObjectCenterDistance / settings.workingVoxelSize ),
				opService );


		/**
		 * Compute object measurements
		 */

		objectMeasurements = new HashMap<>();

		Measurements.measureObjectSumIntensities( objectMeasurements, imgLabeling, image, "channel01" );

		Measurements.measureObjectSumIntensities( objectMeasurements, imgLabeling, skeleton, "skeleton" );

		Measurements.measureObjectPixelSizes( objectMeasurements, imgLabeling );

		Measurements.measureObjectPosition( objectMeasurements, imgLabeling, workingCalibration );

		Measurements.addGlobalBackgroundMeasurement( objectMeasurements, imgLabeling, offset );






		/**
		 * Compute branchpoints per object
		 */


		//final Img< BitType > branchpoints = Branchpoints.branchpoints( skeleton );





		/**
		 * Morphological closing
		 */


		RandomAccessibleInterval< BitType > closed = mask; //createClosedImage( mask );



		Utils.log( "Distance transform..." );

		final RandomAccessibleInterval< DoubleType > doubleBinary = Converters.convert( closed, ( i, o ) -> o.set( i.get() ? Double.MAX_VALUE : 0 ), new DoubleType() );

		final RandomAccessibleInterval< DoubleType > distance = ArrayImgs.doubles( Intervals.dimensionsAsLongArray( doubleBinary ) );

		DistanceTransform.transform( doubleBinary, distance, DistanceTransform.DISTANCE_TYPE.EUCLIDIAN, 1.0D );

		if ( settings.showIntermediateResults )
			show( distance, "distance transform", null, workingCalibration, false );

		/**
		 * Watershed seeds
		 */

		final ImgLabeling< Integer, IntType > seedsImgLabeling = createWatershedSeeds( workingCalibration, distance, closed );

		/**
		 * Watershed
		 */

		Utils.log( "Watershed..." );

		// prepare result label image
		final Img< IntType > watershedLabelImg = ArrayImgs.ints( Intervals.dimensionsAsLongArray( mask ) );
		final ImgLabeling< Integer, IntType > watershedImgLabeling = new ImgLabeling<>( watershedLabelImg );

		opService.image().watershed(
				watershedImgLabeling,
				Utils.invertedView( image ),
				seedsImgLabeling,
				false,
				false );



		if ( settings.splitTouchingObjects )
		{




			Utils.applyMask( watershedLabelImg, closed );

			if ( settings.showIntermediateResults )
				show( watershedLabelImg, "watershed", null, workingCalibration, false );

		}


		/**
		 * Generate output image
		 */


		resultImages = new ArrayList<>();
		resultImages.add( Utils.getEnlargedRai( image ) );
		resultImages.add( Utils.getEnlargedRai( ( RandomAccessibleInterval ) imgLabeling.getSource() ) );
		resultImages.add( Utils.getEnlargedRai( ( RandomAccessibleInterval ) skeleton ) );
	}

	public static HashMap< Integer, Integer > getNumObjectsFromSkeleton( ImgLabeling< Integer, IntType > imgLabeling, RandomAccessibleInterval< BitType > skeleton, MicrogliaMorphometrySettings settings )
	{
		HashMap< Integer, Map< String, Object > > skeletonMeasurements = new HashMap<>();
		Measurements.measureObjectSumIntensities( skeletonMeasurements, imgLabeling, skeleton, "skeleton" );
		HashMap< Integer, Integer > numObjects = new HashMap<>();
		for ( int label : skeletonMeasurements.keySet() )
		{
			final double skeletonLength = settings.workingVoxelSize * (long) skeletonMeasurements.get( label ).get( Measurements.SUM_INTENSITY + "_skeleton" );
			int n = (int) ( Math.ceil( skeletonLength / settings.skeletonMaxLength ) );
			numObjects.put( label, n );
		}
		return numObjects;
	}


	public ImgLabeling< Integer, IntType > createWatershedSeeds( double[] registrationCalibration,
																 RandomAccessibleInterval< DoubleType > distance,
																 RandomAccessibleInterval< BitType > mask )
	{
		Utils.log( "Seeds for watershed...");

		double globalDistanceThreshold = Math.pow( settings.watershedSeedsGlobalDistanceThreshold / settings.workingVoxelSize, 2 );
		double localMaximaDistanceThreshold = Math.pow( settings.watershedSeedsLocalMaximaDistanceThreshold / settings.workingVoxelSize, 2 );

		final RandomAccessibleInterval< BitType >  seeds = Algorithms.createSeeds(
				distance,
				new HyperSphereShape( 1 ),
				globalDistanceThreshold,
				localMaximaDistanceThreshold );

		final ImgLabeling< Integer, IntType > seedsLabelImg = Utils.asImgLabeling( seeds );

		if ( settings.showIntermediateResults ) show( Utils.asIntImg( seedsLabelImg ), "watershed seeds", null, registrationCalibration, false );

		return seedsLabelImg;
	}

}
