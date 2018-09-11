package de.embl.cba.morphometry.microglia;

import de.embl.cba.morphometry.*;
import de.embl.cba.morphometry.geometry.CentroidsParameters;
import de.embl.cba.morphometry.geometry.CoordinatesAndValues;
import de.embl.cba.morphometry.spindle.SpindleMorphometrySettings;
import net.imagej.ops.OpService;
import net.imglib2.*;
import net.imglib2.algorithm.morphology.Closing;
import net.imglib2.algorithm.morphology.distance.DistanceTransform;
import net.imglib2.algorithm.neighborhood.HyperSphereShape;
import net.imglib2.algorithm.neighborhood.Neighborhood;
import net.imglib2.algorithm.neighborhood.Shape;
import net.imglib2.converter.Converters;
import net.imglib2.histogram.Histogram1d;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static de.embl.cba.morphometry.Constants.X;
import static de.embl.cba.morphometry.Constants.Z;
import static de.embl.cba.morphometry.Transforms.getScalingFactors;
import static de.embl.cba.morphometry.viewing.BdvImageViewer.show;
import static java.lang.Math.toRadians;


public class MicrogliaMorphometry< T extends RealType< T > & NativeType< T > >
{

	final MicrogliaMorphometrySettings settings;
	final OpService opService;

	public MicrogliaMorphometry( MicrogliaMorphometrySettings settings, OpService opService )
	{
		this.settings = settings;
		this.opService = opService;
	}

	public RandomAccessibleInterval< T > run()
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

		PositionAndValue mode = intensityHistogram.getMode();

		final PositionAndValue rightHandHalfMaximum = intensityHistogram.getRightHandHalfMaximum();

		double threshold = ( rightHandHalfMaximum.position - mode.position ) * settings.thresholdInUnitsOfBackgroundPeakHalfWidth;

		Utils.log( "Offset: " + mode.position );
		Utils.log( "Threshold: " + ( threshold + mode.position ) );

		/**
		 * Create mask
		 */

		RandomAccessibleInterval< BitType > mask = createMask( image, threshold );

		if ( settings.showIntermediateResults ) show( mask, "mask", null, workingCalibration, false );

		/**
		 * Remove small objects from mask
		 */

		mask = Algorithms.removeSmallObjectsAndReturnMask( mask, settings.minimalObjectSize, settings.workingVoxelSize );

		if ( settings.showIntermediateResults ) show( mask, "size filtered mask", null, workingCalibration, false );

		/**
		 * Morphological closing
		 */

		RandomAccessibleInterval< BitType > closed = mask; //createClosedImage( mask );

//		if ( settings.showIntermediateResults ) show( closed, "closed", null, workingCalibration, false );


		/**
		 * Distance transform
		 *
		 * Note: EUCLIDIAN distances are returned as squared distances
		 */

		Utils.log( "Distance transform..." );

		final RandomAccessibleInterval< DoubleType > doubleBinary = Converters.convert( closed, ( i, o ) -> o.set( i.get() ? Double.MAX_VALUE : 0 ), new DoubleType() );

		final RandomAccessibleInterval< DoubleType > distance = ArrayImgs.doubles( Intervals.dimensionsAsLongArray( doubleBinary ) );

		DistanceTransform.transform( doubleBinary, distance, DistanceTransform.DISTANCE_TYPE.EUCLIDIAN, 1.0D );

		if ( settings.showIntermediateResults ) show( distance, "distance transform", null, workingCalibration, false );

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

		Utils.applyMask( watershedLabelImg, closed );

		if ( settings.showIntermediateResults ) show( watershedLabelImg, "watershed", null, workingCalibration, false );

		/**
		 * Generate output image
		 */

		ArrayList< RandomAccessibleInterval< T > > randomAccessibleIntervals = new ArrayList<>();
		randomAccessibleIntervals.add( image );
		randomAccessibleIntervals.add( (RandomAccessibleInterval) seedsImgLabeling.getSource() );
		randomAccessibleIntervals.add( (RandomAccessibleInterval) watershedImgLabeling.getSource() );
		RandomAccessibleInterval< T > stack = Views.stack( randomAccessibleIntervals );

		return stack;
	}


	public void convertToOriginalImagePixelUnits( AffineTransform3D alignmentTransform, double[] spindlePole )
	{
		Utils.divide( spindlePole, settings.workingVoxelSize );
		alignmentTransform.inverse().apply( spindlePole, spindlePole );
	}

	public void drawPoint( RandomAccessibleInterval< T > rai, double[] position, double radius, double calibration )
	{
		Shape shape = new HyperSphereShape( (int) Math.ceil( radius / calibration ) );
		final RandomAccessible< Neighborhood< T > > nra = shape.neighborhoodsRandomAccessible( rai );
		final RandomAccess< Neighborhood< T > > neighborhoodRandomAccess = nra.randomAccess();

		neighborhoodRandomAccess.setPosition( Utils.asLongs( position )  );
		final Neighborhood< T > neighborhood = neighborhoodRandomAccess.get();

		final Cursor< T > cursor = neighborhood.cursor();
		while( cursor.hasNext() )
		{
			try
			{
				cursor.next().setReal( 200 );
			}
			catch ( ArrayIndexOutOfBoundsException e )
			{
				Utils.log( "[ERROR] Draw points out of bounds..." );
				break;
			}
		}
	}

	public static double[] getLeftAndRightMaxLocs( CoordinatesAndValues tubulinProfile, ArrayList< Double > tubulinProfileAbsoluteDerivative )
	{
		double[] rangeMinMax = new double[ 2 ];
		double[] maxLocs = new double[ 2 ];

		rangeMinMax[ 0 ] = 0;
		rangeMinMax[ 1 ] = Double.MAX_VALUE;
		maxLocs[ 0 ] = Utils.computeMaxLoc( tubulinProfile.coordinates, tubulinProfileAbsoluteDerivative, rangeMinMax );

		rangeMinMax[ 0 ] = - Double.MAX_VALUE;
		rangeMinMax[ 1 ] = 0;
		maxLocs[ 1 ] = Utils.computeMaxLoc( tubulinProfile.coordinates, tubulinProfileAbsoluteDerivative, rangeMinMax );
		return maxLocs;
	}

	public RandomAccessibleInterval< BitType > createClosedImage( RandomAccessibleInterval< BitType > mask )
	{
		RandomAccessibleInterval< BitType > closed = Utils.copyAsArrayImg( mask );

		if ( settings.closingRadius > 0 )
		{
			Utils.log( "Morphological closing...");
			Shape closingShape = new HyperSphereShape( ( int ) ( settings.closingRadius / settings.workingVoxelSize ) );
			Closing.close( Views.extendBorder( mask ), Views.iterable( closed ), closingShape, 1 );
		}

		return closed;
	}

	public < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< BitType > createMask( RandomAccessibleInterval< T > downscaled, double threshold )
	{
		Utils.log( "Creating mask...");

		RandomAccessibleInterval< BitType > mask =
				Converters.convert( downscaled, ( i, o )
						-> o.set( i.getRealDouble() > threshold ? true : false ), new BitType() );

		mask = opService.morphology().fillHoles( mask );

		return mask;
	}

	public < T extends RealType< T > & NativeType< T > > double getThreshold( RandomAccessibleInterval< T > downscaled )
	{

		double threshold = 0;

		if ( settings.thresholdModality.equals( SpindleMorphometrySettings.HUANG_AUTO_THRESHOLD ) )
		{
			final Histogram1d< T > histogram = opService.image().histogram( Views.iterable( downscaled ) );

			double huang = opService.threshold().huang( histogram ).getRealDouble();
			double yen = opService.threshold().yen( histogram ).getRealDouble();

			threshold = huang;
		}
		else
		{
			threshold= settings.thresholdInUnitsOfBackgroundPeakHalfWidth;
		}
		return threshold;
	}

	public ImgLabeling< Integer, IntType > createWatershedSeeds( double[] registrationCalibration,
																 RandomAccessibleInterval< DoubleType > distance,
																 RandomAccessibleInterval< BitType > mask )
	{
		Utils.log( "Seeds for watershed...");

		double globalDistanceThreshold = Math.pow( settings.watershedSeedsGlobalDistanceThreshold / settings.workingVoxelSize, 2 );
		double localMaximaDistanceThreshold = Math.pow( settings.watershedSeedsLocalMaximaDistanceThreshold / settings.workingVoxelSize, 2 );

		final RandomAccessibleInterval< BitType >  seeds = Utils.createSeeds(
				distance,
				new HyperSphereShape( 1 ),
				globalDistanceThreshold,
				localMaximaDistanceThreshold );

		final ImgLabeling< Integer, IntType > seedsLabelImg = Algorithms.createImgLabeling( seeds );

		if ( settings.showIntermediateResults ) show( Utils.asIntImg( seedsLabelImg ), "watershed seeds", null, registrationCalibration, false );

		return seedsLabelImg;
	}

	public AffineTransform3D computeOrientationTransform( RandomAccessibleInterval yawAlignedMask, RandomAccessibleInterval yawAlignedIntensities, double calibration )
	{
		final CoordinatesAndValues coordinatesAndValues = Utils.computeAverageIntensitiesAlongAxis( yawAlignedIntensities, yawAlignedMask, X, calibration );

		if ( settings.showIntermediateResults ) Plots.plot( coordinatesAndValues.coordinates, coordinatesAndValues.values, "x", "average intensity" );

		double maxLoc = Utils.computeMaxLoc( coordinatesAndValues.coordinates, coordinatesAndValues.values, null );

		AffineTransform3D affineTransform3D = new AffineTransform3D();

		if ( maxLoc < 0 ) affineTransform3D.rotate( Z, toRadians( 180.0D ) );

		return affineTransform3D;
	}

	public ArrayList< RealPoint > createTransformedCentroidPointList( CentroidsParameters centroidsParameters, AffineTransform3D rollTransform )
	{
		final ArrayList< RealPoint > transformedRealPoints = new ArrayList<>();

		for ( RealPoint realPoint : centroidsParameters.centroids )
		{
			final RealPoint transformedRealPoint = new RealPoint( 0, 0, 0 );
			rollTransform.apply( realPoint, transformedRealPoint );
			transformedRealPoints.add( transformedRealPoint );
		}
		return transformedRealPoints;
	}

	public static AffineTransform3D computeRollTransform( CentroidsParameters centroidsParameters, SpindleMorphometrySettings settings )
	{
		final double rollAngle = computeRollAngle( centroidsParameters, settings.rollAngleMinDistanceToAxis, settings.rollAngleMinDistanceToCenter, settings.rollAngleMaxDistanceToCenter );

		Utils.log( "Roll angle " + rollAngle );

		AffineTransform3D rollTransform = new AffineTransform3D();

		rollTransform.rotate( X, - toRadians( rollAngle ) );

		return rollTransform;
	}

	public static double computeRollAngle( CentroidsParameters centroidsParameters, double minDistanceToAxis, double minDistanceToCenter, double maxDistanceToCenter )
	{
		final int n = centroidsParameters.axisCoordinates.size();

		List< Double> offCenterAngles = new ArrayList<>(  );

		for ( int i = 0; i < n; ++i )
		{
			if ( ( centroidsParameters.distances.get( i ) > minDistanceToAxis ) &&
					( Math.abs(  centroidsParameters.axisCoordinates.get( i ) ) > minDistanceToCenter ) &&
					( Math.abs(  centroidsParameters.axisCoordinates.get( i ) ) < maxDistanceToCenter ))
			{
				offCenterAngles.add( centroidsParameters.angles.get( i ) );
			}
		}

		Collections.sort( offCenterAngles );

		double medianAngle = Utils.median( offCenterAngles );

		return medianAngle;
	}

	private AffineTransform3D createFinalTransform( double[] inputCalibration, AffineTransform3D registration, double[] registrationCalibration )
	{
		final AffineTransform3D transform =
				Transforms.getScalingTransform( inputCalibration, settings.workingVoxelSize )
				.preConcatenate( registration )
				.preConcatenate( Transforms.getScalingTransform( registrationCalibration, settings.outputResolution ) );

		return transform;
	}


	private double[] getRegistrationCalibration()
	{
		double[] registrationCalibration = new double[ 3 ];
		Arrays.fill( registrationCalibration, settings.workingVoxelSize );
		return registrationCalibration;
	}

	//
	// Useful code snippets
	//

	/**
	 *  Distance transformAllChannels

	 Hi Christian

	 yes, it seems that you were doing the right thing (as confirmed by your
	 visual inspection of the result). One thing to note: You should
	 probably use a DoubleType image with 1e20 and 0 values, to make sure
	 that f(q) is larger than any possible distance in your image. If you
	 choose 255, your distance is effectively bounded at 255. This can be an
	 issue for big images with sparse foreground objects. With squared
	 Euclidian distance, 255 is already reached if a background pixels is
	 further than 15 pixels from a foreground pixel! If you use
	 Converters.convert to generate your image, the memory consumption
	 remains the same.


	 Phil

	 final RandomAccessibleInterval< UnsignedByteType > binary = Converters.convert(
	 downscaled, ( i, o ) -> o.set( i.getRealDouble() > settings.thresholdInUnitsOfBackgroundPeakHalfWidth ? 255 : 0 ), new UnsignedByteType() );

	 if ( settings.showIntermediateResults ) show( binary, "binary", null, calibration, false );


	 final RandomAccessibleInterval< DoubleType > distance = ArrayImgs.doubles( Intervals.dimensionsAsLongArray( binary ) );

	 DistanceTransform.transformAllChannels( binary, distance, DistanceTransform.DISTANCE_TYPE.EUCLIDIAN );


	 final double maxDistance = Algorithms.getMaximumValue( distance );

	 final RandomAccessibleInterval< IntType > invertedDistance = Converters.convert( distance, ( i, o ) -> {
	 o.set( ( int ) ( maxDistance - i.get() ) );
	 }, new IntType() );

	 if ( settings.showIntermediateResults ) show( invertedDistance, "distance", null, calibration, false );

	 */


	/**
	 * Convert ImgLabelling to Rai

	 final RandomAccessibleInterval< IntType > labelMask =
	 Converters.convert( ( RandomAccessibleInterval< LabelingType< Integer > > ) watershedImgLabeling,
	 ( i, o ) -> {
	 o.set( i.getIndex().getInteger() );
	 }, new IntType() );

	 */




}
