package de.embl.cba.morphometry.drosophila.shavenbaby;

import de.embl.cba.morphometry.*;
import de.embl.cba.morphometry.geometry.CentroidsParameters;
import de.embl.cba.morphometry.geometry.CoordinatesAndValues;
import de.embl.cba.morphometry.geometry.EllipsoidParameters;
import de.embl.cba.morphometry.geometry.Ellipsoids;
import net.imagej.ops.OpService;
import net.imglib2.*;
import net.imglib2.algorithm.morphology.distance.DistanceTransform;
import net.imglib2.algorithm.neighborhood.HyperSphereShape;
import net.imglib2.converter.Converters;
import net.imglib2.histogram.Histogram1d;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.interpolation.randomaccess.NearestNeighborInterpolatorFactory;
import net.imglib2.outofbounds.OutOfBoundsConstantValueFactory;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.roi.labeling.LabelRegion;
import net.imglib2.roi.labeling.LabelRegions;
import net.imglib2.roi.labeling.LabelingType;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.integer.UnsignedByteType;
import net.imglib2.type.numeric.integer.UnsignedIntType;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.algorithm.neighborhood.Shape;
import net.imglib2.algorithm.morphology.Closing;

import net.imglib2.util.Intervals;
import net.imglib2.view.Views;

import java.util.*;

import static de.embl.cba.morphometry.Constants.*;
import static de.embl.cba.morphometry.Transforms.getScalingFactors;
import static de.embl.cba.morphometry.drosophila.dapi.DapiRegistration.createXAxisRollTransform;
import static de.embl.cba.morphometry.viewing.BdvImageViewer.show;
import static java.lang.Math.toRadians;


public class ShavenBabyRegistration
{

	final ShavenBabyRegistrationSettings settings;
	final OpService opService;
	private RandomAccessibleInterval< BitType > centralObjectMask;

	public ShavenBabyRegistration( ShavenBabyRegistrationSettings settings, OpService opService )
	{
		this.settings = settings;
		this.opService = opService;
	}

	public RandomAccessibleInterval< BitType > getCentralObjectMask()
	{
		return centralObjectMask;
	}

	public < T extends RealType< T > & NativeType< T > >
	AffineTransform3D computeRegistration( RandomAccessibleInterval< T > svb, // Shavenbaby
										   RandomAccessibleInterval< T > ama, // Amnioserosa
										   double[] inputCalibration  )
	{

		AffineTransform3D registration = new AffineTransform3D();

		double[] registrationCalibration = Utils.get3dDoubleArray( settings.registrationResolution );

		Utils.log( "Refractive index scaling correction..." );

		RefractiveIndexMismatchCorrections.correctCalibration( inputCalibration, settings.refractiveIndexScalingCorrectionFactor );


		/**
		 *  Down-sampling to registration resolution
		 */

		// TODO: during downsampling saturated pixels become zero...
		
		Utils.log( "Down-sampling to registration resolution..." );

		final RandomAccessibleInterval< T > downscaledSvb = Algorithms.createIsotropicArrayImg( svb, getScalingFactors( inputCalibration, settings.registrationResolution ) );
		final RandomAccessibleInterval< T > downscaledAma = Algorithms.createIsotropicArrayImg( ama, getScalingFactors( inputCalibration, settings.registrationResolution ) );

		if ( settings.showIntermediateResults ) show( downscaledSvb, "at registration resolution", null, registrationCalibration, false );


		/**
		 *  Compute offset
		 */


		Utils.log( "Computing offset and threshold..." );

		final IntensityHistogram rawDataIntensityHistogram = new IntensityHistogram( downscaledSvb, 65535.0, 5.0 );

		CoordinateAndValue mode = rawDataIntensityHistogram.getMode();

		Utils.log( "Offset: " + mode.position );


		/**
		 *  Refractive index corrections
		 */
		
		Utils.log( "Refractive index intensity correction..." );

		final RandomAccessibleInterval< T > intensityCorrectedSvb = Utils.copyAsArrayImg( downscaledSvb );
		final RandomAccessibleInterval< T > intensityCorrectedAma = Utils.copyAsArrayImg( downscaledAma );

		RefractiveIndexMismatchCorrections.correctIntensity( intensityCorrectedSvb, registrationCalibration[ Z ], mode.position, settings.refractiveIndexIntensityCorrectionDecayLength );
		RefractiveIndexMismatchCorrections.correctIntensity( intensityCorrectedAma, registrationCalibration[ Z ], mode.position, settings.refractiveIndexIntensityCorrectionDecayLength );

		if ( settings.showIntermediateResults ) show( intensityCorrectedSvb, "intensity corrected svb", null, registrationCalibration, false );
		//if ( settings.showIntermediateResults ) show( intensityCorrectedAma, "intensity corrected ama", null, registrationCalibration, false );



		/**
		 *  Compute threshold
		 */

		final IntensityHistogram correctedIntensityHistogram =
				new IntensityHistogram(
				intensityCorrectedSvb,
				65535.0,
				5.0 );


		final double huang = opService.threshold().huang( opService.image().histogram( Views.iterable( intensityCorrectedSvb ) ) ).getRealDouble();
		final double otsu = opService.threshold().otsu( opService.image().histogram( Views.iterable( intensityCorrectedSvb ) ) ).getRealDouble();
		final double yen = opService.threshold().yen( opService.image().histogram( Views.iterable( intensityCorrectedSvb ) ) ).getRealDouble();

		// Utils.log( "Intensity corrected threshold: " + thresholdAfterIntensityCorrection );


		double thresholdAfterIntensityCorrection = huang; // TODO

		Utils.log( "Intensity corrected threshold: " + thresholdAfterIntensityCorrection );



		/**
		 * Create mask
		 */

		RandomAccessibleInterval< BitType > mask = createMask( intensityCorrectedSvb, thresholdAfterIntensityCorrection );

		mask = opService.morphology().fillHoles( mask );

//		if ( settings.showIntermediateResults ) show( mask, "mask", null, registrationCalibration, false );

		/**
		 * Remove small objects
		 */

		mask = Algorithms.removeSmallObjectsAndReturnMask( mask, settings.minimalObjectSize, settings.registrationResolution );

		if ( settings.showIntermediateResults ) show( mask, "closed and object size filtered", null, registrationCalibration, false );

		/**
		 * Distance transform
		 *
		 * Note: EUCLIDIAN distances are returned as squared distances
		 */

		Utils.log( "Distance transform..." );

		final RandomAccessibleInterval< DoubleType > doubleBinary = Converters.convert( mask, ( i, o ) -> o.set( i.get() ? Double.MAX_VALUE : 0 ), new DoubleType() );

		final RandomAccessibleInterval< DoubleType > distance = ArrayImgs.doubles( Intervals.dimensionsAsLongArray( doubleBinary ) );

		DistanceTransform.transform( doubleBinary, distance, DistanceTransform.DISTANCE_TYPE.EUCLIDIAN, 1.0D );

//		if ( settings.showIntermediateResults ) show( distance, "distance transform", null, registrationCalibration, false );

		/**
		 * Watershed seeds
		 */

		final ImgLabeling< Integer, IntType > seedsLabelImg = createWatershedSeeds( registrationCalibration, distance, mask );


		/**
		 * Watershed
		 */

		Utils.log( "Watershed..." );

		// prepare result label image
		final Img< IntType > watershedLabelImg = ArrayImgs.ints( Intervals.dimensionsAsLongArray( mask ) );
		final ImgLabeling< Integer, IntType > watershedLabeling = new ImgLabeling<>( watershedLabelImg );

		opService.image().watershed(
				watershedLabeling,
				Utils.invertedView( distance ),
				seedsLabelImg,
				false,
				false );

		Utils.applyMask( watershedLabelImg, mask );

		if ( settings.showIntermediateResults ) show( watershedLabelImg, "watershed", null, registrationCalibration, false );


		/**
		 * Get central embryo
		 */

		Utils.log( "Get central embryo..." );

		final LabelRegion< Integer > centralObjectRegion = getCentralObjectLabelRegion( watershedLabeling );

		if ( centralObjectRegion == null )
		{
			return null;
		}

		centralObjectMask = Algorithms.createMaskFromLabelRegion( centralObjectRegion, Intervals.dimensionsAsLongArray( downscaledSvb ) );

		centralObjectMask = close( centralObjectMask, ( int ) ( 20 / settings.registrationResolution ) );

		centralObjectMask = opService.morphology().fillHoles( centralObjectMask );

		if ( settings.showIntermediateResults ) show( centralObjectMask, "central object", null, registrationCalibration, false );


		/**
		 * Compute ellipsoid (probably mainly yaw) alignment
		 */

		Utils.log( "Fit ellipsoid..." );

		final EllipsoidParameters ellipsoidParameters = Ellipsoids.computeParametersFromBinaryImage( centralObjectMask );

		registration.preConcatenate( Ellipsoids.createAlignmentTransform( ellipsoidParameters ) );

		final RandomAccessibleInterval yawAlignedMask = Utils.copyAsArrayImg( Transforms.createTransformedView( centralObjectMask, registration, new NearestNeighborInterpolatorFactory() ) );

		final RandomAccessibleInterval yawAlignedIntensities = Utils.copyAsArrayImg( Transforms.createTransformedView( downscaledSvb, registration ) );


		/**
		 *  Long axis orientation
		 */

		Utils.log( "Computing long axis orientation..." );

		final AffineTransform3D orientationTransform = computeFlippingTransform( yawAlignedMask, yawAlignedIntensities, settings.registrationResolution );

		registration = registration.preConcatenate( orientationTransform );

		final RandomAccessibleInterval< BitType > yawAndOrientationAlignedMask = Utils.copyAsArrayImg( Transforms.createTransformedView( centralObjectMask, registration, new NearestNeighborInterpolatorFactory() ) );

		if ( settings.showIntermediateResults ) show( yawAndOrientationAlignedMask, "long axis aligned svb", null, registrationCalibration, false );


		/**
		 *  Roll transform
		 */

		registration = computeRollTransform( registration, registrationCalibration, intensityCorrectedSvb, intensityCorrectedAma, yawAndOrientationAlignedMask );


		/**
		 * Show aligned input image at registration resolution
		 */

		if ( settings.showIntermediateResults ) show( Transforms.createTransformedView( intensityCorrectedSvb, registration ), "aligned input data ( " + settings.outputResolution + " um )", origin(), registrationCalibration, false );

		/**
		 * Compute final registration
		 */

		registration = createFinalTransform( inputCalibration, registration, registrationCalibration );

		return registration;

	}

	public < T extends RealType< T > & NativeType< T > > AffineTransform3D computeRollTransform( AffineTransform3D registration, double[] registrationCalibration, RandomAccessibleInterval< T > intensityCorrectedSvb, RandomAccessibleInterval< T > intensityCorrectedAma, RandomAccessibleInterval< BitType > yawAndOrientationAlignedMask )
	{
		Utils.log( "Computing roll transform, using method: " + settings.rollAngleComputationMethod );

		if ( settings.rollAngleComputationMethod.equals( ShavenBabyRegistrationSettings.AMNIOSEROSA ) )
		{
			final RandomAccessibleInterval yawAndOrientationAlignedAma = Utils.copyAsArrayImg( Transforms.createTransformedView( intensityCorrectedAma, registration, new NearestNeighborInterpolatorFactory() ) );
			final RandomAccessibleInterval yawAndOrientationAlignedSvb = Utils.copyAsArrayImg( Transforms.createTransformedView( intensityCorrectedSvb, registration, new NearestNeighborInterpolatorFactory() ) );

			final AffineTransform3D intensityBasedRollTransform = computeIntensityBasedRollTransform(
					yawAndOrientationAlignedAma,
					settings.amaProjectionXMin,
					settings.amaProjectionXMax,
					settings.amaProjectionBlurSigma );

			registration = registration.preConcatenate( intensityBasedRollTransform );

		}
		else if ( settings.rollAngleComputationMethod.equals( ShavenBabyRegistrationSettings.CENTROID_SHAPE ) )
		{
			final CentroidsParameters centroidsParameters = Utils.computeCentroidsParametersAlongXAxis( yawAndOrientationAlignedMask, settings.registrationResolution, settings.rollAngleMaxDistanceToCenter );

			if ( settings.showIntermediateResults )
				Plots.plot( centroidsParameters.axisCoordinates, centroidsParameters.angles, "x", "angle" );
			if ( settings.showIntermediateResults )
				Plots.plot( centroidsParameters.axisCoordinates, centroidsParameters.distances, "x", "distance" );
			if ( settings.showIntermediateResults )
				Plots.plot( centroidsParameters.axisCoordinates, centroidsParameters.numVoxels, "x", "numVoxels" );
			if ( settings.showIntermediateResults )
				show( yawAndOrientationAlignedMask, "yaw and orientation aligned mask", centroidsParameters.centroids, registrationCalibration, false );

			final AffineTransform3D rollTransform = computeCentroidBasedRollTransform( centroidsParameters, settings );

			registration = registration.preConcatenate( rollTransform );

		}
		else if ( settings.rollAngleComputationMethod.equals( ShavenBabyRegistrationSettings.PROJECTION_SHAPE ) )
		{

			final RandomAccessibleInterval< UnsignedIntType > intMask = Converters.convert( yawAndOrientationAlignedMask, ( i, o ) -> o.set( i.getRealDouble() > 0 ? 1000 : 0 ), new UnsignedIntType() );

			final AffineTransform3D intensityBasedRollTransform = computeIntensityBasedRollTransform(
					intMask,
					intMask.min( X ) * settings.registrationResolution,
					intMask.max( X ) * settings.registrationResolution,
					12.0 );

			registration = registration.preConcatenate( intensityBasedRollTransform );
		}
		return registration;
	}

	public RandomAccessibleInterval< BitType > close(
			RandomAccessibleInterval< BitType > mask,
			int closingRadius )
	{
		// TODO: Bug(?!) in imglib2 Closing.close makes this necessary
		RandomAccessibleInterval< BitType > closed = ArrayImgs.bits( Intervals.dimensionsAsLongArray( mask ) );
		final RandomAccessibleInterval< BitType > enlargedMask = Utils.getEnlargedRai2( mask, closingRadius );
		final RandomAccessibleInterval< BitType > enlargedClosed = Utils.getEnlargedRai2( closed, closingRadius );

		if ( closingRadius > 0 )
		{
			Utils.log( "Morphological closing...");
			Shape closingShape = new HyperSphereShape( closingRadius );
			Closing.close( Views.extendZero( enlargedMask ), Views.iterable( enlargedClosed ), closingShape, 1 );
		}

		return Views.interval( enlargedClosed, mask );
	}

	public < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< BitType > createMask( RandomAccessibleInterval< T > downscaled, double threshold )
	{
		Utils.log( "Creating mask...");

		RandomAccessibleInterval< BitType > mask = Converters.convert( downscaled, ( i, o ) -> o.set( i.getRealDouble() > threshold ? true : false ), new BitType() );

		return mask;
	}

	public < T extends RealType< T > & NativeType< T > > double getThreshold( RandomAccessibleInterval< T > downscaled )
	{

		double threshold = 0;

		if ( settings.thresholdModality.equals( ShavenBabyRegistrationSettings.HUANG_AUTO_THRESHOLD ) )
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

	public ImgLabeling< Integer, IntType > createWatershedSeeds( double[] registrationCalibration, RandomAccessibleInterval< DoubleType > distance, RandomAccessibleInterval< BitType > mask )
	{
		Utils.log( "Seeds for watershed...");

		double globalDistanceThreshold = Math.pow( settings.watershedSeedsGlobalDistanceThreshold / settings.registrationResolution, 2 );
		double localMaximaDistanceThreshold = Math.pow( settings.watershedSeedsLocalMaximaDistanceThreshold / settings.registrationResolution, 2 );

		// TODO: remove local maxima detection
		final RandomAccessibleInterval< BitType >  seeds = Algorithms.createCenterAndBoundarySeeds(
				distance,
				new HyperSphereShape( 1 ),
				globalDistanceThreshold,
				localMaximaDistanceThreshold );

		final ImgLabeling< Integer, IntType > seedsLabelImg = Utils.asImgLabeling( seeds );

		if ( settings.showIntermediateResults ) show( Utils.asIntImg( seedsLabelImg ), "watershed seeds", null, registrationCalibration, false );
		return seedsLabelImg;
	}

	public AffineTransform3D computeFlippingTransform( RandomAccessibleInterval yawAlignedMask, RandomAccessibleInterval yawAlignedIntensities, double calibration )
	{
		final CoordinatesAndValues coordinatesAndValues = Utils.computeAverageIntensitiesAlongAxis( yawAlignedIntensities, yawAlignedMask, X, calibration );

		if ( settings.showIntermediateResults ) Plots.plot( coordinatesAndValues.coordinates, coordinatesAndValues.values, "x", "average intensity" );

		double maxLoc = Utils.computeMaxLoc( coordinatesAndValues.coordinates, coordinatesAndValues.values, null );

		AffineTransform3D affineTransform3D = new AffineTransform3D();

		if ( maxLoc < 0 ) affineTransform3D.rotate( Z, toRadians( 180.0D ) );

		return affineTransform3D;
	}

	public < T extends RealType< T > & NativeType< T > >
	AffineTransform3D computeIntensityBasedRollTransform(
			RandomAccessibleInterval rai,
			double xMin,
			double xMax,
			double blurSigma )
	{
		final RandomAccessibleInterval< T > longAxisProjection = Utils.createAverageProjectionAlongAxis(
				rai,
				X,
				xMin,
				xMax,
				settings.registrationResolution );

		// if ( settings.showIntermediateResults ) show( longAxisProjection, "amnioserosa projection", null, Utils.get3dDoubleArray( calibration ) , false );

		final RandomAccessibleInterval< T > blurred = Utils.createBlurredRai(
				longAxisProjection,
				blurSigma,
				settings.registrationResolution );

		final Point maximum = Algorithms.findMaximumLocation( blurred, Utils.get2dDoubleArray( settings.registrationResolution ));
		final List< RealPoint > realPoints = Utils.asRealPointList( maximum );
		realPoints.add( new RealPoint( new double[]{ 0, 0 } ) );

		if ( settings.showIntermediateResults ) show( blurred, "perpendicular projection - blurred ", realPoints, Utils.get2dDoubleArray( settings.registrationResolution ), false );

		final AffineTransform3D xAxisRollTransform = createXAxisRollTransform( maximum );
		xAxisRollTransform.rotate( X, Math.PI );

		return xAxisRollTransform;
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

	public static AffineTransform3D computeCentroidBasedRollTransform( CentroidsParameters centroidsParameters, ShavenBabyRegistrationSettings settings )
	{
		final double rollAngle = computeRollAngle( centroidsParameters, settings.rollAngleMinDistanceToAxis, settings.rollAngleMinDistanceToCenter, settings.rollAngleMaxDistanceToCenter );

		Utils.log( "Roll angle " + rollAngle );

		AffineTransform3D rollTransform = new AffineTransform3D();

		rollTransform.rotate( X, toRadians( rollAngle ) );

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

		return -1 * medianAngle;
	}

	public static ArrayList< RealPoint > origin()
	{
		final ArrayList< RealPoint > origin = new ArrayList<>();
		origin.add( new RealPoint( new double[]{ 0, 0, 0 } ) );
		return origin;
	}

	private AffineTransform3D createFinalTransform( double[] inputCalibration, AffineTransform3D registration, double[] registrationCalibration )
	{
		final AffineTransform3D transform =
				Transforms.getScalingTransform( inputCalibration, settings.registrationResolution )
				.preConcatenate( registration )
				.preConcatenate( Transforms.getScalingTransform( registrationCalibration, settings.outputResolution ) );

		return transform;
	}

	private Img< BitType > createMaskFromLabelRegion( LabelRegion< Integer > centralObjectRegion, long[] dimensions )
	{
		final Img< BitType > centralObjectImg = ArrayImgs.bits( dimensions );

		final Cursor< Void > regionCursor = centralObjectRegion.cursor();
		final net.imglib2.RandomAccess< BitType > access = centralObjectImg.randomAccess();
		while ( regionCursor.hasNext() )
		{
			regionCursor.fwd();
			access.setPosition( regionCursor );
			access.get().set( true );
		}
		return centralObjectImg;
	}

	private Img< UnsignedByteType > createUnsignedByteTypeMaskFromLabelRegion( LabelRegion< Integer > centralObjectRegion, long[] dimensions )
	{
		final Img< UnsignedByteType > centralObjectImg = ArrayImgs.unsignedBytes( dimensions );

		final Cursor< Void > regionCursor = centralObjectRegion.cursor();
		final net.imglib2.RandomAccess< UnsignedByteType > access = centralObjectImg.randomAccess();
		while ( regionCursor.hasNext() )
		{
			regionCursor.fwd();
			access.setPosition( regionCursor );
			access.get().set( 255 );
		}
		return centralObjectImg;
	}


	private LabelRegion< Integer > getCentralObjectLabelRegion( ImgLabeling< Integer, IntType > labeling )
	{
		int centralLabel = getCentralLabel( labeling );

		if ( centralLabel == -1 )
		{
			return null;
		}

		final LabelRegions< Integer > labelRegions = new LabelRegions<>( labeling );

		return labelRegions.getLabelRegion( centralLabel );
	}

	private static int getCentralLabel( ImgLabeling< Integer, IntType > labeling )
	{
		final net.imglib2.RandomAccess< LabelingType< Integer > > labelingRandomAccess = labeling.randomAccess();
		for ( int d : XYZ ) labelingRandomAccess.setPosition( labeling.dimension( d ) / 2, d );
		int centralIndex = labelingRandomAccess.get().getIndex().getInteger();

		if ( centralIndex > 0 )
		{
			return labeling.getMapping().labelsAtIndex( centralIndex ).iterator().next();
		}
		else
		{
			return -1;
		}
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
