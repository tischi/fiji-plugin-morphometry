package de.embl.cba.morphometry.spindle;

import bdv.util.Bdv;
import bdv.util.BdvFunctions;
import bdv.util.BdvOptions;
import de.embl.cba.morphometry.*;
import de.embl.cba.morphometry.geometry.CentroidsParameters;
import de.embl.cba.morphometry.geometry.CoordinatesAndValues;
import de.embl.cba.morphometry.geometry.EllipsoidParameters;
import de.embl.cba.morphometry.geometry.Ellipsoids;
import ij.IJ;
import net.imagej.ops.OpService;
import net.imglib2.*;
import net.imglib2.algorithm.morphology.Closing;
import net.imglib2.algorithm.neighborhood.HyperSphereShape;
import net.imglib2.algorithm.neighborhood.Neighborhood;
import net.imglib2.algorithm.neighborhood.Shape;
import net.imglib2.converter.Converters;
import net.imglib2.histogram.Histogram1d;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.interpolation.randomaccess.NearestNeighborInterpolatorFactory;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.roi.labeling.LabelRegion;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.util.Intervals;
import net.imglib2.util.LinAlgHelpers;
import net.imglib2.view.Views;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import static de.embl.cba.morphometry.Algorithms.angleInRadians;
import static de.embl.cba.morphometry.Constants.*;
import static de.embl.cba.morphometry.Transforms.getScalingFactors;
import static de.embl.cba.morphometry.viewing.BdvImageViewer.show;
import static java.lang.Math.toRadians;


public class SpindleMorphometry  < T extends RealType< T > & NativeType< T > >
{

	final SpindleMorphometrySettings settings;
	final OpService opService;

	public SpindleMorphometry( SpindleMorphometrySettings settings, OpService opService )
	{
		this.settings = settings;
		this.opService = opService;
	}

	public void run()
	{

		/**
		 *  Make isotropic
		 */

		Utils.log( "Create isotropic image..." );

		final double[] workingCalibration = Utils.get3dDoubleArray( settings.workingVoxelSize );

		final RandomAccessibleInterval< T > dapi = Algorithms.createIsotropicArrayImg( settings.dapi, getScalingFactors( settings.inputCalibration, settings.workingVoxelSize ) );
		final RandomAccessibleInterval< T > tubulin = Algorithms.createIsotropicArrayImg( settings.tubulin, getScalingFactors( settings.inputCalibration, settings.workingVoxelSize ) );

		if ( settings.showIntermediateResults ) show( dapi, "image isotropic resolution", null, workingCalibration, false );
		if ( settings.showIntermediateResults ) show( tubulin, "tubulin isotropic resolution", null, workingCalibration, false );


		/**
		 *  Compute offset and threshold
		 */

		final RandomAccessibleInterval< T > dapi3um = Algorithms.createIsotropicArrayImg( settings.dapi, getScalingFactors( settings.inputCalibration, 3.0 ) );
		final double maximumValue = Algorithms.getMaximumValue( dapi3um );
		double threshold = maximumValue / 2.0;

		Utils.log( "Dapi threshold: " + threshold );


		/**
		 * Create mask
		 */

		RandomAccessibleInterval< BitType > mask = createMask( dapi, threshold );

		if ( settings.showIntermediateResults ) show( mask, "mask", null, workingCalibration, false );

		/**
		 * Morphological closing
		 */

//		RandomAccessibleInterval< BitType > closed = createClosedImage( mask );
		RandomAccessibleInterval< BitType > closed =  mask ;

//		if ( settings.showIntermediateResults ) show( closed, "closed", null, workingCalibration, false );


		/**
		 * Extract metaphase plate object
		 */

		Utils.log( "Extracting metaphase plate object..." );

		final ImgLabeling< Integer, IntType > labelImg = Algorithms.createImgLabeling( closed );

		final LabelRegion< Integer > largestObject = Algorithms.getLargestObject( labelImg );

		final Img< BitType > dapiMask = Algorithms.createMaskFromLabelRegion( largestObject, Intervals.dimensionsAsLongArray( labelImg ) );

		if ( settings.showIntermediateResults ) show( dapiMask, "meta-phase object", null, workingCalibration, false );

		/**
		 * Compute ellipsoid alignment
		 */

		// TODO: maybe fit intensities rather than mask, but how to mask the intensities of the this cell

		Utils.log( "Determining metaphase plate axes..." );

		final EllipsoidParameters ellipsoidParameters = Ellipsoids.computeParametersFromBinaryImage( dapiMask );

		Utils.log( "Creating aligned images..." );

		final AffineTransform3D alignmentTransform = Ellipsoids.createAlignmentTransform( ellipsoidParameters );
		final RandomAccessibleInterval aligendTubulin = Utils.copyAsArrayImg( Transforms.createTransformedView( tubulin, alignmentTransform, new NearestNeighborInterpolatorFactory() ) );
		final RandomAccessibleInterval alignedDapi = Utils.copyAsArrayImg( Transforms.createTransformedView( dapi, alignmentTransform ) );

		if ( settings.showIntermediateResults ) show( alignedDapi, "aligned image", null, workingCalibration, false );

		RandomAccessibleInterval< T > interestPoints = Utils.copyAsArrayImg( dapi );
		Algorithms.setValues( interestPoints, 0.0 );

		final double[] dnaCenter = new double[ 3 ];
		alignmentTransform.inverse().apply( new double[]{0,0,0}, dnaCenter );
		drawPoint( interestPoints, dnaCenter, settings.interestPointsRadius, settings.workingVoxelSize );


		/**
		 * Compute metaphase plate width
		 */

		Utils.log( "Determining widths..." );

		final CoordinatesAndValues dapiProfile = Utils.computeAverageIntensitiesAlongAxis( alignedDapi, settings.maxShortAxisDist, 2, settings.workingVoxelSize );
		if ( settings.showIntermediateResults ) Plots.plot( dapiProfile.coordinates, dapiProfile.values, "distance to center", "image intensity" );

		/**
		 * Compute spindle length
		 */

		final CoordinatesAndValues tubulinProfile = Utils.computeMaximumIntensitiesAlongAxis( aligendTubulin, settings.maxShortAxisDist, 2, settings.workingVoxelSize );
		if ( settings.showIntermediateResults ) Plots.plot( tubulinProfile.coordinates, tubulinProfile.values, "distance to center", "tubulin intensity" );

		final ArrayList< Double > tubulinProfileAbsoluteDerivative = Algorithms.computeAbsoluteDerivatives( tubulinProfile.values, ( int ) ( 1.0 / settings.workingVoxelSize ) );
		if ( settings.showIntermediateResults ) Plots.plot( tubulinProfile.coordinates, tubulinProfileAbsoluteDerivative, "distance to center", "tubulin intensity absolute derivative" );

		double[] maxLocs = getLeftAndRightMaxLocs( tubulinProfile, tubulinProfileAbsoluteDerivative );

		addSpindlePolesAsInterestPoints( alignmentTransform, interestPoints, maxLocs );

		final ArrayList< RealPoint > spindleLengthPoints = new ArrayList<>(  );
		spindleLengthPoints.add( new RealPoint( new double[]{ 0.0, 0.0, maxLocs[ 0 ] } ));
		spindleLengthPoints.add( new RealPoint( new double[]{ 0.0, 0.0, maxLocs[ 1 ] } ));

		if ( settings.showIntermediateResults ) show( aligendTubulin, "aligned tubulin", spindleLengthPoints, workingCalibration, false );


		AffineTransform3D rotation = alignmentTransform.copy();
		rotation.setTranslation( new double[ 3 ] );

		final double[] zAxis = { 0, 0, 1 };
		final double[] xAxis = { 1, 0, 0 };

		final double[] shortDnaAxisOCS = new double[ 3 ];

		rotation.inverse().apply( zAxis, shortDnaAxisOCS );

		final double[] sideViewDirectionOCS = new double[ 3 ];
		LinAlgHelpers.cross( shortDnaAxisOCS, zAxis, sideViewDirectionOCS );

		final double dot = LinAlgHelpers.dot( xAxis, sideViewDirectionOCS );
		final double angleInRadians = angleInRadians( xAxis, sideViewDirectionOCS );

		rotation = new AffineTransform3D();
		rotation.rotate( 2,  - angleInRadians );

		final RandomAccessibleInterval transformedDapiView = Transforms.createTransformedView( dapi, rotation );
		final RandomAccessibleInterval transformedTubulinView = Transforms.createTransformedView( tubulin, rotation );
		final RandomAccessibleInterval transformedInterestPointView = Transforms.createTransformedView( interestPoints, rotation );

		Bdv bdv = null;
		if ( settings.showIntermediateResults ) BdvFunctions.show( interestPoints, "" ).getBdvHandle();
		if ( settings.showIntermediateResults ) bdv = BdvFunctions.show( transformedDapiView, "" ).getBdvHandle();
		if ( settings.showIntermediateResults ) BdvFunctions.show( transformedInterestPointView, "", BdvOptions.options().addTo( bdv ) );

		Utils.log( "Saving result images..." );
		saveMaximumProjections( transformedDapiView, "image" );
		saveMaximumProjections( transformedTubulinView, "tubulin" );
		saveMaximumProjections( transformedInterestPointView, "points" );


//
//		if ( settings.showIntermediateResults )  show( transformedView, "side view image", null, workingCalibration, false );
//
//		int a = 1;


	}

	public void saveMaximumProjections( RandomAccessibleInterval transformedDapiView, String name )
	{
		for ( int d = 0; d < 3; ++d )
		{
			final String title = name + "_view" + d;

			final String path = settings.outputDirectory.getAbsolutePath()
					+ File.separator
					+ settings.inputDataSetName
					+ File.separator
					+ title + ".tiff";

			new File(path).getParentFile().mkdirs();

			IJ.saveAsTiff( ImageJFunctions.wrap( new Projection( transformedDapiView, d ).maximum(), title ), path );
		}
	}

	public void addSpindlePolesAsInterestPoints( AffineTransform3D alignmentTransform, RandomAccessibleInterval< T > interestPoints, double[] maxLocs )
	{
		double[] spindlePole = new double[ 3 ];
		spindlePole = new double[]{ 0.0, 0.0, maxLocs[ 0 ] };
		convertToOriginalImagePixelUnits( alignmentTransform, spindlePole );
		drawPoint( interestPoints, spindlePole, settings.interestPointsRadius, settings.workingVoxelSize );

		spindlePole = new double[]{ 0.0, 0.0, maxLocs[ 1 ] };
		convertToOriginalImagePixelUnits( alignmentTransform, spindlePole );
		drawPoint( interestPoints, spindlePole, settings.interestPointsRadius, settings.workingVoxelSize );
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

		RandomAccessibleInterval< BitType > mask = Converters.convert( downscaled, ( i, o ) -> o.set( i.getRealDouble() > threshold ? true : false ), new BitType() );

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

		final RandomAccessibleInterval< BitType >  seeds = Algorithms.createSeeds(
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
