package de.embl.cba.morphometry.spindle;

import de.embl.cba.morphometry.*;
import de.embl.cba.morphometry.geometry.CoordinatesAndValues;
import de.embl.cba.morphometry.geometry.CurveAnalysis;
import de.embl.cba.morphometry.geometry.ellipsoids.EllipsoidVectors;
import de.embl.cba.morphometry.geometry.ellipsoids.Ellipsoids3DImageSuite;
import de.embl.cba.morphometry.measurements.Measurements;
import de.embl.cba.morphometry.regions.Regions;
import de.embl.cba.transforms.utils.Transforms;
import ij.IJ;
import ij.ImagePlus;
import ij.measure.Calibration;
import net.imagej.ops.OpService;
import net.imglib2.*;
import net.imglib2.algorithm.labeling.ConnectedComponents;
import net.imglib2.algorithm.neighborhood.HyperSphereShape;
import net.imglib2.algorithm.neighborhood.Neighborhood;
import net.imglib2.algorithm.neighborhood.Shape;
import net.imglib2.converter.Converters;
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
import net.imglib2.util.Intervals;
import net.imglib2.util.LinAlgHelpers;
import net.imglib2.view.Views;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import static de.embl.cba.morphometry.Angles.angleOfSpindleAxisToXAxisInRadians;
import static de.embl.cba.morphometry.viewing.BdvViewer.show;
import static de.embl.cba.transforms.utils.Scalings.createRescaledArrayImg;
import static de.embl.cba.transforms.utils.Transforms.getScalingFactors;


public class SpindleMorphometry  < T extends RealType< T > & NativeType< T > >
{

	public static final String SPINDLE_LENGTH = "Spindle_Length";
	public static final String SPINDLE_WIDTH = "Spindle_Width";
	public static final String DNA_WIDTH = "DNA_Width";
	public static final String DNA_SPINDLE_CENTER_DISTANCE = "Dna_Center_To_Spindle_Center_Distance";
	public static final String SPINDLE_AXIS_TO_COVERSLIP_PLANE_ANGLE_DEGREES = "Spindle_Axis_To_Coverslip_Plane_Angle_Degrees";
	public static final String LENGTH_UNIT = "um";
	public static final String DNA_VOLUME = "Dna_Volume";
	public static final String VOLUME_UNIT = "um3";
	public static final String DNA_LENGTH = "Dna_Length";
	public static final int ALIGNED_DNA_AXIS = 2;
	public static final String SEP = "_";

	final SpindleMorphometrySettings settings;
	final OpService opService;

	private HashMap< Integer, Map< String, Object > > objectMeasurements;
	private RandomAccessibleInterval< BitType > mask;
	private RandomAccessibleInterval transformedDapiView;
	private RandomAccessibleInterval transformedTubulinView;
	private RandomAccessibleInterval transformedInterestPointView;

	public SpindleMorphometry( SpindleMorphometrySettings settings, OpService opService )
	{
		this.settings = settings;
		this.opService = opService;
	}

	public void run()
	{
		/**
		 *  Initialise measurements
		 */

		objectMeasurements = new HashMap<>();

		/**
		 *  Make isotropic
		 */

		Logger.log( "Create isotropic image..." );

		final double[] workingCalibration = Utils.as3dDoubleArray( settings.workingVoxelSize );

		final RandomAccessibleInterval< T > dna = createRescaledArrayImg( settings.dapiImage, getScalingFactors( settings.inputCalibration, settings.workingVoxelSize ) );
		final RandomAccessibleInterval< T > tubulin = createRescaledArrayImg( settings.tubulinImage, getScalingFactors( settings.inputCalibration, settings.workingVoxelSize ) );

		if ( settings.showIntermediateResults ) show( dna, "dapi isotropic resolution", null, workingCalibration, false );
		if ( settings.showIntermediateResults ) show( tubulin, "tubulin isotropic resolution", null, workingCalibration, false );


		/**
		 *  Compute offset and threshold
		 */

		final RandomAccessibleInterval< T > dnaDownscaledToSpindleWidth = createRescaledArrayImg(
				settings.dapiImage,
				getScalingFactors( settings.inputCalibration, 3.0 ) );
		final double maximumValue = Algorithms.getMaximumValue( dnaDownscaledToSpindleWidth );
		double threshold = maximumValue / 2.0;

		Logger.log( "Dapi threshold: " + threshold );


		/**
		 * Create mask
		 */

		RandomAccessibleInterval< BitType > dnaMask = createMask( dna, threshold );

		if ( settings.showIntermediateResults ) show( dnaMask, "dapi mask", null, workingCalibration, false );


		/**
		 * Extract metaphase plate object
		 */

		Logger.log( "Extracting metaphase plate object..." );

		Regions.removeSmallRegionsInMask( dnaMask, settings.minimalMetaphasePlateVolumeInCalibratedUnits, settings.workingVoxelSize );

		final ImgLabeling< Integer, IntType > labelImg = Utils.asImgLabeling( dnaMask, ConnectedComponents.StructuringElement.FOUR_CONNECTED );

		final long radius = (long) ( settings.centralObjectRegionToleranceInCalibratedUnits / settings.workingVoxelSize );
		final LabelRegion< Integer > metaphasePlateRegion = Regions.getCentralRegion( labelImg, radius );

		final Img< BitType > metaphasePlateMask = Algorithms.createMaskFromLabelRegion( metaphasePlateRegion, Intervals.dimensionsAsLongArray( labelImg ) );

		if ( settings.showIntermediateResults ) show( metaphasePlateMask, "meta-phase object", null, workingCalibration, false );

		/**
		 * Compute meta-phase plate volume
		 */

		Measurements.addMeasurement(
				objectMeasurements,
				0,
				DNA_VOLUME + SEP + VOLUME_UNIT,
				metaphasePlateRegion.size() * Math.pow( settings.workingVoxelSize, 3 ) );

		/**
		 * Morphological filtering
		 *
		 * - it appears that the principal axes determination is not working robustly
		 * if the meta-phase plate is too thick
		 */

		final RandomAccessibleInterval< BitType > processedMetaPhasePlate = createProcessedMetaPhasePlate( dnaMask, metaphasePlateMask );

		if ( settings.showIntermediateResults ) show( processedMetaPhasePlate, "processed metaphase plate", null, workingCalibration, false );

		/**
		 * Compute ellipsoid alignment
		 */

		Logger.log( "Determining meta-phase plate axes..." );

		//final EllipsoidMLJ ellipsoidParameters = EllipsoidsMLJ.computeParametersFromBinaryImage( processedMetaPhasePlate );

		final EllipsoidVectors ellipsoidVectors = Ellipsoids3DImageSuite.fitEllipsoid( Utils.asImagePlus( processedMetaPhasePlate, "" ) );

		Logger.log( "Creating aligned images..." );

		//final AffineTransform3D alignmentTransform = EllipsoidsMLJ.createAlignmentTransform( ellipsoidParameters );

		final AffineTransform3D alignmentTransform = Ellipsoids3DImageSuite.createAlignmentTransform( ellipsoidVectors );
		final RandomAccessibleInterval tubulinAlignedAlongShortestAxisOfDNA = Utils.copyAsArrayImg( Transforms.createTransformedView( tubulin, alignmentTransform, new NearestNeighborInterpolatorFactory() ) );
		final RandomAccessibleInterval alignedDNA = Utils.copyAsArrayImg( Transforms.createTransformedView( dna, alignmentTransform ) );
		final RandomAccessibleInterval alignedProcessedMetaphasePlate = Utils.copyAsArrayImg( Transforms.createTransformedView( processedMetaPhasePlate, alignmentTransform ) );

		if ( settings.showIntermediateResults ) show( alignedDNA, "aligned dapi", Transforms.origin(), workingCalibration, false );
		if ( settings.showIntermediateResults ) show( alignedProcessedMetaphasePlate, "aligned processed meta-phase plate", Transforms.origin(), workingCalibration, false );
//		if ( settings.showIntermediateResults ) Viewer3D.show3D( alignedProcessedMetaphasePlate );

		/**
		 * Compute metaphase plate length
		 */

		Logger.log( "Measuring meta-phase plate length..." );

		final CoordinatesAndValues dnaProfile = Utils.computeAverageIntensitiesAlongAxis( alignedDNA, settings.maxShortAxisDist, ALIGNED_DNA_AXIS, settings.workingVoxelSize );
		if ( settings.showIntermediateResults ) Plots.plot( dnaProfile.coordinates, dnaProfile.values, "distance to center", "dapi intensity along shortest axis" );

		// TODO: replace also by gradient and maximum?
		final double dnaLength = CurveAnalysis.computeFWHM( dnaProfile );

		Measurements.addMeasurement(
				objectMeasurements,
				0,
				DNA_LENGTH + SEP + LENGTH_UNIT,
				dnaLength);


		/**
		 * Measure metaphase plate width
		 */

		Logger.log( "Measuring meta-phase plate width..." );

		Projection projection = new Projection<>( alignedDNA, 2 );

		final RandomAccessibleInterval< T > dnaProjection = projection.maximum();

		final double dnaWidth = measureWidthByRadialProfile( dnaProjection, "dna width" );

		Measurements.addMeasurement(
				objectMeasurements,
				0,
				DNA_WIDTH + SEP + LENGTH_UNIT,
				dnaWidth );

		/**
		 * Compute spindle end-points along shortest DNA axis
		 */

		final CoordinatesAndValues tubulinProfile = Utils.computeMaximumIntensitiesAlongAxis( tubulinAlignedAlongShortestAxisOfDNA, settings.maxShortAxisDist, ALIGNED_DNA_AXIS, settings.workingVoxelSize );
		if ( settings.showIntermediateResults ) Plots.plot( tubulinProfile.coordinates, tubulinProfile.values, "center distance [um]", "tubulin max along shortest DNA axis" );

		final CoordinatesAndValues tubulinProfileDerivative =
				CurveAnalysis.computeDerivatives(
						tubulinProfile,
						(int) Math.ceil( settings.derivativeDelta / settings.workingVoxelSize ) );

		if ( settings.showIntermediateResults ) Plots.plot( tubulinProfileDerivative.coordinates, tubulinProfileDerivative.values, "distance to center", "d/dx tubulin max along shortest DNA axis" );

		double[] dnaAxisBasedSpindlePoleCoordinates = getLeftMaxAndRightMinLoc( tubulinProfileDerivative.coordinates, tubulinProfileDerivative.values );

		final ArrayList< RealPoint > spindleLengthPoints = new ArrayList<>(  );

		spindleLengthPoints.add( new RealPoint( new double[]{ 0.0, 0.0, 0 } ));
		spindleLengthPoints.add( new RealPoint( new double[]{ 0.0, 0.0, dnaAxisBasedSpindlePoleCoordinates[ 0 ] } ));
		spindleLengthPoints.add( new RealPoint( new double[]{ 0.0, 0.0, dnaAxisBasedSpindlePoleCoordinates[ 1 ] } ));


		/**
		 * Refine spindle pole locations, by finding off axis maximum
		 * in a plane
		 */

		final ArrayList< double[] > spindlePoles = new ArrayList<>();

		for ( int pole = 0; pole < 2; pole++ )
		{
			RandomAccessibleInterval< T > tubulinIntensitySlice =
					Views.hyperSlice(
							tubulinAlignedAlongShortestAxisOfDNA,
							2,
							(long) ( dnaAxisBasedSpindlePoleCoordinates[ pole ] / settings.workingVoxelSize ) );

			final double[] xy = Utils.computeMaximumLocation( tubulinIntensitySlice, settings.maxShortAxisDist / settings.workingVoxelSize );
			spindlePoles.add( new double[] {
							xy[ 0 ] * settings.workingVoxelSize,
							xy[ 1 ] * settings.workingVoxelSize,
							dnaAxisBasedSpindlePoleCoordinates[ pole ] } );
		}

		/**
		 * Add spindle length (pole to pole distance) to measurements
		 */

		final double distance = LinAlgHelpers.distance( spindlePoles.get( 0 ), spindlePoles.get( 1 ) );

		Measurements.addMeasurement(
				objectMeasurements,
				0,
				SPINDLE_LENGTH + SEP + LENGTH_UNIT,
				distance );


		/**
		 * Determine spindle center
		 */

		final double[] poleToPoleVector = new double[ 3 ];
		LinAlgHelpers.subtract( spindlePoles.get( 0 ), spindlePoles.get( 1 ), poleToPoleVector );
		for ( int i = 0; i < 3; i++ )
		{
			poleToPoleVector[ i ] *= 0.5;
		}
		final double[] spindleCenter = new double[ 3 ];
		LinAlgHelpers.add( spindlePoles.get( 1 ), poleToPoleVector, spindleCenter );

		final double dnaCenterToSpindleCenterDistance = LinAlgHelpers.distance( new double[]{ 0, 0, 0}, spindleCenter) * settings.workingVoxelSize;

		Measurements.addMeasurement(
				objectMeasurements,
				0,
				DNA_SPINDLE_CENTER_DISTANCE + SEP + LENGTH_UNIT,
				dnaCenterToSpindleCenterDistance );


		if ( settings.showIntermediateResults )
		{
			final ArrayList< RealPoint > realPoints = new ArrayList<>();
			realPoints.add( new RealPoint( spindlePoles.get( 0 ) ) );
			realPoints.add( new RealPoint( spindlePoles.get( 1 ) ) );
			realPoints.add( new RealPoint( spindleCenter ) );
			show(
					tubulinAlignedAlongShortestAxisOfDNA, "tubulin aligned along shortest DNA axis",
					realPoints,
					workingCalibration,
					false );
		}

		/**
		 * Measure angle between spindle axis and coverslip plane
		 */

		final double[] poleToPoleVectorInCBCS = transformToCoverslipBasedCoordinateSystem( alignmentTransform, poleToPoleVector );
		final double angleSpindleAxisToCoverslipPlaneInDegrees = 90.0 - Math.abs( 180.0 / Math.PI * Transforms.getAngle( new double[]{ 0, 0, 1 }, poleToPoleVectorInCBCS ) );

		Measurements.addMeasurement(
				objectMeasurements,
				0,
				SPINDLE_AXIS_TO_COVERSLIP_PLANE_ANGLE_DEGREES,
				angleSpindleAxisToCoverslipPlaneInDegrees );


		/**
		 * Measure spindle width
		 */

		// Realign spindle along pole to pole axis
		//
		final double[] poleToPoleAxis = new double[ 3 ];
		LinAlgHelpers.subtract( spindlePoles.get( 0 ), spindlePoles.get( 1 ), poleToPoleAxis );
		LinAlgHelpers.normalize( poleToPoleAxis );

		AffineTransform3D spindleTransform = new AffineTransform3D();
		spindleTransform.translate( spindleCenter );
		spindleTransform = spindleTransform.inverse();
		AffineTransform3D poleToPoleAxisRotation = Transforms.getRotationTransform3D( new double[]{ 0, 0, 1 }, poleToPoleAxis );
		spindleTransform = spindleTransform.preConcatenate( poleToPoleAxisRotation );

		final RandomAccessibleInterval tubulinAlignedAlongSpindlePoleToPoleAxis
				= Transforms.createTransformedView( tubulinAlignedAlongShortestAxisOfDNA, spindleTransform );

		final double[] newPole00 = new double[ 3 ];
		spindleTransform.apply( spindlePoles.get( 0 ), newPole00 );
		final double[] newPole01 = new double[ 3 ];
		spindleTransform.apply( spindlePoles.get( 1 ), newPole01 );
		final double[] newCenter = new double[ 3 ];
		spindleTransform.apply( spindleCenter, newCenter );

		if ( settings.showIntermediateResults )
		{
			final ArrayList< RealPoint > realPoints = new ArrayList<>();
			realPoints.add( new RealPoint( newPole00 ) );
			realPoints.add( new RealPoint( newPole01 ) );
			realPoints.add( new RealPoint( newCenter ) );
			show(
					tubulinAlignedAlongSpindlePoleToPoleAxis, "tubulin aligned pole to pole",
					realPoints,
					workingCalibration,
					false );
		}

		// Measure width in projection along spindle axis
		//
		projection = new Projection<>( tubulinAlignedAlongSpindlePoleToPoleAxis, 2 );

		final RandomAccessibleInterval< T > spindleProjection = projection.maximum();

		final double spindleWidth = measureWidthByRadialProfile( spindleProjection, "tubulin perpendicular to spindle axis" );

		Measurements.addMeasurement(
				objectMeasurements,
				0,
				SPINDLE_WIDTH + SEP + LENGTH_UNIT,
				spindleWidth );

//		final Histogram1d< T > spindleProjectionHistogram = opService.image().histogram( Views.iterable( spindleProjection ) );
//		final double spindleProjectionThreshold = opService.threshold().yen( spindleProjectionHistogram ).getRealDouble();
//		RandomAccessibleInterval< BitType > spindleProjectionMask =
//				Algorithms.createMask( spindleProjection, spindleProjectionThreshold );
//
//		final long spindleProjectionAreaInPixels = Measurements.measureSize( spindleProjectionMask );
//		double spindleWidth = 2 * Math.sqrt( spindleProjectionAreaInPixels / Math.PI ) * settings.workingVoxelSize;
//
//		Measurements.addMeasurement(
//				objectMeasurements,
//				0,
//				SPINDLE_WIDTH + SEP + LENGTH_UNIT,
//				spindleWidth );
//
//		ImageJFunctions.show( spindleProjection, "spindleProjection" );
//		ImageJFunctions.show( spindleProjectionMask, "spindleProjectionMask" );


		/**
		 * Create interest point image in coverslip coordinate system
		 */

		RandomAccessibleInterval< T > interestPointsImage = createInterestPointImage( dna, alignmentTransform, spindlePoles );

		/**
		 * Create output image
		 */

		createOutputImage( dna, tubulin, alignmentTransform, interestPointsImage );


	}

	public void createOutputImage( RandomAccessibleInterval< T > dna, RandomAccessibleInterval< T > tubulin, AffineTransform3D alignmentTransform, RandomAccessibleInterval< T > interestPointsImage )
	{
		AffineTransform3D rotation = alignmentTransform.copy(); // "rotation" rotates the data such that the spindle axis is the zAxis
		rotation.setTranslation( new double[ 3 ] );

		final double[] zAxis = { 0, 0, 1 };

		// compute the spindle axis in the coordinate system of the input data
		final double[] spindleAxis = new double[ 3 ];
		rotation.inverse().apply( zAxis, spindleAxis );

		// set z-component of spindleAxis to zero, because we want to rotate parallel to coverslip
		spindleAxis[ ALIGNED_DNA_AXIS ] = 0;

		// compute rotation vector
		// note: we do not know in which direction of the spindle the spindleAxis vector points,
		// but we think it does not matter for alignment to the x-axis
		final double angleBetweenXAxisAndSpindleWithVector = angleOfSpindleAxisToXAxisInRadians( spindleAxis );

		AffineTransform3D rotationTransform = new AffineTransform3D();
		rotationTransform.rotate( ALIGNED_DNA_AXIS,  angleBetweenXAxisAndSpindleWithVector );

		Logger.log( "Rotating input data around z-axis by [degrees]: " + 180 / Math.PI * angleBetweenXAxisAndSpindleWithVector );

		transformedDapiView = Transforms.createTransformedView( dna, rotationTransform );
		transformedTubulinView = Transforms.createTransformedView( tubulin, rotationTransform );
		transformedInterestPointView = Transforms.createTransformedView( interestPointsImage, rotationTransform, new NearestNeighborInterpolatorFactory() );
	}

	public RandomAccessibleInterval< T > createInterestPointImage(
			RandomAccessibleInterval< T > dna,
			AffineTransform3D alignmentTransform,
			ArrayList< double[] > spindlePoles )
	{
		RandomAccessibleInterval< T > interestPointsImage = Utils.copyAsArrayImg( dna );
		Utils.setValues( interestPointsImage, 0.0 );

		final double[] origin = { 0, 0, 0 };
		drawTransformedPoint( alignmentTransform, interestPointsImage, origin, 200 );

		for ( int p = 0; p < 2; p++ )
		{
			drawTransformedPoint( alignmentTransform, interestPointsImage, spindlePoles.get( p ), 100 );
		}
		return interestPointsImage;
	}

	public double measureWidthByRadialProfile( RandomAccessibleInterval< T > image, final String name )
	{
		final double[] center = { 0, 0 };
		final double spacing = settings.workingVoxelSize;
		double maxDistanceInMicrometer = 15;
		final CoordinatesAndValues radialProfile = Algorithms.computeRadialProfile( image, center, spacing, maxDistanceInMicrometer );
		final CoordinatesAndValues radialProfileDerivative = CurveAnalysis.computeDerivatives( radialProfile, (int) Math.ceil( settings.derivativeDelta / settings.workingVoxelSize ) );

		if ( settings.showIntermediateResults ) Plots.plot( radialProfile, "center distance [um]", name + " intensity" );
		if ( settings.showIntermediateResults ) Plots.plot( radialProfileDerivative, "center distance [um]", "d/dx "+ name + " intensity" );

		return CurveAnalysis.minLocCoordinate( radialProfileDerivative );
	}

	public ImagePlus getOutputImage( )
	{
		final ArrayList< RandomAccessibleInterval< T > > list = new ArrayList<>();
		list.add( transformedDapiView );
		list.add( transformedTubulinView );
		list.add( transformedInterestPointView );

		RandomAccessibleInterval< T > image = Views.stack( list );
		image = Views.permute( image, 2, 3 );

		final ImagePlus output = ImageJFunctions.wrap( image, "output" );

		final Calibration calibration = new Calibration();
		calibration.pixelHeight = settings.workingVoxelSize;
		calibration.pixelWidth = settings.workingVoxelSize;
		calibration.pixelDepth = settings.workingVoxelSize;
		output.setCalibration( calibration );

		return output;
	}


	public void drawTransformedPoint( AffineTransform3D alignmentTransform, RandomAccessibleInterval< T > interestPointsImage, double[] point, int value )
	{
		final double[] transformedPoint = transformToCoverslipBasedCoordinateSystem( alignmentTransform, point );
		drawPoint( interestPointsImage, transformedPoint, settings.interestPointsRadius, settings.workingVoxelSize, value );
	}

	public HashMap< Integer, Map< String, Object > > getObjectMeasurements()
	{
		return objectMeasurements;
	}

	public RandomAccessibleInterval< BitType > createProcessedMetaPhasePlate( RandomAccessibleInterval< BitType > mask, Img< BitType > metaphasePlate )
	{
		Logger.log( "Perform morphological filtering on dapi mask..." );

		RandomAccessibleInterval< BitType > filtered;

		if ( settings.erosionOfDapiMaskInCalibratedUnits > 0 )
		{
			filtered = Algorithms.erode( metaphasePlate, ( int ) Math.ceil( settings.erosionOfDapiMaskInCalibratedUnits / settings.workingVoxelSize ) );
		}
		else
		{
			filtered = mask;
		}
		return filtered;
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


	private double[] transformToCoverslipBasedCoordinateSystem( AffineTransform3D alignmentTransform, double[] point )
	{
		final double[] transformedPoint = new double[ point.length ];
		for ( int i = 0; i < point.length; i++ )
		{
			transformedPoint[ i ] = point[ i ];
		}

		Utils.divide( transformedPoint, settings.workingVoxelSize );
		alignmentTransform.inverse().apply( transformedPoint, transformedPoint );

		return transformedPoint;
	}

	public void drawPoint( RandomAccessibleInterval< T > rai, double[] position, double calibratedRadius, double calibration, int value )
	{
		Shape shape = new HyperSphereShape( (int) Math.ceil( calibratedRadius / calibration ) );
		final RandomAccessible< Neighborhood< T > > nra = shape.neighborhoodsRandomAccessible( rai );
		final RandomAccess< Neighborhood< T > > neighborhoodRandomAccess = nra.randomAccess();

		neighborhoodRandomAccess.setPosition( Utils.asLongs( position )  );
		final Neighborhood< T > neighborhood = neighborhoodRandomAccess.get();

		final Cursor< T > cursor = neighborhood.cursor();
		while( cursor.hasNext() )
		{
			try
			{
				cursor.next().setReal( value );
			}
			catch ( ArrayIndexOutOfBoundsException e )
			{
				Logger.log( "[ERROR] Draw points out of bounds..." );
				break;
			}
		}
	}

	public static double[] getLeftAndRightMaxLocs( CoordinatesAndValues tubulinProfile, ArrayList< Double > tubulinProfileAbsoluteDerivative )
	{
		double[] rangeMinMax = new double[ 2 ];
		double[] maxLocs = new double[ 2 ];

		// left
		rangeMinMax[ 0 ] = - Double.MAX_VALUE;
		rangeMinMax[ 1 ] = 0;
		maxLocs[ 1 ] = Utils.computeMaxLoc( tubulinProfile.coordinates, tubulinProfileAbsoluteDerivative, rangeMinMax );

		// right
		rangeMinMax[ 0 ] = 0;
		rangeMinMax[ 1 ] = Double.MAX_VALUE;
		maxLocs[ 0 ] = Utils.computeMaxLoc( tubulinProfile.coordinates, tubulinProfileAbsoluteDerivative, rangeMinMax );

		return maxLocs;
	}

	public static double[] getLeftMaxAndRightMinLoc( ArrayList< Double > coordinates, ArrayList< Double > derivative )
	{
		double[] rangeMinMax = new double[ 2 ];
		double[] maxLocs = new double[ 2 ];

		// left
		rangeMinMax[ 0 ] = - Double.MAX_VALUE;
		rangeMinMax[ 1 ] = 0;
		maxLocs[ 1 ] = Utils.computeMaxLoc( coordinates, derivative, rangeMinMax );

		// right
		rangeMinMax[ 0 ] = 0;
		rangeMinMax[ 1 ] = Double.MAX_VALUE;
		maxLocs[ 0 ] = Utils.computeMinLoc( coordinates, derivative, rangeMinMax );


		return maxLocs;
	}

	public < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< BitType > createMask( RandomAccessibleInterval< T > downscaled, double threshold )
	{
		RandomAccessibleInterval< BitType > mask = Converters.convert( downscaled, ( i, o ) -> o.set( i.getRealDouble() > threshold ? true : false ), new BitType() );

		mask = opService.morphology().fillHoles( mask );

		return mask;
	}



}
