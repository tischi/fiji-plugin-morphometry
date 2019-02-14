package de.embl.cba.morphometry.spindle;

import de.embl.cba.morphometry.*;
import de.embl.cba.morphometry.geometry.CoordinatesAndValues;
import de.embl.cba.morphometry.geometry.CurveAnalysis;
import de.embl.cba.morphometry.geometry.ellipsoids.EllipsoidVectors;
import de.embl.cba.morphometry.geometry.ellipsoids.Ellipsoids3DImageSuite;
import de.embl.cba.morphometry.measurements.Measurements;
import de.embl.cba.morphometry.regions.Regions;
import de.embl.cba.transforms.utils.Transforms;
import ij.CompositeImage;
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

	public static final String SPINDLE_AXIAL_EXTEND = "Spindle_Length";
	public static final String SPINDLE_LATERAL_EXTEND = "Spindle_Width";
	public static final String DNA_AXIAL_EXTEND = "DNA_Width";
	public static final String DNA_LATERAL_EXTEND = "DNA_Length";
	public static final String DNA_VOLUME = "DNA_Volume";
	public static final String SPINDLE_VOLUME = "Spindle_Volume";

	public static final String DNA_RELATIVE_CENTRAL_INTENSITY = "DNA_Normalised_Central_Intensity";
	public static final String DNA_SPINDLE_CENTER_DISTANCE = "Dna_Center_To_Spindle_Center_Distance";
	public static final String SPINDLE_AXIS_TO_COVERSLIP_PLANE_ANGLE_DEGREES = "Spindle_Axis_To_Coverslip_Plane_Angle_Degrees";
	public static final String LENGTH_UNIT = "um";
	public static final String VOLUME_UNIT = "um3";
	public static final int ALIGNED_DNA_AXIS = 2;
	public static final String SEP = "_";

	final SpindleMorphometrySettings settings;
	final OpService opService;

	private HashMap< Integer, Map< String, Object > > objectMeasurements;
	private RandomAccessibleInterval< BitType > mask;
	private RandomAccessibleInterval transformedDNAView;
	private RandomAccessibleInterval transformedSpindleView;
	private RandomAccessibleInterval transformedInterestPointView;
	private RandomAccessibleInterval transformedDNAVolumeView;
	private RandomAccessibleInterval transformedSpindleVolumeView;

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

		final RandomAccessibleInterval< T > dna = createRescaledArrayImg( settings.dnaImage, getScalingFactors( settings.inputCalibration, settings.workingVoxelSize ) );
		final RandomAccessibleInterval< T > tubulin = createRescaledArrayImg( settings.tubulinImage, getScalingFactors( settings.inputCalibration, settings.workingVoxelSize ) );

		if ( settings.showIntermediateResults ) show( dna, "DNA isotropic resolution", null, workingCalibration, false );
		if ( settings.showIntermediateResults ) show( tubulin, "tubulin isotropic resolution", null, workingCalibration, false );


		/**
		 *  Compute DNA offset and threshold
		 */

		final RandomAccessibleInterval< T > dnaDownscaledToSpindleWidth = createRescaledArrayImg(
				settings.dnaImage,
				getScalingFactors( settings.inputCalibration, 3.0 ) );
		final double maximumValue = Algorithms.getMaximumValue( dnaDownscaledToSpindleWidth );
		double dnaThreshold = maximumValue / 2.0;
		Logger.log( "DNA threshold: " + dnaThreshold );


		/**
		 * Extract metaphase plate object
		 */

		final RandomAccessibleInterval< BitType > metaphasePlateMask = createCentralObjectMask( dna, dnaThreshold );

		if ( settings.showIntermediateResults ) show( metaphasePlateMask, "meta-phase plate", null, workingCalibration, false );


		/**
		 * Morphological filtering
		 *
		 * - it appears that the principal axes determination is not working robustly
		 * if the meta-phase plate is too thick
		 */

		final RandomAccessibleInterval< BitType > processedMetaPhasePlate = createProcessedMetaPhasePlate( metaphasePlateMask );

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

		if ( settings.showIntermediateResults ) show( alignedDNA, "aligned DNA", Transforms.origin(), workingCalibration, false );
		if ( settings.showIntermediateResults ) show( alignedProcessedMetaphasePlate, "aligned processed meta-phase plate", Transforms.origin(), workingCalibration, false );
//		if ( settings.showIntermediateResults ) Viewer3D.show3D( alignedProcessedMetaphasePlate );

		/**
		 * Compute meta-phase plate axial extend
		 */

		Logger.log( "Measuring meta-phase plate axial extend..." );

		final CoordinatesAndValues dnaProfileAlongDnaAxis = Utils.computeAverageIntensitiesAlongAxis( alignedDNA, settings.maxShortAxisDist, ALIGNED_DNA_AXIS, settings.workingVoxelSize );
		if ( settings.showIntermediateResults ) Plots.plot( dnaProfileAlongDnaAxis, "distance to center", "DNA intensity along shortest axis" );

		final CoordinatesAndValues dnaProfileAlongDnaAxisDerivative = CurveAnalysis.derivative( dnaProfileAlongDnaAxis, ( int ) Math.ceil( settings.derivativeDelta / settings.workingVoxelSize ) );
		if ( settings.showIntermediateResults ) Plots.plot( dnaProfileAlongDnaAxisDerivative, "distance to center", "d/dx DNA intensity along shortest axis" );

		final double[] dnaAxialBoundaries = CurveAnalysis.leftMaxAndRightMinLoc( dnaProfileAlongDnaAxisDerivative );

		final double dnaAxialExtend = dnaAxialBoundaries[ 1 ] - dnaAxialBoundaries[ 0 ];

		Measurements.addMeasurement(
				objectMeasurements,
				0,
				DNA_AXIAL_EXTEND + SEP + LENGTH_UNIT,
				dnaAxialExtend);

		/**
		 * Compute DNA lateral extend
		 */

		Logger.log( "Measuring meta-phase lateral extend..." );

		Projection projection = new Projection( alignedDNA, 2 );

		final RandomAccessibleInterval< T > dnaProjectionAlongDnaAxis = projection.maximum();

		final RadialWidthAndProfile dnaLateralExtendAndProfile =
				measureExtendByRadialProfile( dnaProjectionAlongDnaAxis, "dna lateral" );

		final double dnaLateralIntensityValueAtWidth = dnaLateralExtendAndProfile.profile.values.get( dnaLateralExtendAndProfile.widthIndex );

		Measurements.addMeasurement(
				objectMeasurements,
				0,
				DNA_LATERAL_EXTEND + SEP + LENGTH_UNIT,
				dnaLateralExtendAndProfile.width );


		/**
		 * Extract DNA object at the threshold determined by the lateral maximal gradient
		 */

		Logger.log( "DNA volume threshold: " + dnaLateralIntensityValueAtWidth );

		final RandomAccessibleInterval< BitType > dnaVolumeMask = createCentralObjectMask( dna, dnaLateralIntensityValueAtWidth );

		if ( settings.showIntermediateResults ) show( dnaVolumeMask, "dna volume mask", null, workingCalibration, false );


		/**
		 * Compute DNA volume
		 */

		final long dnaVolumeInPixels = Measurements.measureSizeInPixels( dnaVolumeMask );

		Measurements.addMeasurement(
				objectMeasurements,
				0,
				DNA_VOLUME + SEP + VOLUME_UNIT,
				dnaVolumeInPixels * Math.pow( settings.workingVoxelSize, 3 ) );

		/**
		 * Measure central DNA hole
		 */

		final double dnaLateralRadialProfileMaxIntensity = CurveAnalysis.maximumIndexAndValue( dnaLateralExtendAndProfile.profile ).value;
		final double dnaCenterMaxIntensity = dnaLateralExtendAndProfile.profile.values.get( 0 );
		final double dnaRelativeCentralIntensity = dnaCenterMaxIntensity / dnaLateralRadialProfileMaxIntensity;

		Measurements.addMeasurement(
				objectMeasurements,
				0,
				DNA_RELATIVE_CENTRAL_INTENSITY,
				dnaRelativeCentralIntensity );

		/**
		 * Compute spindle end-points along shortest DNA axis
		 */

		final CoordinatesAndValues tubulinProfile = Utils.computeMaximumIntensitiesAlongAxis( tubulinAlignedAlongShortestAxisOfDNA, settings.maxShortAxisDist, ALIGNED_DNA_AXIS, settings.workingVoxelSize );
		if ( settings.showIntermediateResults ) Plots.plot( tubulinProfile.coordinates, tubulinProfile.values, "center distance [um]", "tubulin max along shortest DNA axis" );

		final CoordinatesAndValues tubulinProfileDerivative =
				CurveAnalysis.derivative(
						tubulinProfile,
						(int) Math.ceil( settings.derivativeDelta / settings.workingVoxelSize ) );

		if ( settings.showIntermediateResults ) Plots.plot( tubulinProfileDerivative.coordinates, tubulinProfileDerivative.values, "distance to center", "d/dx tubulin max along shortest DNA axis" );

		double[] dnaAxisBasedSpindlePoleCoordinates =
				CurveAnalysis.leftMaxAndRightMinLoc( tubulinProfileDerivative );

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
				SPINDLE_AXIAL_EXTEND + SEP + LENGTH_UNIT,
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
					tubulinAlignedAlongShortestAxisOfDNA, "tubulin aligned along DNA axis",
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

		/**
		 * Measure spindle lateral extend in projection along spindle axis
		 */

		projection = new Projection<>( tubulinAlignedAlongSpindlePoleToPoleAxis, 2 );

		final RandomAccessibleInterval< T > spindleProjection = projection.maximum();

		final RadialWidthAndProfile tubulinLateralWidthAndProfile = measureExtendByRadialProfile( spindleProjection, "tubulin lateral" );

		Measurements.addMeasurement(
				objectMeasurements,
				0,
				SPINDLE_LATERAL_EXTEND + SEP + LENGTH_UNIT,
				tubulinLateralWidthAndProfile.width );


		/**
		 * Extract spindle object at the threshold determined by the lateral maximal gradient
		 */


		final Double tubulinVolumeThreshold = tubulinLateralWidthAndProfile.profile.values.get( tubulinLateralWidthAndProfile.widthIndex );

		Logger.log( "Spindle volume threshold: " + tubulinVolumeThreshold );

		// TODO: the spindle might "fall apart" if the threshold is too high and we might loose some bits...
		final RandomAccessibleInterval< BitType > spindleVolumeMask = createCentralObjectMask( tubulin, tubulinVolumeThreshold );

		if ( settings.showIntermediateResults ) show( spindleVolumeMask, "tubulin volume mask", null, workingCalibration, false );

		/**
		 * Compute spindle volume
		 */

		final double spindleVolume = Measurements.measureSizeInPixels( spindleVolumeMask ) * Math.pow( settings.workingVoxelSize, 3 );

		Measurements.addMeasurement(
				objectMeasurements,
				0,
				SPINDLE_VOLUME + SEP + VOLUME_UNIT,
				spindleVolume );


		/**
		 * Create interest point image in coverslip coordinate system
		 */

		RandomAccessibleInterval< T > interestPointsImage = createInterestPointImage( dna, alignmentTransform, spindlePoles );

		/**
		 * Create output image
		 */

		createOutputImage( dna, tubulin, dnaVolumeMask, spindleVolumeMask, alignmentTransform, interestPointsImage );


	}

	private RandomAccessibleInterval< BitType > createCentralObjectMask( RandomAccessibleInterval< T > image, double threshold )
	{
		RandomAccessibleInterval< BitType > mask = createMask( image, threshold );
		Regions.removeSmallRegionsInMask( mask, settings.minimalMetaphasePlateVolumeInCalibratedUnits, settings.workingVoxelSize );
		final ImgLabeling< Integer, IntType > labelImg = Utils.asImgLabeling( mask, ConnectedComponents.StructuringElement.FOUR_CONNECTED );
		final long radius = (long) ( settings.centralObjectRegionToleranceInCalibratedUnits / settings.workingVoxelSize );
		final LabelRegion< Integer > centralRegion = Regions.getCentralRegion( labelImg, radius );
		return Algorithms.createMaskFromLabelRegion( centralRegion, Intervals.dimensionsAsLongArray( labelImg ) );
	}

	private void createOutputImage(
			RandomAccessibleInterval< T > dna,
			RandomAccessibleInterval< T > spindle,
			RandomAccessibleInterval< BitType > dnaVolumeMask,
			RandomAccessibleInterval< BitType > spindleVolumeMask,
			AffineTransform3D alignmentTransform,
			RandomAccessibleInterval< T > interestPointsImage )
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

		transformedDNAView = Transforms.createTransformedView( dna, rotationTransform );
		transformedSpindleView = Transforms.createTransformedView( spindle, rotationTransform );
		transformedDNAVolumeView = Transforms.createTransformedView( dnaVolumeMask, rotationTransform );
		transformedSpindleVolumeView = Transforms.createTransformedView( spindleVolumeMask, rotationTransform );
		transformedInterestPointView = Transforms.createTransformedView( interestPointsImage, rotationTransform, new NearestNeighborInterpolatorFactory() );

	}

	private RandomAccessibleInterval< T > createInterestPointImage(
			RandomAccessibleInterval< T > dna,
			AffineTransform3D alignmentTransform,
			ArrayList< double[] > spindlePoles )
	{
		RandomAccessibleInterval< T > interestPointsImage = Utils.copyAsArrayImg( dna );
		Utils.setValues( interestPointsImage, 0.0 );

		final double[] origin = { 0, 0, 0 };
		drawTransformedPoint( alignmentTransform, interestPointsImage, origin, 255 );

		for ( int p = 0; p < 2; p++ )
		{
			drawTransformedPoint( alignmentTransform, interestPointsImage, spindlePoles.get( p ), 255 );
		}
		return interestPointsImage;
	}

	class RadialWidthAndProfile
	{
		Double width;
		CoordinatesAndValues profile;
		int widthIndex;
	}

	private RadialWidthAndProfile measureExtendByRadialProfile(
			final RandomAccessibleInterval< T > image,
			final String name
	)
	{
		// config
		final double[] center = { 0, 0 };
		final double spacing = settings.workingVoxelSize;
		final double maxDistanceInMicrometer = 15;

		// measure
		final RadialWidthAndProfile widthAndProfile = new RadialWidthAndProfile();
		widthAndProfile.profile = Algorithms.computeRadialProfile( image, center, spacing, maxDistanceInMicrometer );

		final CoordinatesAndValues radialProfileDerivative = CurveAnalysis.derivative( widthAndProfile.profile, (int) Math.ceil( settings.derivativeDelta / settings.workingVoxelSize ) );
		widthAndProfile.width = CurveAnalysis.minLoc( radialProfileDerivative );

		widthAndProfile.widthIndex = (int) (widthAndProfile.width / settings.workingVoxelSize);

		// plot
 		if ( settings.showIntermediateResults ) Plots.plot( widthAndProfile.profile, "center distance [um]", name + " intensity" );
		if ( settings.showIntermediateResults ) Plots.plot( radialProfileDerivative, "center distance [um]", "d/dx "+ name + " intensity" );

		return widthAndProfile;
	}

	public CompositeImage getOutputImage( )
	{
		final ArrayList< RandomAccessibleInterval< T > > list = new ArrayList<>();
		list.add( transformedDNAView );
		list.add( transformedDNAVolumeView );
		list.add( transformedSpindleView );
		list.add( transformedSpindleVolumeView );
		list.add( transformedInterestPointView );

		RandomAccessibleInterval< T > image = Views.stack( list );
		image = Views.permute( image, 2, 3 );

		final ImagePlus imp = ImageJFunctions.wrap( image, "output" );

		final Calibration calibration = new Calibration();
		calibration.pixelHeight = settings.workingVoxelSize;
		calibration.pixelWidth = settings.workingVoxelSize;
		calibration.pixelDepth = settings.workingVoxelSize;
		imp.setCalibration( calibration );

		final CompositeImage compositeImage = new CompositeImage( imp );

		compositeImage.setC(1);
		IJ.run(compositeImage, "Blue", "");

		compositeImage.setC(2);
		IJ.run(compositeImage, "Blue", "");
		compositeImage.setDisplayRange( 0, 2 );

		compositeImage.setC(3);
		IJ.run(compositeImage, "Green", "");

		compositeImage.setC(4);
		IJ.run(compositeImage, "Green", "");
		compositeImage.setDisplayRange( 0, 2);

		compositeImage.setC(5);
		IJ.run(compositeImage, "Cyan", "");

		compositeImage.setDisplayMode( CompositeImage.COMPOSITE );

		return compositeImage;
	}


	private void drawTransformedPoint( AffineTransform3D alignmentTransform, RandomAccessibleInterval< T > interestPointsImage, double[] point, int value )
	{
		final double[] transformedPoint = transformToCoverslipBasedCoordinateSystem( alignmentTransform, point );
		drawPoint( interestPointsImage, transformedPoint, settings.interestPointsRadius, settings.workingVoxelSize, value );
	}

	public HashMap< Integer, Map< String, Object > > getObjectMeasurements()
	{
		return objectMeasurements;
	}

	private RandomAccessibleInterval< BitType > createProcessedMetaPhasePlate( RandomAccessibleInterval< BitType > metaphasePlate )
	{
		Logger.log( "Perform morphological filtering on dna mask..." );

		RandomAccessibleInterval< BitType > filtered;

		if ( settings.erosionOfDnaMaskInCalibratedUnits > 0 )
		{
			filtered = Algorithms.erode( metaphasePlate, ( int ) Math.ceil( settings.erosionOfDnaMaskInCalibratedUnits / settings.workingVoxelSize ) );
		}
		else
		{
			filtered = metaphasePlate;
		}

		return filtered;
	}

	private void saveMaximumProjections( RandomAccessibleInterval transformedDNAView, String name )
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

			IJ.saveAsTiff( ImageJFunctions.wrap( new Projection( transformedDNAView, d ).maximum(), title ), path );
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

	private void drawPoint( RandomAccessibleInterval< T > rai, double[] position, double calibratedRadius, double calibration, int value )
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

	private < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< BitType > createMask( RandomAccessibleInterval< T > downscaled, double threshold )
	{
		RandomAccessibleInterval< BitType > mask = Converters.convert( downscaled, ( i, o ) -> o.set( i.getRealDouble() > threshold ? true : false ), new BitType() );

		mask = opService.morphology().fillHoles( mask );

		return mask;
	}



}
