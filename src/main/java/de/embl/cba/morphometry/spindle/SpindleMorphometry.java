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
import net.imglib2.RandomAccess;
import net.imglib2.algorithm.labeling.ConnectedComponents;
import net.imglib2.algorithm.neighborhood.HyperSphereShape;
import net.imglib2.algorithm.neighborhood.Neighborhood;
import net.imglib2.algorithm.neighborhood.Shape;
import net.imglib2.converter.Converters;
import net.imglib2.img.Img;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.interpolation.randomaccess.NLinearInterpolatorFactory;
import net.imglib2.interpolation.randomaccess.NearestNeighborInterpolatorFactory;
import net.imglib2.realtransform.AffineTransform2D;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.realtransform.RealViews;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.roi.labeling.LabelRegion;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.util.Intervals;
import net.imglib2.util.LinAlgHelpers;
import net.imglib2.view.IntervalView;
import net.imglib2.view.Views;

import java.io.File;
import java.util.*;

import static de.embl.cba.morphometry.Angles.angleOfSpindleAxisToXAxisInRadians;
import static de.embl.cba.morphometry.viewing.BdvViewer.show;
import static de.embl.cba.transforms.utils.Scalings.createRescaledArrayImg;
import static de.embl.cba.transforms.utils.Transforms.getScalingFactors;


public class SpindleMorphometry  < T extends RealType< T > & NativeType< T > >
{

	public static final String SEP = "_";

	final SpindleMorphometrySettings< T > settings;
	final OpService opService;

	private HashMap< Integer, Map< String, Object > > objectMeasurements;
	private RandomAccessibleInterval< BitType > mask;
	private RandomAccessibleInterval transformedDNAView;
	private RandomAccessibleInterval transformedSpindleView;
	private RandomAccessibleInterval transformedInterestPointView;
	private RandomAccessibleInterval transformedDNAVolumeView;
	private RandomAccessibleInterval transformedSpindleVolumeView;
	private AffineTransform3D dnaAxesAlignmentTransform;
	private RandomAccessibleInterval dnaAlignedDna;
	private RandomAccessibleInterval dnaAlignedDnaMask;
	private RandomAccessibleInterval dnaAlignedTubulin;
	private RandomAccessibleInterval< T > dna;
	private RandomAccessibleInterval< T > tubulin;
	private ArrayList< double[] > spindlePoles;
	private RandomAccessibleInterval< BitType > dnaVolumeMask;
	private RandomAccessibleInterval< BitType > spindleVolumeMask;
	private double[] workingCalibration;
	private RandomAccessibleInterval< BitType > segmentedDna;
	private EllipsoidVectors dnaEllipsoidVectors;
	private ProfileAndRadius dnaLateralProfileAndRadius;
	private double[] spindlePoleToPoleVector;
	private double[] spindleCenter;
	private ProfileAndRadius spindleLateralRadiusAndProfile;
	private RandomAccessibleInterval< T > poleToPoleAlignedSpindleRai;
	private double spindleThreshold;
	private double dnaLateralExtend;
	private double spindleAxialExtend;
	private SpindleMeasurements spindleMeasurements;

	public SpindleMorphometry( SpindleMorphometrySettings settings, OpService opService )
	{
		this.settings = settings;
		this.opService = opService;
	}

	public String run()
	{
		objectMeasurements = new HashMap<>();

		spindleMeasurements = new SpindleMeasurements( objectMeasurements );

		try
		{
			spindleMeasurements.log = measure();
		}
		catch ( Exception e )
		{
			e.printStackTrace();
			spindleMeasurements.log = "Exception during computation.";
		}

		spindleMeasurements.setObjectMeasurements();

		spindleMeasurements.log = SpindleMeasurements.ANALYSIS_FINISHED;

		return spindleMeasurements.log;
	}

	private String measure()
	{

		/**
		 * TODO:
		 * - maybe smooth the dna and tubulin image to reduce noise?
		 * - for low dynamic ranges, the smoothing should be done in a doubletype image
		 */

		createIsotropicallyResampledImages();

		double dnaThreshold = determineDnaThreshold();

		if ( dnaThreshold < settings.minimalDynamicRange )
			return SpindleMeasurements.ANALYSIS_INTERRUPTED_LOW_DYNAMIC_DNA;

		segmentedDna = segmentDna( dna, dnaThreshold );

		dnaEllipsoidVectors = determineDnaAxes( segmentedDna );

		computeDnaAlignmentTransformAndAlignImages(
				dna, tubulin, segmentedDna, dnaEllipsoidVectors );

		measureDnaAxialExtend( dnaAlignedDna );

		dnaLateralProfileAndRadius = measureDnaLateralExtend( dnaAlignedDna );

		dnaVolumeMask = measureDnaVolume( dna, dnaLateralProfileAndRadius );

		measureDnaHole( dnaLateralProfileAndRadius );

		spindlePoles = measureSpindlePoleLocations( dnaAlignedTubulin );

		spindlePoleToPoleVector = getSpindleAxisVector( spindlePoles );

		spindleCenter = getSpindleCenter( spindlePoles, spindlePoleToPoleVector );

		poleToPoleAlignedSpindleRai =
				createSpindlePolesAlignedRai( spindlePoles, spindleCenter );

		spindleThreshold =
				measureSpindleThreshold( poleToPoleAlignedSpindleRai );

		if ( spindleThreshold < settings.minimalDynamicRange )
			return SpindleMeasurements.ANALYSIS_INTERRUPTED_LOW_DYNAMIC_TUBULIN;

		spindleVolumeMask =
				measureSpindleVolume( poleToPoleAlignedSpindleRai, spindleThreshold );

		measureSpindleLateralExtends( spindleVolumeMask );

		measureDnaCenterToSpindleCenterDistance(
				workingCalibration, spindlePoles, spindleCenter );

		measureSpindleAxisToCoverslipPlaneAngle( spindlePoleToPoleVector );

		return "No comment";
	}

	private void measureSpindleLateralExtends(
			RandomAccessibleInterval< BitType > alignedSpindleMask )
	{
		Logger.log( "Creating projection of spindle mask along spindle axis..." );

		// TODO: only include pixels within the spindle range

		final RandomAccessibleInterval< BitType > projectedMask =
				new Projection<>(
					alignedSpindleMask,
					2 ).maximum();

		Regions.onlyKeepLargestRegion( projectedMask, ConnectedComponents.StructuringElement.EIGHT_CONNECTED  );

		if ( settings.showIntermediateResults )
			show( projectedMask,
					"Spindle mask projection along pole to pole axis",
					null,
					workingCalibration,
					false);

		final ArrayList< Long > widths = measureRadialWidthsInPixels( projectedMask );
		Collections.sort( widths );

		final double[] calibratedMinMaxWidths = {
				widths.get( 0 ) * workingCalibration[ 0 ],
				widths.get( widths.size() - 1 ) * workingCalibration[ 1 ]
		};

		spindleMeasurements.spindleWidthMin = calibratedMinMaxWidths[ 0 ];
		spindleMeasurements.spindleWidthMax = calibratedMinMaxWidths[ 1 ];

	}

	public RandomAccessibleInterval< BitType > measureSpindleVolume(
			RandomAccessibleInterval< T > poleToPoleAlignedSpindle,
			double spindleThreshold )
	{
		final RandomAccessibleInterval< BitType > spindleMask =
				createSpindleMask( poleToPoleAlignedSpindle, spindleThreshold );

		if ( settings.showIntermediateResults )
			show( spindleMask,
					"spindle volume mask", null,
					workingCalibration, false );

		/**
		 * Compute spindle volume
		 */

		spindleMeasurements.spindleVolume =
				Measurements.measureSizeInPixels( spindleMask )
						* Math.pow( settings.workingVoxelSize, 3 );

		return spindleMask;
	}


	/**
	 *
	 * @param poleToPoleAlignedSpindleRai
	 * @return
	 */
	public double measureSpindleThreshold(
			RandomAccessibleInterval< T > poleToPoleAlignedSpindleRai )
	{
		final RandomAccessibleInterval< T > spindleProjection =
				createMaximumProjectionAlongSpindleAxis( poleToPoleAlignedSpindleRai );

		final ArrayList< Double > thresholdCandidates =
				measureRadialThresholds( spindleProjection );

		spindleMeasurements.spindleThreshold = Utils.median( thresholdCandidates );

		return spindleMeasurements.spindleThreshold ;
	}

	private RandomAccessibleInterval< T > createMaximumProjectionAlongSpindleAxis(
			RandomAccessibleInterval< T > poleToPoleAlignedSpindleRai )
	{
		Logger.log( "Computing maximum projection of spindle along spindle axis..." );

		final FinalInterval spindleInterval = createSpindleInterval();

		Projection projection = new Projection<>(
				Views.interval(  poleToPoleAlignedSpindleRai, spindleInterval ),
				2 );

		final RandomAccessibleInterval< T > maximum = projection.maximum();

		if ( settings.showIntermediateResults )
			show( maximum, "Spindle maximum projection along pole to pole axis",
					null, workingCalibration, false);

		return maximum;
	}

	private FinalInterval createSpindleInterval()
	{
		long[] min = new long[3];
		long[] max = new long[3];

		max[ 0 ] = (long) Math.ceil( dnaLateralExtend / 2.0 / settings.workingVoxelSize );
		max[ 1 ] = (long) Math.ceil( dnaLateralExtend / 2.0 / settings.workingVoxelSize );
		max[ 2 ] = (long) Math.ceil( 1.2 * spindleAxialExtend / 2.0 / settings.workingVoxelSize );

		for ( int d = 0; d < 3; d++ )
			min[ d ] = - max[ d ];

		return new FinalInterval( min, max );
	}


	/**
	 *
	 * Align spindle along pole to pole axis
	 *
	 * @param spindlePoles
	 * @param spindleCenter
	 * @return
	 */
	private RandomAccessibleInterval< T > createSpindlePolesAlignedRai(
			ArrayList< double[] > spindlePoles,
			double[] spindleCenter )
	{
		final double[] poleToPoleAxis = new double[ 3 ];
		LinAlgHelpers.subtract( spindlePoles.get( 0 ), spindlePoles.get( 1 ), poleToPoleAxis );
		LinAlgHelpers.normalize( poleToPoleAxis );

		AffineTransform3D spindleTransform = new AffineTransform3D();
		spindleTransform.translate( spindleCenter );
		spindleTransform = spindleTransform.inverse();
		AffineTransform3D poleToPoleAxisRotation =
				Transforms.getRotationTransform3D( new double[]{ 0, 0, 1 }, poleToPoleAxis );
		spindleTransform = spindleTransform.preConcatenate( poleToPoleAxisRotation );

		final RandomAccessibleInterval< T > spindleAlignedAlongSpindlePoleToPoleAxis
				= Transforms.createTransformedView( dnaAlignedTubulin, spindleTransform );

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

			show( spindleAlignedAlongSpindlePoleToPoleAxis,
					"spindle aligned pole to pole",
					realPoints,
					workingCalibration,
					false );
		}

		return spindleAlignedAlongSpindlePoleToPoleAxis;
	}

	@Deprecated
	public void measureRadialWidthsOld( RandomAccessibleInterval< T > rai )
	{
		RealRandomAccessible< T > rra =
				Views.interpolate(
						Views.extendZero( rai ), new NLinearInterpolatorFactory<>() );

		double dAngle = Math.PI / 3;

		for ( double angle = 0; angle < Math.PI; angle += dAngle )
		{
			final AffineTransform2D transform2D = new AffineTransform2D();
			transform2D.rotate( angle );

			IntervalView< T > rotated =
					Views.interval(
							Views.raster(
									RealViews.transform( rra, transform2D ) ), rai );

			final CoordinatesAndValues coordinatesAndValues =
					Utils.computeAverageIntensitiesAlongAxis(
							rotated,
							1.0,
							0,
							settings.workingVoxelSize );

			final CoordinatesAndValues derivative =
					CurveAnalysis.derivative(
							coordinatesAndValues,
							( int ) Math.ceil( settings.spindleDerivativeDelta
									/ settings.workingVoxelSize ) );


//			if ( settings.showIntermediateResults )
			{
//				Plots.plot(
//						coordinatesAndValues, "Center distance [um]",
//						"Spindle lateral intensity, angle = " + angle );
				Plots.plot( derivative, "Center distance [um]",
						"d/dx Spindle lateral intensity, angle = " + angle  );
			}

		}

	}

	public ArrayList< Long > measureRadialWidthsInPixels(
			RandomAccessibleInterval< BitType > mask )
	{
		RealRandomAccessible< BitType > rra =
				Views.interpolate(
						Views.extendZero( mask ), new NearestNeighborInterpolatorFactory<>() );

		double dAngle = Math.PI / 18;

		final ArrayList< Long > numPixels = new ArrayList<>();

		for ( double angle = 0; angle < Math.PI; angle += dAngle )
		{
			final AffineTransform2D transform2D = new AffineTransform2D();
			transform2D.rotate( angle );

			IntervalView< BitType > rotated =
					Views.interval(
							Views.raster(
									RealViews.transform( rra, transform2D ) ), mask );

			numPixels.add( Utils.countNonZeroPixelsAlongAxis( rotated, 0 ) );


////			if ( settings.showIntermediateResults )
//			{
////				Plots.plot(
////						coordinatesAndValues, "Center distance [um]",
////						"Spindle lateral intensity, angle = " + angle );
//				Plots.plot( derivative, "Center distance [um]",
//						"d/dx Spindle lateral intensity, angle = " + angle  );
//			}

		}

		return numPixels;

	}

	public ArrayList< Double > measureRadialThresholds( RandomAccessibleInterval< T > rai )
	{

		final ArrayList< Double > thresholdCandidates = new ArrayList<>();

		RealRandomAccessible< T > rra =
				Views.interpolate(
						Views.extendZero( rai ), new NLinearInterpolatorFactory<>() );

		double dAngle = Math.PI / 5;

		for ( double angle = 0; angle < Math.PI; angle += dAngle )
		{
			final AffineTransform2D transform2D = new AffineTransform2D();
			transform2D.rotate( angle );

			IntervalView< T > rotated =
					Views.interval(
							Views.raster(
									RealViews.transform( rra, transform2D ) ), rai );

			final CoordinatesAndValues intensities =
					Utils.computeAverageIntensitiesAlongAxis(
							rotated,
							1.0,
							0,
							settings.workingVoxelSize );

			final CoordinatesAndValues derivatives =
					CurveAnalysis.derivative(
							intensities,
							( int ) Math.ceil( settings.spindleDerivativeDelta
									/ settings.workingVoxelSize ) );

			final ArrayList< CoordinateAndValue > extrema
					= CurveAnalysis.leftMaxAndRightMinLoc( derivatives );


			final double leftLoc = extrema.get( 0 ).coordinate - settings.spindleDerivativeDelta / 2.0;

			thresholdCandidates.add( CurveAnalysis.getValueAtCoordinate(
					intensities,
					leftLoc ) );

			final double rightLoc = extrema.get( 1 ).coordinate + settings.spindleDerivativeDelta / 2.0;

			thresholdCandidates.add( CurveAnalysis.getValueAtCoordinate(
					intensities,
					rightLoc ) );

			if ( settings.showIntermediateResults )
			{
				Plots.plot( intensities, "Center distance [um]",
						"Spindle lateral intensity, angle = " + angle  );
				Plots.plot( derivatives, "Center distance [um]",
						"d/dx Spindle lateral intensity, angle = " + angle  );
			}

		}

		return thresholdCandidates;
	}


	public void measureSpindleAxisToCoverslipPlaneAngle( double[] poleToPoleVector )
	{
		final double[] poleToPoleVectorInCSCS =
				transformToPixelUnitsCoverslipCoordinateSystem( dnaAxesAlignmentTransform, poleToPoleVector );

		spindleMeasurements.angleSpindleAxisToCoverslipPlaneInDegrees =
				90.0 - Math.abs( 180.0 / Math.PI *
						Transforms.getAngle( new double[]{ 0, 0, 1 }, poleToPoleVectorInCSCS ) );
	}

	public void measureDnaCenterToSpindleCenterDistance(
			double[] workingCalibration,
			ArrayList< double[] > spindlePoles,
			double[] spindleCenter )
	{
		spindleMeasurements.dnaCenterToSpindleCenterDistance =
				LinAlgHelpers.distance( new double[]{ 0, 0, 0}, spindleCenter) * settings.workingVoxelSize;

		if ( settings.showIntermediateResults )
		{
			final ArrayList< RealPoint > realPoints = new ArrayList<>();
			realPoints.add( new RealPoint( spindlePoles.get( 0 ) ) );
			realPoints.add( new RealPoint( spindlePoles.get( 1 ) ) );
			realPoints.add( new RealPoint( spindleCenter ) );
			show(
					dnaAlignedTubulin, "spindle aligned along DNA axis",
					realPoints,
					workingCalibration,
					false );
		}
	}

	public double[] getSpindleCenter( ArrayList< double[] > spindlePoles, double[] poleToPoleVector )
	{
		final double[] spindleCenter = new double[ 3 ];

		for ( int i = 0; i < 3; i++ )
			poleToPoleVector[ i ] *= 0.5;

		LinAlgHelpers.add( spindlePoles.get( 1 ), poleToPoleVector, spindleCenter );

		return spindleCenter;
	}

	public double[] getSpindleAxisVector( ArrayList< double[] > spindlePoles )
	{
		final double[] poleToPoleVector = new double[ 3 ];

		LinAlgHelpers.subtract( spindlePoles.get( 0 ), spindlePoles.get( 1 ), poleToPoleVector );

		return poleToPoleVector;
	}

	public void createIsotropicallyResampledImages()
	{
		Logger.log( "Create isotropic images..." );

		workingCalibration = Utils.as3dDoubleArray( settings.workingVoxelSize );

		dna = createRescaledArrayImg(
				settings.dnaImage,
				getScalingFactors( settings.inputCalibration, settings.workingVoxelSize ) );

		tubulin = createRescaledArrayImg(
				settings.tubulinImage,
				getScalingFactors( settings.inputCalibration, settings.workingVoxelSize ) );

		if ( settings.showIntermediateResults )
			show( dna, "DNA isotropic voxel size", null, workingCalibration, false );

		if ( settings.showIntermediateResults )
			show( tubulin, "spindle isotropic voxel size", null, workingCalibration, false );
	}

	public double determineDnaThreshold()
	{
		final RandomAccessibleInterval< T > dnaDownscaledToMetaphasePlateWidth =
				createRescaledArrayImg(
						dna,
						getScalingFactors(
								new double[]{
										settings.workingVoxelSize,
										settings.workingVoxelSize,
										settings.workingVoxelSize },
								settings.dnaThresholdResolution ) );

		final double maximumValue =
				Algorithms.getMaximumValue( dnaDownscaledToMetaphasePlateWidth );

//		Viewers.showRai3dWithImageJ( dnaDownscaledToMetaphasePlateWidth, "DNA Threshold" );

		Logger.log( "DNA downscaled maximum value: " + maximumValue );
		Logger.log( "DNA threshold factor: " + settings.dnaThresholdFactor );
		double dnaThreshold = maximumValue * settings.dnaThresholdFactor;
		Logger.log( "DNA threshold: " + dnaThreshold );

		return dnaThreshold;
	}

	public RandomAccessibleInterval< BitType > segmentDna(
			RandomAccessibleInterval< T > dna,
			double dnaThreshold )
	{
		final RandomAccessibleInterval< BitType > dnaMask =
				createCentralObjectsMask( dna, dnaThreshold );

		if ( settings.showIntermediateResults )
			show( dnaMask, "dna mask", null,
					workingCalibration, false );

		/**
		 * Morphological filtering
		 *
		 * it appears that the principal axes determination
		 * is not working robustly
		 * if the meta-phase plate is too thick
		 *
		 * This was mainly an issue using MLJ, now with 3D Image Suite, this issue seems solved
		 */

//		final RandomAccessibleInterval< BitType > processedDnaMask = createProcessedMetaPhasePlate( dnaMask );

//		if ( settings.showIntermediateResults )
//			show( processedDnaMask, "eroded DNA mask", null, workingCalibration, false );

		return dnaMask;
	}

	public void computeDnaAlignmentTransformAndAlignImages(
			RandomAccessibleInterval< T > dna,
			RandomAccessibleInterval< T > tubulin,
			RandomAccessibleInterval< BitType > dnaMask,
			EllipsoidVectors ellipsoidVectors )
	{
		Logger.log( "Creating aligned images..." );

		dnaAxesAlignmentTransform =
				Ellipsoids3DImageSuite.createShortestAxisAlignmentTransform( ellipsoidVectors );

		dnaAlignedTubulin = Utils.copyAsArrayImg(
				Transforms.createTransformedView( tubulin,
						dnaAxesAlignmentTransform, new NearestNeighborInterpolatorFactory() ) );

		dnaAlignedDna = Utils.copyAsArrayImg(
				Transforms.createTransformedView( dna, dnaAxesAlignmentTransform ) );

		dnaAlignedDnaMask = Utils.copyAsArrayImg(
				Transforms.createTransformedView( dnaMask, dnaAxesAlignmentTransform ) );

		if ( settings.showIntermediateResults )
			show( dnaAlignedTubulin, "tubulin aligned along DNA shortest axis", Transforms.origin(),
					workingCalibration, false );

		if ( settings.showIntermediateResults )
			show( dnaAlignedDna, "DNA aligned along DNA shortest axis", Transforms.origin(),
					workingCalibration, false );

		if ( settings.showIntermediateResults )
			show( dnaAlignedDnaMask, "DNA mask aligned along DNA shortest axis",
					Transforms.origin(), workingCalibration, false );

//		if ( settings.showIntermediateResults ) Viewer3D.show3D( dnaAlignedDnaMask );
	}

	public EllipsoidVectors determineDnaAxes( RandomAccessibleInterval< BitType > dnaMask )
	{
		Logger.log( "Determining DNA axes..." );

		return Ellipsoids3DImageSuite.fitEllipsoid( Utils.getAsImagePlusMovie( dnaMask, "" ) );
	}

	public void measureDnaAxialExtend( RandomAccessibleInterval alignedDNA )
	{
		Logger.log( "Measuring DNA axial extend..." );

		final CoordinatesAndValues dnaProfileAlongDnaAxis =
				Utils.computeAverageIntensitiesAlongAxis(
						alignedDNA,
						settings.maxDnaAxisDist,
						2,
						settings.workingVoxelSize );

		final CoordinatesAndValues dnaProfileAlongDnaAxisDerivative =
				CurveAnalysis.derivative(
						dnaProfileAlongDnaAxis,
						( int ) Math.ceil( settings.derivativeDelta
								/ settings.workingVoxelSize ) );

		final ArrayList< CoordinateAndValue > dnaAxialBoundaries =
				CurveAnalysis.leftMaxAndRightMinLoc( dnaProfileAlongDnaAxisDerivative );

		spindleMeasurements.dnaAxialExtend =
				dnaAxialBoundaries.get( 1 ).coordinate -
						dnaAxialBoundaries.get( 0 ).coordinate;

		logDnaAxialThreshold( dnaProfileAlongDnaAxis, dnaAxialBoundaries );

		if ( settings.showIntermediateResults )
			Plots.plot( dnaProfileAlongDnaAxis,
					"distance to center", "DNA intensity along DNA axis" );

		if ( settings.showIntermediateResults )
			Plots.plot( dnaProfileAlongDnaAxisDerivative,
					"distance to center", "d/dx DNA intensity along DNA axis" );

	}

	public void logDnaAxialThreshold(
			CoordinatesAndValues dnaProfileAlongDnaAxis,
			ArrayList< CoordinateAndValue > dnaAxialBoundaries )
	{
		for ( int i = 0; i < 2; i++ )
		{
			final Double valueAtCoordinate =
					CurveAnalysis.getValueAtCoordinate( dnaProfileAlongDnaAxis,
							dnaAxialBoundaries.get( i ).coordinate );
			Logger.log( "DNA axial threshold " + i + ": " + valueAtCoordinate );
		}
	}

	public ProfileAndRadius measureDnaLateralExtend(
			RandomAccessibleInterval< T > alignedDNA )
	{
		Logger.log( "Measuring DNA lateral extend..." );

		Projection projection = new Projection( alignedDNA, 2 );

		final RandomAccessibleInterval< T > dnaProjectionAlongDnaAxis = projection.maximum();

		final ProfileAndRadius dnaLateralProfileAndRadius =
				measureRadialProfileAndRadius(
						dnaProjectionAlongDnaAxis,
						"dna lateral",
						settings.derivativeDelta );

		spindleMeasurements.dnaLateralExtend = 2.0 * dnaLateralProfileAndRadius.radius;

		return dnaLateralProfileAndRadius;
	}

	public ArrayList< double[] > measureSpindlePoleLocations(
			RandomAccessibleInterval tubulinAlignedAlongShortestAxisOfDNA )
	{
		final ArrayList< double[] > dnaAxisBasedSpindlePoles =
				getSpindlePolesAlongDnaAxis( tubulinAlignedAlongShortestAxisOfDNA );

		final ArrayList< double[] > spindlePoles = getRefinedSpindlePoles(
				tubulinAlignedAlongShortestAxisOfDNA,
				dnaAxisBasedSpindlePoles );

		addPoleRefinementDistanceMeasurements( dnaAxisBasedSpindlePoles, spindlePoles );

		spindleMeasurements.spindleAxialExtend = LinAlgHelpers.distance(
				spindlePoles.get( 0 ), spindlePoles.get( 1 ) );

		return spindlePoles;
	}


	public RandomAccessibleInterval< BitType > measureDnaVolume(
			RandomAccessibleInterval< T > dna,
			ProfileAndRadius dnaLateralExtendAndProfile )
	{
		spindleMeasurements.dnaVolumeThreshold =
				dnaLateralExtendAndProfile.profile.values.get(
						dnaLateralExtendAndProfile.radiusIndex );

		final RandomAccessibleInterval< BitType > dnaVolumeMask =
				createCentralObjectsMask( dna, spindleMeasurements.dnaVolumeThreshold );

		final long dnaVolumeInPixels =
				Measurements.measureSizeInPixels( dnaVolumeMask );

		spindleMeasurements.dnaVolumeCalibrated =
				dnaVolumeInPixels * Math.pow( settings.workingVoxelSize, 3 );


		if ( settings.showIntermediateResults )
			show( dnaVolumeMask, "dna volume mask",
					null, workingCalibration, false );

		return dnaVolumeMask;
	}

	public ArrayList< double[] > getSpindlePolesAlongDnaAxis(
			RandomAccessibleInterval tubulinAlignedAlongShortestAxisOfDNA )
	{
		final CoordinatesAndValues tubulinProfile =
				Utils.computeMaximumIntensitiesAlongAxis(
						tubulinAlignedAlongShortestAxisOfDNA,
						settings.maxDnaAxisDist,
						SpindleMeasurements.ALIGNED_DNA_AXIS,
						settings.workingVoxelSize );

		final CoordinatesAndValues tubulinProfileDerivative =
				CurveAnalysis.derivative(
						tubulinProfile,
						(int) Math.ceil( settings.derivativeDelta
								/ settings.workingVoxelSize ) );

		ArrayList< CoordinateAndValue > tubulinExtrema =
				CurveAnalysis.leftMaxAndRightMinLoc( tubulinProfileDerivative );

		final ArrayList< double[] > dnaAxisBasedSpindlePoles = new ArrayList<>();
		dnaAxisBasedSpindlePoles.add( new double[]{ 0, 0, tubulinExtrema.get( 0 ).coordinate });
		dnaAxisBasedSpindlePoles.add( new double[]{ 0, 0, tubulinExtrema.get( 1 ).coordinate });

		for ( int p = 0; p < 2; p++ )
		{
			final Double value =
					CurveAnalysis.getValueAtCoordinate(
							tubulinProfile,
							tubulinExtrema.get( p ).coordinate );
			Logger.log( "Spindle axial threshold " + p + ": " + value );
		}

		if ( settings.showIntermediateResults )
			Plots.plot(
					tubulinProfile.coordinates,
					tubulinProfile.values,
					"center distance [um]",
					"spindle max along shortest DNA axis" );

		if ( settings.showIntermediateResults )
			Plots.plot(
					tubulinProfileDerivative.coordinates,
					tubulinProfileDerivative.values,
					"distance to center",
					"d/dx spindle max along shortest DNA axis" );

		return dnaAxisBasedSpindlePoles;
	}

	public void measureDnaHole( ProfileAndRadius dnaLateralProfileAndRadius )
	{
		final double dnaLateralRadialProfileMaxIntensity =
				CurveAnalysis.maximumIndexAndValue( dnaLateralProfileAndRadius.profile ).value;

		final double dnaCenterMaxIntensity = dnaLateralProfileAndRadius.profile.values.get( 0 );

		spindleMeasurements.dnaRelativeCentralIntensity =
				dnaCenterMaxIntensity / dnaLateralRadialProfileMaxIntensity;

	}

	public void addPoleRefinementDistanceMeasurements(
			ArrayList< double[] > dnaAxisBasedSpindlePoles,
			ArrayList< double[] > spindlePoles )
	{
		spindleMeasurements.spindlePoleARefinementDistance =
				LinAlgHelpers.distance(
						dnaAxisBasedSpindlePoles.get( 0 ), spindlePoles.get( 0 ) );

		spindleMeasurements.spindlePoleBRefinementDistance =
				LinAlgHelpers.distance(
						dnaAxisBasedSpindlePoles.get( 1 ), spindlePoles.get( 1 ) );

	}

	private ArrayList< double[] > getRefinedSpindlePoles(
			RandomAccessibleInterval tubulinAlignedAlongShortestAxisOfDNA,
			ArrayList< double[] > dnaAxisBasedSpindlePoleCoordinates )
	{
		final ArrayList< double[] > spindlePoles = new ArrayList<>();

		for ( int pole = 0; pole < 2; pole++ )
		{
			RandomAccessibleInterval< T > tubulinIntensitySlice =
					Views.hyperSlice(
							tubulinAlignedAlongShortestAxisOfDNA,
							2,
							(long) ( dnaAxisBasedSpindlePoleCoordinates.get( pole )[ 2 ]
									/ settings.workingVoxelSize ) );

			final double[] xy = Utils.computeMaximumLocation(
					tubulinIntensitySlice,
					settings.maxSpindlePoleRefinementDistance
							/ settings.workingVoxelSize );

			spindlePoles.add( new double[] {
							xy[ 0 ] * settings.workingVoxelSize,
							xy[ 1 ] * settings.workingVoxelSize,
							dnaAxisBasedSpindlePoleCoordinates.get( pole )[ 2 ] } );
		}

		return spindlePoles;
	}

	private RandomAccessibleInterval< BitType > createCentralObjectsMask(
			RandomAccessibleInterval< T > image,
			double threshold )
	{
		RandomAccessibleInterval< BitType > mask = createMask( image, threshold );

		Regions.removeSmallRegionsInMask(
				mask,
				settings.minimalDnaFragementsVolume,
				settings.workingVoxelSize );

		final ImgLabeling< Integer, IntType > labelImg =
				Utils.asImgLabeling( mask, ConnectedComponents.StructuringElement.FOUR_CONNECTED );

		final long radius = (long) ( settings.maxCentralObjectRegionsDistance / settings.workingVoxelSize );

		final Set< LabelRegion< Integer > > centralRegions =
				Regions.getCentralRegions( labelImg, Transforms.getCenter( labelImg ), radius );

		final Img< BitType > maskFromLabelRegions = Algorithms.createMaskFromLabelRegions(
				centralRegions,
				Intervals.dimensionsAsLongArray( labelImg ) );

		return maskFromLabelRegions;
	}


	private RandomAccessibleInterval< BitType > createSpindleMask(
			RandomAccessibleInterval< T > poleToPoleAlignedSpindle,
			double spindleThreshold )
	{

		// TODO: do not consider pixels to zero that are not within the spindle range:
		// 1. further away from the spindle axis than the dnaLateralExtend/2
		// 2. further away from the center than the spindleAxialExtend/2

		final RandomAccessibleInterval< BitType > mask = Utils.createEmptyMask( poleToPoleAlignedSpindle );

		final Cursor< T > cursor = Views.iterable( poleToPoleAlignedSpindle ).localizingCursor();
		final RandomAccess< BitType > access = mask.randomAccess();

		long[] position = new long[ 3 ];

		final double maximumPerpendicularAxisDistanceSquared = Math.pow( dnaLateralExtend / 2, 2 );
		final double maximumAlongAxisDistance = 1.1 * spindleAxialExtend / 2;

		while ( cursor.hasNext() )
		{
			cursor.next();
			cursor.localize( position );

			// check position

			final double perpendicularAxisDistanceSquared =
					Math.pow( position[ 0 ] * settings.workingVoxelSize, 2 ) +
							Math.pow( position[ 1 ] * settings.workingVoxelSize, 2 );

			if ( perpendicularAxisDistanceSquared > maximumPerpendicularAxisDistanceSquared )
				continue;

			final double alongAxisDistance = Math.abs( position[ 2 ] * settings.workingVoxelSize );

			if ( alongAxisDistance > maximumAlongAxisDistance )
				continue;

			// check threshold

			if ( cursor.get().getRealDouble() > spindleThreshold )
			{
				access.setPosition( cursor );
				access.get().set( true );
			}

		}

		return mask;

//
//		RandomAccessibleInterval< BitType > mask = createMask( poleToPoleAlignedSpindle, threshold );
//
//		Regions.removeSmallRegionsInMask(
//				mask,
//				settings.minimalDnaFragementsVolume,
//				settings.workingVoxelSize );
//
//		final ImgLabeling< Integer, IntType > labelImg =
//				Utils.asImgLabeling( mask, ConnectedComponents.StructuringElement.FOUR_CONNECTED );
//
//		final long radius = (long) ( settings.maxCentralObjectRegionsDistance / settings.workingVoxelSize );
//
//		final Set< LabelRegion< Integer > > centralRegions =
//				Regions.getCentralRegions( labelImg, Transforms.getCenter( labelImg ), radius );
//
//		final Img< BitType > maskFromLabelRegions = Algorithms.createMaskFromLabelRegions(
//				centralRegions,
//				Intervals.dimensionsAsLongArray( labelImg ) );
//
//		return maskFromLabelRegions;
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
		spindleAxis[ SpindleMeasurements.ALIGNED_DNA_AXIS ] = 0;

		// compute rotation vector
		// note: we do not know in which direction of the spindle the spindleAxis vector points,
		// but we think it does not matter for alignment to the x-axis
		final double angleBetweenXAxisAndSpindleWithVector = angleOfSpindleAxisToXAxisInRadians( spindleAxis );

		AffineTransform3D rotationTransform = new AffineTransform3D();
		rotationTransform.rotate( SpindleMeasurements.ALIGNED_DNA_AXIS,  angleBetweenXAxisAndSpindleWithVector );

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

	class ProfileAndRadius
	{
		CoordinatesAndValues profile;
		Double radius;
		int radiusIndex;
	}

	private ProfileAndRadius measureRadialProfileAndRadius(
			final RandomAccessibleInterval< T > image,
			final String name,
			double derivativeDelta )
	{
		// config
		final double[] center = { 0, 0 };
		final double spacing = settings.workingVoxelSize;
		final double maxDistanceInMicrometer = 15;

		// measure
		final ProfileAndRadius profileAndRadius = new ProfileAndRadius();
		profileAndRadius.profile =
				Algorithms.computeRadialProfile( image, center, spacing, maxDistanceInMicrometer );

		final CoordinatesAndValues radialProfileDerivative =
				CurveAnalysis.derivative(
						profileAndRadius.profile,
						(int) Math.ceil( derivativeDelta / settings.workingVoxelSize ) );

		profileAndRadius.radius = CurveAnalysis.minimum( radialProfileDerivative ).coordinate;
		profileAndRadius.radiusIndex = (int) ( profileAndRadius.radius / spacing );

 		if ( settings.showIntermediateResults )
 			Plots.plot( profileAndRadius.profile,
					"center distance [um]", name + " intensity" );

 		if ( settings.showIntermediateResults )
			Plots.plot( radialProfileDerivative,
					"center distance [um]", "d/dx "+ name + " intensity" );

		return profileAndRadius;
	}

	public CompositeImage getOutputImage( )
	{
		RandomAccessibleInterval< T > interestPointsImage =
				createInterestPointImage( dna, dnaAxesAlignmentTransform, spindlePoles );

		createOutputImage(
				dna, tubulin,
				dnaVolumeMask, spindleVolumeMask,
				dnaAxesAlignmentTransform, interestPointsImage );

		RandomAccessibleInterval< T > image = getRandomAccessibleIntervalOutput();

		image = Views.permute( image, 2, 3 );

		final ImagePlus imp = ImageJFunctions.wrap( image, "output" );

		final Calibration calibration = new Calibration();
		calibration.pixelHeight = settings.workingVoxelSize;
		calibration.pixelWidth = settings.workingVoxelSize;
		calibration.pixelDepth = settings.workingVoxelSize;
		calibration.setUnit( "micrometer" );
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

	public RandomAccessibleInterval< T > getRandomAccessibleIntervalOutput()
	{
		final ArrayList< RandomAccessibleInterval< T > > list = new ArrayList<>();
		list.add( transformedDNAView );
		list.add( transformedDNAVolumeView );
		list.add( transformedSpindleView );
		list.add( transformedSpindleVolumeView );
		list.add( transformedInterestPointView );
		return Views.stack( list );
	}


	private void drawTransformedPoint(
			AffineTransform3D alignmentTransform,
			RandomAccessibleInterval< T > interestPointsImage,
			double[] point,
			int value )
	{
		final double[] transformedPoint = transformToPixelUnitsCoverslipCoordinateSystem( alignmentTransform, point );

		drawPoint(
				interestPointsImage,
				transformedPoint,
				settings.interestPointsRadius,
				settings.workingVoxelSize,
				value );
	}

	public HashMap< Integer, Map< String, Object > > getObjectMeasurements()
	{
		return objectMeasurements;
	}

	private RandomAccessibleInterval< BitType > createProcessedMetaPhasePlate( RandomAccessibleInterval< BitType > metaphasePlate )
	{
		Logger.log( "Perform morphological filtering on DNA mask..." );

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


	private double[] transformToPixelUnitsCoverslipCoordinateSystem(
			AffineTransform3D alignmentTransform,
			double[] point )
	{
		final double[] transformedPoint = Utils.copy( point );

		Utils.divide( transformedPoint, settings.workingVoxelSize ); // Scale to pixel units

		alignmentTransform.inverse().apply( transformedPoint, transformedPoint );

		return transformedPoint;
	}

	private void drawPoint(
			RandomAccessibleInterval< T > rai,
			double[] position,
			double calibratedRadius, double calibration, int value )
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
