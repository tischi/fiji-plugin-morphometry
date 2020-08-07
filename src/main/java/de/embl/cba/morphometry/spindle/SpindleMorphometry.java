package de.embl.cba.morphometry.spindle;

import bdv.util.AxisOrder;
import bdv.util.BdvFunctions;
import bdv.util.BdvOptions;
import de.embl.cba.cats.CATS;
import de.embl.cba.cats.results.ResultExportSettings;
import de.embl.cba.cats.results.ResultImage;
import de.embl.cba.morphometry.*;
import de.embl.cba.morphometry.geometry.CoordinatesAndValues;
import de.embl.cba.morphometry.geometry.CurveAnalysis;
import de.embl.cba.morphometry.geometry.ellipsoids.EllipsoidVectors;
import de.embl.cba.morphometry.geometry.ellipsoids.Ellipsoids3DImageSuite;
import de.embl.cba.morphometry.measurements.Measurements;
import de.embl.cba.morphometry.regions.Regions;
import de.embl.cba.morphometry.viewing.Viewers;
import de.embl.cba.neighborhood.RectangleShape2;
import de.embl.cba.transforms.utils.Transforms;
import ij.CompositeImage;
import ij.IJ;
import ij.ImagePlus;
import ij.measure.Calibration;
import net.imagej.ops.OpService;
import net.imglib2.*;
import net.imglib2.RandomAccess;
import net.imglib2.algorithm.gauss3.Gauss3;
import net.imglib2.algorithm.labeling.ConnectedComponents;
import net.imglib2.algorithm.neighborhood.HyperSphereShape;
import net.imglib2.algorithm.neighborhood.Neighborhood;
import net.imglib2.algorithm.neighborhood.Shape;
import net.imglib2.converter.Converters;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.interpolation.randomaccess.ClampingNLinearInterpolatorFactory;
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
import net.imglib2.type.numeric.real.FloatType;
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

public class SpindleMorphometry  < R extends RealType< R > & NativeType< R > >
{
	final SpindleMorphometrySettings< R > settings;
	final OpService opService;

	private HashMap< Integer, Map< String, Object > > objectMeasurements;
	private AffineTransform3D rescaledToDnaAlignmentTransform;
	private RandomAccessibleInterval dnaAlignedDna;
	private RandomAccessibleInterval dnaAlignedInitialDnaMask;
	private RandomAccessibleInterval dnaAlignedTubulin;
	private RandomAccessibleInterval< R > dna;
	private RandomAccessibleInterval< R > tubulin;
	private ArrayList< double[] > dnaAlignedSpindlePoles;
	private RandomAccessibleInterval< BitType > dnaAlignedDnaMask;
	private RandomAccessibleInterval< BitType > spindleAlignedSpindleMask;
	private double[] workingCalibration;
	private RandomAccessibleInterval< BitType > initialDnaMask;
	private EllipsoidVectors dnaEllipsoidVectors;
	private ProfileAndRadius dnaLateralProfileAndRadius;
	private double[] dnaAlignedSpindlePoleToPoleVector;
	private double[] dnaAlignedSpindleCenter;
	private RandomAccessibleInterval< R > spindleAlignedTublin;

	private SpindleMeasurements measurements;
	private RandomAccessibleInterval< R > raiXYCZ;
	private ArrayList< RandomAccessibleInterval< R > > rescaledVolumes;
	private RandomAccessibleInterval< BitType > spindleVolumeMask;
	private AffineTransform3D dnaAlignedToSpindleAlignedTransform;
	private ArrayList< double[] > spindleAlignedSpindlePoles;
	private RandomAccessibleInterval< R > spindleAlignedDna;
	private RandomAccessibleInterval< BitType > spindleAlignedDnaMask;
	private RandomAccessibleInterval< BitType > dnaAlignedSpindleMask;

	public SpindleMorphometry( SpindleMorphometrySettings settings, OpService opService )
	{
		this.settings = settings;
		this.opService = opService;
	}

	public String run( RandomAccessibleInterval<R> raiXYCZ )
	{
		this.raiXYCZ = raiXYCZ;

		objectMeasurements = new HashMap<>();

		measurements = new SpindleMeasurements( objectMeasurements );

		try
		{
			measurements.log = measure();
		}
		catch ( Exception e )
		{
			e.printStackTrace();
			measurements.log = "Exception during computation: \n" + e.toString();
		}

		measurements.setObjectMeasurements();

		return measurements.log;
	}

	public SpindleMeasurements getMeasurements()
	{
		return measurements;
	}

	private String measure()
	{
		measurements.version = settings.version;

		/**
		 * TODO:
		 * - maybe smooth the dna and tubulin image to reduce noise?
		 * - for low dynamic ranges, the smoothing should be done in a doubletype image
		 */

		createIsotropicallyResampledImages();

 //		if ( settings.useCATS )
//		{
//			createDnaMaskWithCATS();
//			return "";
//		}

		if ( settings.cellCenterDetectionMethod.equals( SpindleMorphometrySettings.CellCenterDetectionMethod.None ) )
		{
			// Do nothing
		}
		else
		{
			// TODO
//			final RealPoint cellCentreMicrometer = findCellCentre();
//			cropAroundCellCentre( cellCentreMicrometer );
		}

		measurements.dnaInitialThreshold = determineDnaThreshold();

		if ( measurements.dnaInitialThreshold < settings.minimalDynamicRange )
			return SpindleMeasurements.ANALYSIS_INTERRUPTED_LOW_DYNAMIC_DNA;

		initialDnaMask = segmentDna( dna, measurements.dnaInitialThreshold );

		dnaEllipsoidVectors = determineDnaAxes( initialDnaMask );

		rescaledToDnaAlignmentTransform = computeDnaAlignmentTransformAndAlignImages( dna, tubulin, initialDnaMask, dnaEllipsoidVectors );

		measureDnaAxialExtend( dnaAlignedDna );

		dnaLateralProfileAndRadius = measureDnaLateralExtend( dnaAlignedDna );

		dnaAlignedDnaMask = measureDnaVolume( dnaAlignedDna, dnaLateralProfileAndRadius );

		measureDnaHole( dnaLateralProfileAndRadius );

		measurements.spindleThreshold = measureSpindleThresholdAtDNAPeriphery( dnaAlignedTubulin, dnaAlignedDnaMask );

		dnaAlignedSpindleMask = createSpindleMask( dnaAlignedTubulin, measurements.spindleThreshold );

		// optional step to smoothen out the spindle outline
		//openFast( dnaAlignedSpindleMask );

		if ( settings.showIntermediateResults )
			show( dnaAlignedSpindleMask,
					"spindle volume mask", null,
					workingCalibration, false );


		// dnaAlignedSpindlePoles = measureTubulinIntensityBasedSpindlePoleLocations( dnaAlignedTubulin );

		dnaAlignedSpindlePoles = measureSpindlePoleLocationsFromSpindleMaskAndRefineWithTubulinIntensity( dnaAlignedSpindleMask, dnaAlignedTubulin );

		final ArrayList< double[] > rawAlignedSpindlePoles = new ArrayList<>();
		for ( double[] dnaAlignedSpindlePole : dnaAlignedSpindlePoles )
		{
			rawAlignedSpindlePoles.add( transformToCoverslipCoordinateSystem( rescaledToDnaAlignmentTransform, dnaAlignedSpindlePole) );
		}

		dnaAlignedSpindlePoleToPoleVector = getSpindleAxisVector( dnaAlignedSpindlePoles );

		dnaAlignedSpindleCenter = getMiddle( dnaAlignedSpindlePoles );

		measureDnaCenterToSpindleCenterDistance( dnaAlignedSpindleCenter );

		dnaAlignedToSpindleAlignedTransform = createSpindlePolesAlignedImages( dnaAlignedSpindlePoles, dnaAlignedSpindleCenter );

		if ( measurements.spindleThreshold < settings.minimalDynamicRange )
			return SpindleMeasurements.ANALYSIS_INTERRUPTED_LOW_DYNAMIC_TUBULIN;

		spindleAlignedSpindleMask = createSpindleMask(
				spindleAlignedTublin,
				measurements.spindleThreshold );

		measureSpindleVolume( spindleAlignedSpindleMask );

		measurements.spindleCoefficientOfVariation
				= Utils.measureCoefficientOfVariation(
						spindleAlignedTublin,
						spindleAlignedSpindleMask,
						measurements.spindleThreshold
				);

		measurements.spindleSum = Utils.measureSum(
				spindleAlignedTublin,
				spindleAlignedSpindleMask,
				measurements.spindleThreshold );

		measureSpindleLateralExtends( spindleAlignedSpindleMask );

		measureSpindleAxisToCoverslipPlaneAngle( dnaAlignedSpindlePoleToPoleVector );

		return SpindleMeasurements.ANALYSIS_FINISHED;
	}

	private RandomAccessibleInterval< BitType > openFast(
			RandomAccessibleInterval< BitType > dnaAlignedSpindleMask )
	{
		final RandomAccessibleInterval< BitType > downsampled = Algorithms.createNearestNeighborResampledArrayImg( dnaAlignedSpindleMask, new double[]{ 0.2, 0.2, 0.2 } );

		final RandomAccessibleInterval< BitType > opened = Algorithms.open( downsampled, 2 );

		final RandomAccessibleInterval< BitType > openedUpsampled = Algorithms.createNearestNeighborResampledArrayImg( opened, new double[]{ 1 / 0.2, 1 / 0.2, 1 / 0.2 } );

		return openedUpsampled;
	}

	private double[] getMiddle( ArrayList< double[] > vectors )
	{
		double[] middle = new double[ 3 ];
		LinAlgHelpers.add( vectors.get( 0 ), vectors.get( 1 ), middle );
		for ( int d = 0; d < 3; d++ )
		{
			middle[ d ] = 0.5 * middle[ d ];
		}
		return middle;
	}

	private void cropAroundCellCentre( RealPoint cellCentreMicrometer )
	{
		final FinalInterval cellBoundingBox = FinalInterval.createMinMax(
				getVoxelBoundingBoxMin( cellCentreMicrometer, 0 ), getVoxelBoundingBoxMin( cellCentreMicrometer, 1 ), getVoxelBoundingBoxMin( cellCentreMicrometer, 2 ),
				getVoxelBoundingBoxMax( cellCentreMicrometer, 0 ), getVoxelBoundingBoxMax( cellCentreMicrometer, 1 ), getVoxelBoundingBoxMax( cellCentreMicrometer, 2 )
		);

		final FinalInterval intersect = Intervals.intersect( cellBoundingBox, tubulin );

		tubulin = Views.zeroMin( Views.interval( tubulin, intersect ) );
		dna = Views.zeroMin( Views.interval( dna, intersect ) );

		Viewers.showRai3dWithImageJ( dna, "cropped dna" );
	}

	private int getVoxelBoundingBoxMin( RealPoint cellCentreMicrometer, int d )
	{
		return ( int ) ( cellCentreMicrometer.getDoublePosition( d ) / settings.workingVoxelSize - 1.5 * settings.cellRadius / settings.workingVoxelSize );
	}

	private int getVoxelBoundingBoxMax( RealPoint cellCentreMicrometer, int d )
	{
		return ( int ) ( cellCentreMicrometer.getDoublePosition( d ) / settings.workingVoxelSize + 1.5 * settings.cellRadius / settings.workingVoxelSize );
	}

	private RealPoint findCellCentre()
	{
		Logger.log( "Finding cell centre..." );

		RandomAccessibleInterval< R > cellCenterDetectionImage;
		switch ( settings.cellCenterDetectionMethod )
		{
			case BlurredTubulinImage:
				cellCenterDetectionImage = tubulin;
				break;
			case BlurredDnaImage:
				cellCenterDetectionImage = dna;
				break;
			default:
				cellCenterDetectionImage = dna;
				break;
		}

		RandomAccessibleInterval< R> blurred = new ArrayImgFactory( new FloatType( ) ).create( cellCenterDetectionImage );
		blurred = Transforms.getWithAdjustedOrigin( cellCenterDetectionImage, blurred );

		final double[] sigmas = new double[ 3 ];
		for ( int d = 0; d < 3; d++ )
		{
			sigmas[ d ] = settings.cellRadius / settings.workingVoxelSize;
		}

		// blur input image and write into output image
		Gauss3.gauss( sigmas, Views.extendBorder( cellCenterDetectionImage ), blurred ) ;

		final RealPoint cellCentre = Algorithms.getMaximumLocation( blurred, Utils.as3dDoubleArray( settings.workingVoxelSize ) );

		Viewers.showRai3dWithImageJ( blurred, "blurred cell centre detection image" );

		return cellCentre;
	}

	private void createDnaMaskWithCATS()
	{
		final CATS cats = new CATS();

		final ArrayList< RandomAccessibleInterval< R > > rais = new ArrayList< RandomAccessibleInterval< R > >();
		rais.add( tubulin );
		rais.add( dna );
		final ImagePlus wrap = ImageJFunctions.wrap( Views.permute( Views.stack( rais ), 2, 3 ), "" );

		//IJ.run( wrap, "8-bit", "");

		if ( settings.showMetaphaseClassification )
		{
			wrap.show();
		}

		if ( false )
		{

			cats.setInputImage( wrap );
			cats.setResultImageRAM();
			long memGB = 50;
			cats.setMaxMemoryBytes( memGB * 1024 * 1024 * 1024 );
			cats.setNumThreads( 4 );
			cats.loadClassifier( settings.classifier );
			cats.applyClassifierWithTiling();

			// configure results export
			final ResultExportSettings resultExportSettings = new ResultExportSettings();
			resultExportSettings.inputImagePlus = cats.getInputImage();
			resultExportSettings.exportType = ResultExportSettings.GET_AS_IMAGEPLUS_ARRAYLIST;
			resultExportSettings.exportNamesPrefix = "";
			resultExportSettings.classNames = cats.getClassNames();

			final ResultImage resultImage = cats.getResultImage();
			final ArrayList< ImagePlus > output = resultImage.exportResults( resultExportSettings );
			output.get( 0 ).show();
			output.get( 1 ).show();
		}

	}

	private void measureSpindleLateralExtends(
			RandomAccessibleInterval< BitType > alignedSpindleMask )
	{
		Logger.log( "Creating projection of spindle mask along spindle axis..." );

		RandomAccessibleInterval< BitType > projectedMask =
				new Projection<>(
					alignedSpindleMask,
					2 ).maximum();

		// remove spurious microtubules sticking out
		projectedMask = Algorithms.open( projectedMask, 2 );

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

		measurements.spindleWidthMin = calibratedMinMaxWidths[ 0 ];
		measurements.spindleWidthMax = calibratedMinMaxWidths[ 1 ];

	}

	public void measureSpindleVolume(
			RandomAccessibleInterval< BitType > spindleMask )
	{

		measurements.spindleVolume =
				Measurements.measureSizeInPixels( spindleMask )
						* Math.pow( settings.workingVoxelSize, 3 );
	}

	/**
	 *
	 * @param poleToPoleAlignedSpindleRai
	 * @return
	 */
	public double measureSpindleThresholdMethodProjectionRadial(
			RandomAccessibleInterval< R > poleToPoleAlignedSpindleRai )
	{
		final RandomAccessibleInterval< R > spindleProjection =
				createMaximumProjectionAlongSpindleAxis( poleToPoleAlignedSpindleRai );

		final ArrayList< Double > thresholdCandidates =
				measureRadialThresholds( spindleProjection );

		return Utils.median( thresholdCandidates ) ;
	}


	/**
	 *
	 * @param dnaAlignedTubulin
	 * @param dnaAlignedDnaMask
	 * @return
	 */
	public double measureSpindleThresholdAtDNAPeriphery(
			RandomAccessibleInterval< R > dnaAlignedTubulin,
			RandomAccessibleInterval< BitType > dnaAlignedDnaMask )
	{
		final double dnaRadius = measurements.dnaLateralExtend / 2.0 ;
		final long lateralExtend = (long) ( ( dnaRadius + 1.0 ) / settings.workingVoxelSize );
		final long axialExtend = (long) ( 1.0 / settings.workingVoxelSize );

		final FinalInterval interval = FinalInterval.createMinMax(
				-lateralExtend,
				-lateralExtend,
				-axialExtend,
				lateralExtend,
				lateralExtend,
				axialExtend );

		final Cursor< R > cursor = Views.iterable(
				Views.interval( dnaAlignedTubulin, interval ) ).cursor();

		double[] lateralVoxelPosition = new double[ 3 ];
		double distSquared = 0;
		final double minDistSquared = Math.pow( dnaRadius, 2);
		final double minInsideDistSquared = Math.pow( dnaRadius - 2.0, 2);
		final double maxDistSquared = Math.pow( dnaRadius + 1.0, 2);

		final ArrayList< Double > cytoplasmicTubulinValues = new ArrayList<>();
		final ArrayList< Double > spindleTubulinValues = new ArrayList<>();

		final RandomAccess< BitType > dnaAlignedDnaMaskAccess = dnaAlignedDnaMask.randomAccess();

		while( cursor.hasNext() )
		{
			cursor.next();

			dnaAlignedDnaMaskAccess.setPosition( cursor );

			if ( dnaAlignedDnaMaskAccess.get().get() )
			{
				// pixels containing DNA exclude tubulin and thus would
				// lead to a too low threshold
				continue;
			}

			cursor.localize( lateralVoxelPosition );

			distSquared = 0;
			for ( int d = 0; d < 2; d++ )
				distSquared += Math.pow( lateralVoxelPosition[ d ] * settings.workingVoxelSize, 2 );

			if ( distSquared < minInsideDistSquared )
			{
				spindleTubulinValues.add( cursor.get().getRealDouble() );
				continue;
			}

			if ( distSquared > maxDistSquared ) continue;
			if ( distSquared < minDistSquared ) continue;

			cytoplasmicTubulinValues.add( cursor.get().getRealDouble() );
		}

		final double meanCytoplasmic = Utils.mean( cytoplasmicTubulinValues );
		final double sdev = Utils.sdev( cytoplasmicTubulinValues, meanCytoplasmic );

		final double meanInsideSpindle = Utils.mean( spindleTubulinValues );

		final double threshold = ( meanCytoplasmic + meanInsideSpindle ) / 2;

		Logger.log( "Spindle threshold: " + threshold );

		return threshold;
	}

	private RandomAccessibleInterval< R > createMaximumProjectionAlongSpindleAxis(
			RandomAccessibleInterval< R > poleToPoleAlignedSpindleRai )
	{
		Logger.log( "Computing maximum projection of spindle along spindle axis..." );

		final FinalInterval spindleInterval = createSpindleInterval();

		Projection projection = new Projection<>(
				Views.interval( poleToPoleAlignedSpindleRai, spindleInterval ),
				2 );

		final RandomAccessibleInterval< R > maximum = projection.maximum();

		if ( settings.showIntermediateResults )
			show( maximum, "Spindle maximum projection along pole to pole axis",
					null, workingCalibration, false);

		return maximum;
	}

	private FinalInterval createSpindleInterval()
	{
		long[] min = new long[3];
		long[] max = new long[3];

		max[ 0 ] = (long) Math.ceil( measurements.dnaLateralExtend
				/ 2.0 / settings.workingVoxelSize );
		max[ 1 ] = (long) Math.ceil( measurements.dnaLateralExtend
				/ 2.0 / settings.workingVoxelSize );
		max[ 2 ] = (long) Math.ceil( 1.2 * measurements.spindleAxialExtend
				/ 2.0 / settings.workingVoxelSize );

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
	private AffineTransform3D createSpindlePolesAlignedImages( ArrayList< double[] > spindlePoles, double[] spindleCenter )
	{
		dnaAlignedToSpindleAlignedTransform = computeTransform( spindlePoles, spindleCenter );

		spindleAlignedTublin = Transforms.createTransformedView( dnaAlignedTubulin, dnaAlignedToSpindleAlignedTransform );
		spindleAlignedDna = Transforms.createTransformedView( dnaAlignedDna, dnaAlignedToSpindleAlignedTransform );
		spindleAlignedDnaMask = Transforms.createTransformedView( dnaAlignedInitialDnaMask, dnaAlignedToSpindleAlignedTransform, new NearestNeighborInterpolatorFactory() );

		spindleAlignedSpindlePoles = new ArrayList<>();
		for ( int i = 0; i < 2; i++ )
		{
			final double[] transformed = transformed( spindlePoles.get( i ), dnaAlignedToSpindleAlignedTransform );
			spindleAlignedSpindlePoles.add( transformed );
		}

		if ( settings.showIntermediateResults )
		{
			final ArrayList< RealPoint > interestPoints = new ArrayList<>();

			for ( int i = 0; i < 2; i++ )
			{
				interestPoints.add( new RealPoint( spindleAlignedSpindlePoles.get( i ) ) );
			}
			interestPoints.add( new RealPoint( 0,0,0 ) );

			show( spindleAlignedTublin,
					"spindle aligned pole to pole",
					interestPoints,
					workingCalibration,
					false );
		}

		return dnaAlignedToSpindleAlignedTransform;
	}

	private double[] transformed( double[] location, AffineTransform3D affineTransform3D )
	{
		final double[] transformed = new double[ 3 ];
		// transformation is in voxel units
		final double[] voxels = toVoxels( location );
		affineTransform3D.apply( voxels, transformed );
		final double[] calibrated = toCalibrated( transformed );
		return calibrated;
	}

	private AffineTransform3D computeTransform( ArrayList< double[] > spindlePoles, double[] spindleCenter )
	{
		final double[] poleToPoleAxis = new double[ 3 ];
		LinAlgHelpers.subtract( spindlePoles.get( 0 ), spindlePoles.get( 1 ), poleToPoleAxis );
		LinAlgHelpers.normalize( poleToPoleAxis );

		// put spindle center in the center
		final double[] spindleCenterVoxelUnits = toVoxels( spindleCenter );
		AffineTransform3D dnaAlignedToSpindleAlignedTransform = new AffineTransform3D();
		dnaAlignedToSpindleAlignedTransform.translate( spindleCenterVoxelUnits );
		dnaAlignedToSpindleAlignedTransform = dnaAlignedToSpindleAlignedTransform.inverse();

		// rotate spindle along z-axis
		AffineTransform3D poleToPoleAxisRotation =
				Transforms.getRotationTransform3D( new double[]{ 0, 0, 1 }, poleToPoleAxis );
		dnaAlignedToSpindleAlignedTransform.preConcatenate( poleToPoleAxisRotation );

		return dnaAlignedToSpindleAlignedTransform;
	}

	private double[] toVoxels( double[] location )
	{
		final double[] spindleCenterVoxelUnits = location.clone();
		Utils.divide( spindleCenterVoxelUnits, settings.workingVoxelSize );
		return spindleCenterVoxelUnits;
	}

	private double[] toCalibrated( double[] location )
	{
		final double[] spindleCenterVoxelUnits = location.clone();
		Utils.divide( spindleCenterVoxelUnits, 1.0 / settings.workingVoxelSize );
		return spindleCenterVoxelUnits;
	}

	@Deprecated
	public void measureRadialWidthsOld( RandomAccessibleInterval< R > rai )
	{
		RealRandomAccessible< R > rra =
				Views.interpolate(
						Views.extendZero( rai ), new NLinearInterpolatorFactory<>() );

		double dAngle = Math.PI / 3;

		for ( double angle = 0; angle < Math.PI; angle += dAngle )
		{
			final AffineTransform2D transform2D = new AffineTransform2D();
			transform2D.rotate( angle );

			IntervalView< R > rotated =
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

	public ArrayList< Double > measureRadialThresholds( RandomAccessibleInterval< R > rai )
	{
		final ArrayList< Double > thresholdCandidates = new ArrayList<>();

		RealRandomAccessible< R > rra =
				Views.interpolate(
						Views.extendZero( rai ),
						new ClampingNLinearInterpolatorFactory<>() );

		double dAngle = Math.PI / 5;

		for ( double angle = 0; angle < Math.PI; angle += dAngle )
		{
			final AffineTransform2D transform2D = new AffineTransform2D();
			transform2D.rotate( angle );

			IntervalView< R > rotated =
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

			final double leftLoc = extrema.get( 0 ).coordinate
					- settings.spindleDerivativeDelta / 2.0;

			thresholdCandidates.add(
					CurveAnalysis.getValueAtCoordinate(
						intensities,
						leftLoc ) );

			final double rightLoc = extrema.get( 1 ).coordinate
					+ settings.spindleDerivativeDelta / 2.0;

			thresholdCandidates.add(
					CurveAnalysis.getValueAtCoordinate(
						intensities,
						rightLoc ) );

			if ( settings.showIntermediateResults )
			{
//				Plots.plot( intensities, "Center distance [um]",
//						"Spindle lateral intensity, angle = " + angle  );
//				Plots.plot( derivatives, "Center distance [um]",
//						"d/dx Spindle lateral intensity, angle = " + angle  );
			}

		}

		return thresholdCandidates;
	}

	public void measureSpindleAxisToCoverslipPlaneAngle( double[] poleToPoleVector )
	{
		final double[] poleToPoleVectorInCSCS =
				transformToCoverslipCoordinateSystem( rescaledToDnaAlignmentTransform, poleToPoleVector );

		measurements.angleSpindleAxisToCoverslipPlaneInDegrees =
				90.0 - Math.abs( 180.0 / Math.PI *
						Transforms.getAngle( new double[]{ 0, 0, 1 }, poleToPoleVectorInCSCS ) );
	}

	public void measureDnaCenterToSpindleCenterDistance( double[] spindleCenter )
	{
		final double[] dnaCentre = { 0, 0, 0 };
		measurements.dnaCenterToSpindleCenterDistance = LinAlgHelpers.distance( dnaCentre, spindleCenter);
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

		long numChannels = raiXYCZ.dimension( 2 );

		rescaledVolumes = new ArrayList<>();

		final double[] scalingFactors = getScalingFactors( settings.inputCalibration, settings.workingVoxelSize );

		for ( int c = 0; c < numChannels; c++ )
		{
			final IntervalView< R > volume = Views.hyperSlice( raiXYCZ, 2, c );
			final RandomAccessibleInterval< R > rescaledVolume =
					createRescaledArrayImg(
						volume,
						scalingFactors );
			rescaledVolumes.add( rescaledVolume );
		}

		dna = rescaledVolumes.get( (int) settings.dnaChannelIndex );
		tubulin = rescaledVolumes.get( (int) settings.tubulinChannelIndex );

		if ( settings.showIntermediateResults )
			show( dna, "DNA isotropic voxel size", null, workingCalibration, false );

		if ( settings.showIntermediateResults )
			show( tubulin, "spindle isotropic voxel size", null, workingCalibration, false );
	}

	public double determineDnaThreshold()
	{
		final RandomAccessibleInterval< R > dnaDownscaledToMetaphasePlateWidth =
				createRescaledArrayImg(
						dna,
						getScalingFactors(
								new double[]{
										settings.workingVoxelSize,
										settings.workingVoxelSize,
										settings.workingVoxelSize },
								settings.dnaThresholdResolution ) );

		final double maximumValue = Algorithms.getMaximumValue( dnaDownscaledToMetaphasePlateWidth );

//		Viewers.showRai3dWithImageJ( dnaDownscaledToMetaphasePlateWidth, "DNA Threshold" );

//		Logger.log( "DNA downscaled maximum value: " + maximumValue );
//		Logger.log( "DNA threshold factor: " + settings.dnaThresholdFactor );
		double dnaThreshold = maximumValue * settings.dnaThresholdFactor;
		Logger.log( "DNA initial threshold: " + dnaThreshold );

		return dnaThreshold;
	}

	public RandomAccessibleInterval< BitType > segmentDna(
			RandomAccessibleInterval< R > dna,
			double dnaThreshold )
	{
		final RandomAccessibleInterval< BitType > dnaMask = createCentralObjectsMask( dna, dnaThreshold );

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

	public AffineTransform3D computeDnaAlignmentTransformAndAlignImages(
			RandomAccessibleInterval< R > dna,
			RandomAccessibleInterval< R > tubulin,
			RandomAccessibleInterval< BitType > dnaMask,
			EllipsoidVectors ellipsoidVectors )
	{
		Logger.log( "Creating aligned images..." );

		rescaledToDnaAlignmentTransform =
				Ellipsoids3DImageSuite.createShortestAxisAlignmentTransform( ellipsoidVectors );

		dnaAlignedTubulin = Utils.copyAsArrayImg(
				Transforms.createTransformedView( tubulin,
						rescaledToDnaAlignmentTransform, new NearestNeighborInterpolatorFactory() ) );

		dnaAlignedDna = Utils.copyAsArrayImg(
				Transforms.createTransformedView( dna, rescaledToDnaAlignmentTransform ) );

		dnaAlignedInitialDnaMask = Utils.copyAsArrayImg(
				Transforms.createTransformedView( dnaMask, rescaledToDnaAlignmentTransform ) );

		if ( settings.showIntermediateResults )
			show( dnaAlignedTubulin, "tubulin aligned along DNA shortest axis", Transforms.origin(),
					workingCalibration, false );

		if ( settings.showIntermediateResults )
			show( dnaAlignedDna, "DNA aligned along DNA shortest axis", Transforms.origin(),
					workingCalibration, false );

		if ( settings.showIntermediateResults )
			show( dnaAlignedInitialDnaMask, "DNA mask aligned along DNA shortest axis",
					Transforms.origin(), workingCalibration, false );

		return rescaledToDnaAlignmentTransform;
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
						settings.maxDnaLateralRadius,
						2,
						settings.workingVoxelSize );

		final CoordinatesAndValues dnaProfileAlongDnaAxisDerivative =
				CurveAnalysis.derivative(
						dnaProfileAlongDnaAxis,
						( int ) Math.ceil( settings.derivativeDelta
								/ settings.workingVoxelSize ) );

		final ArrayList< CoordinateAndValue > dnaAxialBoundaries =
				CurveAnalysis.leftMaxAndRightMinLoc( dnaProfileAlongDnaAxisDerivative );

		measurements.dnaAxialExtend =
				dnaAxialBoundaries.get( 1 ).coordinate -
						dnaAxialBoundaries.get( 0 ).coordinate;

		if ( settings.showIntermediateResults )
			Plots.plot( dnaProfileAlongDnaAxis,
					"distance to center", "DNA intensity along DNA axis" );

		if ( settings.showIntermediateResults )
			Plots.plot( dnaProfileAlongDnaAxisDerivative,
					"distance to center", "d/dx DNA intensity along DNA axis" );

		logDnaAxialThreshold( dnaProfileAlongDnaAxis, dnaAxialBoundaries );
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
			RandomAccessibleInterval< R > alignedDNA )
	{
		Logger.log( "Measuring DNA lateral extend..." );

		Projection projection = new Projection( alignedDNA, 2 );

		final RandomAccessibleInterval< R > dnaProjectionAlongDnaAxis = projection.maximum();

		final ProfileAndRadius dnaLateralProfileAndRadius =
				measureRadialProfileAndRadius(
						dnaProjectionAlongDnaAxis,
						"dna lateral",
						settings.derivativeDelta );

		measurements.dnaLateralExtend = 2.0 * dnaLateralProfileAndRadius.radius;

		return dnaLateralProfileAndRadius;
	}

	/**
	 * Find the spindle poles
	 *
	 * @param dnaAlignedSpindleMask
	 * @param dnaAlignedTubulin
	 * @return
	 */
	public ArrayList< double[] > measureSpindlePoleLocationsFromSpindleMaskAndRefineWithTubulinIntensity(
			final RandomAccessibleInterval< BitType > dnaAlignedSpindleMask,
			final RandomAccessibleInterval< R > dnaAlignedTubulin )
	{
		final ArrayList< double[] > spindlePoles = determineSpindlePolesAlongDnaAxisFromSpindleMask( dnaAlignedSpindleMask );

		final double distance = LinAlgHelpers.distance( spindlePoles.get( 0 ), spindlePoles.get( 1 ) );
		System.out.println("[DEBUG] Spindle_length: " + distance);

		final ArrayList< double[] > refinedSpindlePoles = refineSpindlePoles( dnaAlignedTubulin, spindlePoles, settings.maxSpindlePoleRefinementDistance * 2 );

		final double distanceRefined = LinAlgHelpers.distance( refinedSpindlePoles.get( 0 ), refinedSpindlePoles.get( 1 ) );
		System.out.println("[DEBUG] Refined Spindle_length: " + distanceRefined);

		addPoleRefinementDistanceMeasurements( spindlePoles, refinedSpindlePoles );

		measurements.spindleAxialExtend = LinAlgHelpers.distance( refinedSpindlePoles.get( 0 ), refinedSpindlePoles.get( 1 ) );

		if ( settings.showIntermediateResults )
		{
			final ArrayList< RealPoint > interestPoints = new ArrayList<>();

			for ( int i = 0; i < 2; i++ )
			{
				interestPoints.add( new RealPoint( spindlePoles.get( i ) ) );
			}

			show( dnaAlignedTubulin,
					"dna aligned tubulin with spindle poles",
					interestPoints,
					workingCalibration,
					false );
		}

		if ( settings.showIntermediateResults )
		{
			final ArrayList< RealPoint > interestPoints = new ArrayList<>();

			for ( int i = 0; i < 2; i++ )
			{
				interestPoints.add( new RealPoint( refinedSpindlePoles.get( i ) ) );
			}

			show( dnaAlignedTubulin,
					"dna aligned tubulin with refined spindle poles",
					interestPoints,
					workingCalibration,
					false );
		}


		return refinedSpindlePoles;
	}

	public RandomAccessibleInterval< BitType > measureDnaVolume(
			RandomAccessibleInterval< R > dna,
			ProfileAndRadius dnaLateralExtendAndProfile )
	{
		measurements.dnaVolumeThreshold =
				dnaLateralExtendAndProfile.profile.values.get(
						dnaLateralExtendAndProfile.radiusIndex );

		final RandomAccessibleInterval< BitType > dnaVolumeMask =
				createCentralObjectsMask( dna, measurements.dnaVolumeThreshold );

		final long dnaVolumeInPixels =
				Measurements.measureSizeInPixels( dnaVolumeMask );

		measurements.dnaVolumeCalibrated =
				dnaVolumeInPixels * Math.pow( settings.workingVoxelSize, 3 );

		if ( settings.showIntermediateResults )
			show( dnaVolumeMask, "dna volume mask",
					null, workingCalibration, false );

		return dnaVolumeMask;
	}

	public ArrayList< double[] > determineSpindlePolesAlongDnaAxisFromTubulinIntensity(
			RandomAccessibleInterval tubulinAlignedAlongShortestAxisOfDNA )
	{
		/**
		 * Taking the maximum here can be quite noisy, depending on the microscopy
		 * TODO: Maybe we should smooth the image before doing this?
		 */
		final CoordinatesAndValues tubulinProfile =
				Utils.computeMaximumIntensitiesAlongAxis(
						tubulinAlignedAlongShortestAxisOfDNA,
						settings.maxDnaLateralRadius,
						SpindleMeasurements.ALIGNED_DNA_AXIS,
						settings.workingVoxelSize );

		final int derivativeDeltaVoxels = ( int ) Math.ceil( settings.derivativeDelta / settings.workingVoxelSize );
		final CoordinatesAndValues tubulinProfileDerivative =
				CurveAnalysis.derivative(
						tubulinProfile,
						derivativeDeltaVoxels );

		ArrayList< CoordinateAndValue > tubulinExtrema =
				CurveAnalysis.leftMaxAndRightMinLoc( tubulinProfileDerivative );

		// the actual spindle poles are some
		// distance more "inside" from the derivative extrema,
		// partly due to blur by diffraction limit
		double offset = 0.5; // micrometer

		final ArrayList< double[] > dnaAxisBasedSpindlePoles = new ArrayList<>();
		dnaAxisBasedSpindlePoles.add(
				new double[]{ 0, 0, tubulinExtrema.get( 0 ).coordinate + offset });
		dnaAxisBasedSpindlePoles.add(
				new double[]{ 0, 0, tubulinExtrema.get( 1 ).coordinate - offset });

		if ( settings.showIntermediateResults )
			Plots.plot(
					tubulinProfile.coordinates,
					tubulinProfile.values,
					"center distance [um]",
					"spindle max along shortest DNA axis" );

		if ( settings.showIntermediateResults )
		{
			Plots.plot(
					tubulinProfileDerivative.coordinates,
					tubulinProfileDerivative.values,
					"distance to center",
					"d/dx spindle max along shortest DNA axis" );
			Logger.log( "Spindle max derivative plus offset " + ( tubulinExtrema.get( 0 ).coordinate + offset ) );
			Logger.log( "Spindle max derivative minus offset " + ( tubulinExtrema.get( 1 ).coordinate - offset ) );
		}

		return dnaAxisBasedSpindlePoles;
	}

	public ArrayList< double[] > determineSpindlePolesAlongDnaAxisFromSpindleMask( RandomAccessibleInterval< BitType > dnaAlignedSpindleMask )
	{
		final CoordinatesAndValues tubulinProfile =
				Utils.computeMaximumIntensitiesAlongAxis(
						dnaAlignedSpindleMask,
						settings.maxDnaLateralRadius, // TODO: replace by measured DNA width?!
						SpindleMeasurements.ALIGNED_DNA_AXIS,
						settings.workingVoxelSize );

		final int derivativeDeltaVoxels = 2; // ( int ) Math.ceil( settings.derivativeDelta / settings.workingVoxelSize );
		final CoordinatesAndValues tubulinProfileDerivative =
				CurveAnalysis.derivative(
						tubulinProfile,
						derivativeDeltaVoxels );

		ArrayList< CoordinateAndValue > tubulinExtrema =
				CurveAnalysis.leftMaxAndRightMinLoc( tubulinProfileDerivative );

		double offset = 0.0;

		final ArrayList< double[] > dnaAxisBasedSpindlePoles = new ArrayList<>();
		dnaAxisBasedSpindlePoles.add(
				new double[]{ 0, 0, tubulinExtrema.get( 0 ).coordinate + offset });
		dnaAxisBasedSpindlePoles.add(
				new double[]{ 0, 0, tubulinExtrema.get( 1 ).coordinate - offset });

		if ( settings.showIntermediateResults )
			Plots.plot(
					tubulinProfile.coordinates,
					tubulinProfile.values,
					"center distance [um]",
					"spindle mask avg along shortest DNA axis" );

		if ( settings.showIntermediateResults )
		{
			Plots.plot(
					tubulinProfileDerivative.coordinates,
					tubulinProfileDerivative.values,
					"distance to center",
					"d/dx spindle max along shortest DNA axis" );
			Logger.log( "Spindle max derivative plus offset " + ( tubulinExtrema.get( 0 ).coordinate + offset ) );
			Logger.log( "Spindle max derivative minus offset " + ( tubulinExtrema.get( 1 ).coordinate - offset ) );
		}

		return dnaAxisBasedSpindlePoles;
	}

	public void measureDnaHole( ProfileAndRadius dnaLateralProfileAndRadius )
	{
		final double dnaLateralRadialProfileMaxIntensity =
				CurveAnalysis.maximumIndexAndValue( dnaLateralProfileAndRadius.profile ).value;

		final double dnaCenterMaxIntensity = dnaLateralProfileAndRadius.profile.values.get( 0 );

		measurements.dnaRelativeCentralIntensity =
				dnaCenterMaxIntensity / dnaLateralRadialProfileMaxIntensity;

	}

	public void addPoleRefinementDistanceMeasurements(
			ArrayList< double[] > poles,
			ArrayList< double[] > refinedPoles )
	{
		measurements.spindlePoleARefinementDistance =
				LinAlgHelpers.distance(
						poles.get( 0 ), refinedPoles.get( 0 ) );

		measurements.spindlePoleBRefinementDistance =
				LinAlgHelpers.distance(
						poles.get( 1 ), refinedPoles.get( 1 ) );
	}

	/**
	 * Refine spindle pole positions
	 *
	 * @param dnaAlignedTubulin
	 * @param dnaAlignedSpindlePoles
	 * @param searchRadiusAlongSpindleAxis
	 * @return
	 */
	private ArrayList< double[] > refineSpindlePoles(
			final RandomAccessibleInterval< R > dnaAlignedTubulin,
			final ArrayList< double[] > dnaAlignedSpindlePoles,
			double searchRadiusAlongSpindleAxis )
	{
		final double searchFullWidthPerpendicularToSpindleAxis = 4.0; // micrometer
		final double blurSigma = 0.75; // micrometer

		final RandomAccessibleInterval< R > blurred = Utils.createBlurredRai(
				dnaAlignedTubulin,
				blurSigma,
				settings.workingVoxelSize );

		// wide perpendicular to spindle axis, narrow along spindle axis
		final RectangleShape2 rectangleShape2 = new RectangleShape2( new long[]{
				( long ) ( searchFullWidthPerpendicularToSpindleAxis / settings.workingVoxelSize ),
				( long ) ( searchFullWidthPerpendicularToSpindleAxis / settings.workingVoxelSize ),
				( long ) ( searchRadiusAlongSpindleAxis / settings.workingVoxelSize )
		}, false );

		final RandomAccessible< Neighborhood < R > > neighborhoodsAccessible = rectangleShape2.neighborhoodsRandomAccessible( blurred );

		final RandomAccess< Neighborhood< R > > access = neighborhoodsAccessible.randomAccess();

		final ArrayList< double[] > spindlePoles = new ArrayList<>();

		for ( int pole = 0; pole < 2; pole++ )
		{
			final long[] pixelUnitsPolePosition = getPixelUnitsPolePosition( dnaAlignedSpindlePoles, pole );
			access.setPosition( pixelUnitsPolePosition );
			final Neighborhood< R > neighborhood = access.get();

			final RealPoint maximumLocation = Algorithms.getMaximumLocation(
					neighborhood,
					Utils.as3dDoubleArray( settings.workingVoxelSize ) );

			spindlePoles.add( new double[]{
					maximumLocation.getDoublePosition( 0 ),
					maximumLocation.getDoublePosition( 1 ),
					maximumLocation.getDoublePosition( 2 )} );
		}

		return spindlePoles;
	}

	private long[] getPixelUnitsPolePosition( ArrayList< double[] > dnaAlignedSpindlePoles, int pole )
	{
		final double[] polePosition = toVoxels( dnaAlignedSpindlePoles.get( pole ) );
		return Utils.asLongs( polePosition );
	}

	private RandomAccessibleInterval< BitType > createCentralObjectsMask(
			RandomAccessibleInterval< R > image,
			double threshold )
	{
		RandomAccessibleInterval< BitType > mask = createMask( image, threshold );

		Regions.removeSmallRegionsInMask(
				mask,
				settings.minimalDnaFragmentsVolume,
				settings.workingVoxelSize );

		final ImgLabeling< Integer, IntType > labelImg =
				Regions.asImgLabeling( mask, ConnectedComponents.StructuringElement.FOUR_CONNECTED );

		final long radius = (long) ( settings.maxCentralObjectRegionsDistance / settings.workingVoxelSize );

		final Set< LabelRegion< Integer > > centralRegions =
				Regions.getCentralRegions( labelImg, Transforms.getCenter( labelImg ), radius );

		final RandomAccessibleInterval< BitType > maskFromLabelRegions =
				Regions.asMask(
					centralRegions,
					labelImg );

		return maskFromLabelRegions;
	}


	private RandomAccessibleInterval< BitType > createSpindleMask(
			RandomAccessibleInterval< R > tubulin,
			double spindleThreshold )
	{
		final RandomAccessibleInterval< BitType > mask =
				Algorithms.createMask( tubulin, spindleThreshold );

		final ImgLabeling< Integer, IntType > imgLabeling
				= Regions.asImgLabeling( mask, ConnectedComponents.StructuringElement.FOUR_CONNECTED );

		final Set< LabelRegion< Integer > > centralRegions = Regions.getCentralRegions(
				imgLabeling,
				new double[]{ 0, 0, 0 },
				( long ) ( 3.0 / settings.workingVoxelSize ) );

		final RandomAccessibleInterval< BitType > centralRegionsMask =
				Regions.asMask( centralRegions,
						Intervals.dimensionsAsLongArray( mask ),
						Intervals.minAsLongArray( mask ));

		return centralRegionsMask;
	}


	private RandomAccessibleInterval< BitType > createSpindleMask00(
			RandomAccessibleInterval< R > poleToPoleAlignedSpindle,
			double spindleThreshold )
	{

		final RandomAccessibleInterval< BitType > mask =
				Utils.createEmptyMask( poleToPoleAlignedSpindle );

		final Cursor< R > cursor =
				Views.iterable( poleToPoleAlignedSpindle ).localizingCursor();

		final RandomAccess< BitType > access = mask.randomAccess();

		long[] position = new long[ 3 ];

		final double maximumPerpendicularAxisDistanceSquared = Math.pow(
				measurements.dnaLateralExtend / 2.0, 2 );
		final double maximumAlongAxisDistance = 1.1 *
				measurements.spindleAxialExtend / 2.0;

		while ( cursor.hasNext() )
		{
			cursor.next();
			cursor.localize( position );

			final double perpendicularAxisDistanceSquared =
					Math.pow( position[ 0 ] * settings.workingVoxelSize, 2 ) +
							Math.pow( position[ 1 ] * settings.workingVoxelSize, 2 );

			if ( perpendicularAxisDistanceSquared > maximumPerpendicularAxisDistanceSquared )
				continue;

			final double alongAxisDistance =
					Math.abs( position[ 2 ] * settings.workingVoxelSize );

			if ( alongAxisDistance > maximumAlongAxisDistance )
				continue;

			if ( cursor.get().getRealDouble() > spindleThreshold )
			{
				access.setPosition( cursor );
				access.get().set( true );
			}

		}
		return mask;
	}

	private RandomAccessibleInterval< R > createTransformedMultiChannelOutputImage(
			RandomAccessibleInterval< BitType > spindleAlignedDnaMask,
			RandomAccessibleInterval< BitType > spindleAlignedSpindleMask,
			AffineTransform3D rescaledInputToSpindleAlignedTransform,
			RandomAccessibleInterval< R > spindleAlignedInterestPointsImage,
			FinalInterval crop )
	{
		final ArrayList< RandomAccessibleInterval< R > > alignedVolumes = new ArrayList<>();

		for ( int c = 0; c < rescaledVolumes.size(); c++ )
		{
			final RandomAccessible transformedRA = Transforms.createTransformedRaView( rescaledVolumes.get( c ),
					rescaledInputToSpindleAlignedTransform, new ClampingNLinearInterpolatorFactory() );

			alignedVolumes.add( Views.interval( transformedRA, crop ) );
		}

		alignedVolumes.add( (RandomAccessibleInterval) Views.interval( Views.extendValue( spindleAlignedDnaMask, new BitType( false ) ), crop ) );

		alignedVolumes.add( (RandomAccessibleInterval) Views.interval( Views.extendValue( spindleAlignedSpindleMask, new BitType( false ) ), crop )  );

		alignedVolumes.add( (RandomAccessibleInterval) Views.interval( Views.extendZero( spindleAlignedInterestPointsImage ), crop ) );


		// The volumes are now aligned such that the spindle is along the z-axis
		// We rotate it now to be along the x-axis to make it easier to view in ImageJ
		final ArrayList< RandomAccessibleInterval< R > > alignedAndRotatedVolumes = new ArrayList<>();
		for ( RandomAccessibleInterval< R > rai : alignedVolumes )
			alignedAndRotatedVolumes.add( Views.rotate( rai, 0, 2 ) );

		return Views.stack( alignedAndRotatedVolumes );
	}

	/**
	 *
	 *
	 * @param alignmentTransform
	 * @return
	 */
	private AffineTransform3D createOutputImageAffineTransform3D( AffineTransform3D alignmentTransform )
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
		rotationTransform.rotate( SpindleMeasurements.ALIGNED_DNA_AXIS, angleBetweenXAxisAndSpindleWithVector );

		Logger.log( "Rotating input data around z-axis by [degrees]: " + 180 / Math.PI * angleBetweenXAxisAndSpindleWithVector );

		return rotationTransform;
	}

	private RandomAccessibleInterval< R > createInterestPointImage(
			ArrayList< double[] > spindlePoles,
			RandomAccessibleInterval< R > spindleAlignedDna )
	{
		RandomAccessibleInterval< R > interestPointsImage = Utils.copyAsArrayImg( spindleAlignedDna );
		Utils.setValues( interestPointsImage, 0.0 );

		final double[] origin = { 0, 0, 0 };
		drawPoint(
				interestPointsImage,
				origin,
				settings.interestPointsRadius,
				settings.workingVoxelSize,
				1 );

		final double distance = LinAlgHelpers.distance( spindlePoles.get( 0 ), spindlePoles.get( 1 ) );

		for ( double[] spindlePole : spindlePoles )
			drawPoint(
					interestPointsImage,
					spindlePole,
					settings.interestPointsRadius,
					settings.workingVoxelSize,
					1 );

		return interestPointsImage;
	}



	class ProfileAndRadius
	{
		CoordinatesAndValues profile;
		Double radius;
		int radiusIndex;
	}

	private ProfileAndRadius measureRadialProfileAndRadius(
			final RandomAccessibleInterval< R > image,
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

	public CompositeImage createOutputImage()
	{
		final AffineTransform3D rescaledInputToSpindleAlignedTransform =
				rescaledToDnaAlignmentTransform.preConcatenate( dnaAlignedToSpindleAlignedTransform );

		final int size = (int) ( settings.cellRadius / settings.workingVoxelSize );

		final int[] center = { 0, 0, 0 };
		final FinalInterval crop = FinalInterval.createMinMax(
				center[ 0 ] - 3 * size,
				center[ 1 ] - 3 * size,
				center[ 2 ] - 3 * size, // spindle axis
				center[ 0 ] + 3 * size,
				center[ 1 ] + 3 * size,
				center[ 2 ] + 3 * size // spindle axis
		);

		RandomAccessibleInterval< R > raiXYZC = createTransformedMultiChannelOutputImage(
				spindleAlignedDnaMask,
				spindleAlignedSpindleMask,
				rescaledInputToSpindleAlignedTransform,
				createInterestPointImage( spindleAlignedSpindlePoles, spindleAlignedDna ),
				crop );

		if ( settings.showIntermediateResults )
		{
			BdvFunctions.show( raiXYZC, "output", BdvOptions.options().axisOrder( AxisOrder.XYZC ) );
		}

		final ImagePlus imp = ImageJFunctions.wrap( Views.permute( raiXYZC, 2, 3 ), "output" );

		final Calibration calibration = new Calibration();
		calibration.pixelHeight = settings.workingVoxelSize;
		calibration.pixelWidth = settings.workingVoxelSize;
		calibration.pixelDepth = settings.workingVoxelSize;
		calibration.setUnit( "micrometer" );
		imp.setCalibration( calibration );

		final CompositeImage compositeImage = new CompositeImage( imp );

		final int numChannels = (int) rescaledVolumes.size();

		for ( int c = 0; c < numChannels; c++ )
		{
			compositeImage.setC( c + 1 );
			if ( c == ( int ) settings.dnaChannelIndex )
				IJ.run( compositeImage, "Blue", "" );
			else if ( c == ( int ) settings.tubulinChannelIndex )
 				IJ.run( compositeImage, "Green", "" );
			else
				IJ.run( compositeImage, "Grays", "" );

			if ( compositeImage.getBitDepth() == 8 )
				compositeImage.setDisplayRange( 0, 255 );
			else if ( compositeImage.getBitDepth() == 16 )
				compositeImage.setDisplayRange( 0, 65535 );

		}


		int lastChannelIndex = numChannels;

		// Binary masks

		compositeImage.setC( ++lastChannelIndex );
		IJ.run( compositeImage, "Blue", "" );
		compositeImage.setDisplayRange( 0, 3 );

		compositeImage.setC( ++lastChannelIndex );
		IJ.run( compositeImage, "Green", "" );
		compositeImage.setDisplayRange( 0, 3 );

		// Interest points

		compositeImage.setC( ++lastChannelIndex  );
		IJ.run( compositeImage, "Grays", "" );
		compositeImage.setDisplayRange( 0, 1 );

		compositeImage.setDisplayMode( CompositeImage.COMPOSITE );

		return compositeImage;
	}

	private void drawTransformedPoint(
			AffineTransform3D alignmentTransform,
			RandomAccessibleInterval< R > interestPointsImage,
			double[] point,
			int value )
	{
		final double[] transformedPoint = transformToCoverslipCoordinateSystem( alignmentTransform, point );

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

	private double[] transformToCoverslipCoordinateSystem(
			AffineTransform3D alignmentTransform,
			double[] point )
	{
		final double[] transformedPoint = Utils.copy( point );

		alignmentTransform.inverse().apply( point, transformedPoint ); // transform back

		// Utils.divide( transformedPoint, settings.workingVoxelSize ); // Scale to pixel units

		return transformedPoint;
	}

	private void drawPoint(
			RandomAccessibleInterval< R > rai,
			double[] position,
			double calibratedRadius,
			double calibration,
			int value )
	{
		Shape shape = new HyperSphereShape( (int) Math.ceil( calibratedRadius / calibration ) );
		final RandomAccessible< Neighborhood< R > > nra = shape.neighborhoodsRandomAccessible( rai );
		final RandomAccess< Neighborhood< R > > neighborhoodRandomAccess = nra.randomAccess();

		final double[] pixelPositionDoubles = Arrays.stream( position ).map( x -> x / calibration ).toArray();
		final long[] pixelPositionLongs = Utils.asRoundedLongs( pixelPositionDoubles );
		neighborhoodRandomAccess.setPosition( pixelPositionLongs );
		final Neighborhood< R > neighborhood = neighborhoodRandomAccess.get();

		final long[] pos = new long[ 3 ];
		final Cursor< R > cursor = neighborhood.cursor();
		while( cursor.hasNext() )
		{
			try
			{
				cursor.fwd();
				cursor.localize( pos );
				cursor.get().setReal( value );
			}
			catch ( ArrayIndexOutOfBoundsException e )
			{
				Logger.log( "[ERROR] Draw points out of bounds..." );
				break;
			}
		}
	}

	private < T extends RealType< T > & NativeType< T > >
	RandomAccessibleInterval< BitType > createMask(
			RandomAccessibleInterval< T > rai,
			double threshold )
	{
		RandomAccessibleInterval< BitType > mask
				= Converters.convert( rai, ( i, o ) ->
				o.set( i.getRealDouble() > threshold ? true : false ), new BitType() );

		// "Bug" in ops requires a Views.zeroMin().
		mask = opService.morphology().fillHoles( Views.zeroMin( mask ) );

		mask = Transforms.getWithAdjustedOrigin( rai, mask );

		return mask;
	}

}
