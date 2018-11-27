package de.embl.cba.morphometry.spindle;

import de.embl.cba.morphometry.*;
import de.embl.cba.morphometry.geometry.CoordinatesAndValues;
import de.embl.cba.morphometry.geometry.EllipsoidParameters;
import de.embl.cba.morphometry.geometry.Ellipsoids;
import ij.IJ;
import ij.ImagePlus;
import net.imagej.ops.OpService;
import net.imglib2.*;
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

import java.io.File;
import java.util.ArrayList;

import static de.embl.cba.morphometry.Algorithms.angleInRadians;
import static de.embl.cba.morphometry.Transforms.getScalingFactors;
import static de.embl.cba.morphometry.viewing.BdvViewer.show;


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

		final double[] workingCalibration = Utils.as3dDoubleArray( settings.workingVoxelSize );

		final RandomAccessibleInterval< T > dapi = Algorithms.createRescaledArrayImg( settings.dapiImage, getScalingFactors( settings.inputCalibration, settings.workingVoxelSize ) );
		final RandomAccessibleInterval< T > tubulin = Algorithms.createRescaledArrayImg( settings.tubulinImage, getScalingFactors( settings.inputCalibration, settings.workingVoxelSize ) );

		if ( settings.showIntermediateResults ) show( dapi, "image isotropic resolution", null, workingCalibration, false );
		if ( settings.showIntermediateResults ) show( tubulin, "tubulin isotropic resolution", null, workingCalibration, false );


		/**
		 *  Compute offset and threshold
		 */

		final RandomAccessibleInterval< T > dapiDownscaledToSpindleWidth = Algorithms.createRescaledArrayImg( settings.dapiImage, getScalingFactors( settings.inputCalibration, 3.0 ) );
		final double maximumValue = Algorithms.getMaximumValue( dapiDownscaledToSpindleWidth );
		double threshold = maximumValue / 2.0;

		Utils.log( "Dapi threshold: " + threshold );


		/**
		 * Create mask
		 */

		RandomAccessibleInterval< BitType > mask = createMask( dapi, threshold );

		if ( settings.showIntermediateResults ) show( mask, "dapi mask", null, workingCalibration, false );

		/**
		 * Morphological filtering
		 *
		 * - it appears that the principal axes determination is not working robustly
		 * if the meta-phase plate is too thick
		 */

		Utils.log( "Perform morphological filtering on dapi mask..." );

		RandomAccessibleInterval< BitType > closed = Algorithms.erode( mask, (int) Math.ceil( 1.0 / settings.workingVoxelSize )  );

		if ( settings.showIntermediateResults ) show( closed, "filtered dapi mask", null, workingCalibration, false );


		/**
		 * Extract metaphase plate object
		 */

		Utils.log( "Extracting metaphase plate object..." );

		final ImgLabeling< Integer, IntType > labelImg = Utils.asImgLabeling( closed );

		final LabelRegion< Integer > largestObject = Algorithms.getLargestObject( labelImg );

		final Img< BitType > dapiMask = Algorithms.createMaskFromLabelRegion( largestObject, Intervals.dimensionsAsLongArray( labelImg ) );

		if ( settings.showIntermediateResults ) show( dapiMask, "meta-phase object", null, workingCalibration, false );

		/**
		 * Compute ellipsoid alignment
		 */

		Utils.log( "Determining meta-phase plate axes..." );

		final EllipsoidParameters ellipsoidParameters = Ellipsoids.computeParametersFromBinaryImage( dapiMask );

		Utils.log( "Creating aligned images..." );

		final AffineTransform3D alignmentTransform = Ellipsoids.createAlignmentTransform( ellipsoidParameters );
		final RandomAccessibleInterval aligendTubulin = Utils.copyAsArrayImg( Transforms.createTransformedView( tubulin, alignmentTransform, new NearestNeighborInterpolatorFactory() ) );
		final RandomAccessibleInterval alignedDapi = Utils.copyAsArrayImg( Transforms.createTransformedView( dapi, alignmentTransform ) );
		final RandomAccessibleInterval alignedDapiMask = Utils.copyAsArrayImg( Transforms.createTransformedView( dapiMask, alignmentTransform ) );

		if ( settings.showIntermediateResults ) show( alignedDapi, "aligned dapi", Transforms.origin(), workingCalibration, false );
		if ( settings.showIntermediateResults ) show( alignedDapiMask, "aligned dapi mask", Transforms.origin(), workingCalibration, false );

		RandomAccessibleInterval< T > interestPoints = Utils.copyAsArrayImg( dapi );
		Utils.setValues( interestPoints, 0.0 );

		final double[] dnaCenter = new double[ 3 ];
		alignmentTransform.inverse().apply( new double[]{0,0,0}, dnaCenter );
		drawPoint( interestPoints, dnaCenter, settings.interestPointsRadius, settings.workingVoxelSize );


		/**
		 * Compute metaphase plate width
		 */

		Utils.log( "Determining widths..." );

		final CoordinatesAndValues dapiProfile = Utils.computeAverageIntensitiesAlongAxis( alignedDapi, settings.maxShortAxisDist, 2, settings.workingVoxelSize );
		if ( settings.showIntermediateResults ) Plots.plot( dapiProfile.coordinates, dapiProfile.values, "distance to center", "dapi intensity along shortest axis" );

		// TODO: compute meta-phase width, using FWHM and write to log window

		/**
		 * Compute spindle length
		 */

		final CoordinatesAndValues tubulinProfile = Utils.computeMaximumIntensitiesAlongAxis( aligendTubulin, settings.maxShortAxisDist, 2, settings.workingVoxelSize );
		if ( settings.showIntermediateResults ) Plots.plot( tubulinProfile.coordinates, tubulinProfile.values, "distance to center", "tubulin maximal intensities" );

		final CoordinatesAndValues tubulinProfileDerivative = Algorithms.computeDerivatives( tubulinProfile, (int) Math.ceil( settings.derivativeDelta / settings.workingVoxelSize ) );
		if ( settings.showIntermediateResults ) Plots.plot( tubulinProfileDerivative.coordinates, tubulinProfileDerivative.values, "distance to center", "tubulin intensity derivative" );

		double[] spindlePoles = getLeftMaxAndRightMinLoc( tubulinProfileDerivative.coordinates, tubulinProfileDerivative.values );

		Utils.log( "Left spindle pole found at: " + spindlePoles[ 0 ] );
		Utils.log( "Right spindle pole found at: " + spindlePoles[ 1 ] );

		addSpindlePolesAsInterestPoints( alignmentTransform, interestPoints, spindlePoles );

		final ArrayList< RealPoint > spindleLengthPoints = new ArrayList<>(  );

		spindleLengthPoints.add( new RealPoint( new double[]{ 0.0, 0.0, 0 } ));
		spindleLengthPoints.add( new RealPoint( new double[]{ 0.0, 0.0, spindlePoles[ 0 ] } ));
		spindleLengthPoints.add( new RealPoint( new double[]{ 0.0, 0.0, spindlePoles[ 1 ] } ));

		if ( settings.showIntermediateResults ) show( aligendTubulin, "aligned tubulin", spindleLengthPoints, workingCalibration, false );


		/**
		 * Create output images
		 */

		AffineTransform3D rotation = alignmentTransform.copy();
		rotation.setTranslation( new double[ 3 ] );

		final double[] zAxis = { 0, 0, 1 };
		final double[] xAxis = { 1, 0, 0 };

		// OCS = Orthogonal Coverslip System
		// Create views that do not tilt the coverslip
		// This was biologically interesting to see how the axes
		// are with respect to the coverslip surface

		final double[] shortDnaAxisOCS = new double[ 3 ];

		rotation.inverse().apply( zAxis, shortDnaAxisOCS );

		final double[] sideViewDirectionOCS = new double[ 3 ];
		LinAlgHelpers.cross( shortDnaAxisOCS, zAxis, sideViewDirectionOCS );

		final double dot = LinAlgHelpers.dot( xAxis, sideViewDirectionOCS );
		final double angleInRadians = angleInRadians( xAxis, sideViewDirectionOCS );

		rotation = new AffineTransform3D();

		// TODO: maybe the sign of the angle depends on its magnitude??
		rotation.rotate( 2,  +angleInRadians );

		final RandomAccessibleInterval transformedDapiView = Transforms.createTransformedView( dapi, rotation );
		final RandomAccessibleInterval transformedTubulinView = Transforms.createTransformedView( tubulin, rotation );
		final RandomAccessibleInterval transformedInterestPointView = Transforms.createTransformedView( interestPoints, rotation );

//		Bdv bdv = null;
//		if ( settings.showIntermediateResults ) BdvFunctions.show( interestPoints, "" ).getBdvHandle();
//		if ( settings.showIntermediateResults ) bdv = BdvFunctions.show( transformedDapiView, "" ).getBdvHandle();
//		if ( settings.showIntermediateResults ) BdvFunctions.show( transformedInterestPointView, "", BdvOptions.options().addTo( bdv ) );

		Utils.log( "Saving result images..." );

		saveImagePlus( Utils.asImagePlus( alignedDapiMask ) );
		saveMaximumProjections( transformedDapiView, "dapi" );
		saveMaximumProjections( transformedTubulinView, "tubulin" );
		saveMaximumProjections( transformedInterestPointView, "points" );

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

	public void saveImagePlus( ImagePlus imagePlus )
	{

		final String path = settings.outputDirectory.getAbsolutePath()
					+ File.separator
					+ settings.inputDataSetName
					+ File.separator
					+ imagePlus.getTitle() + ".tiff";

		new File(path).getParentFile().mkdirs();

		IJ.saveAsTiff( imagePlus, path );

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
		Utils.log( "Creating mask...");

		RandomAccessibleInterval< BitType > mask = Converters.convert( downscaled, ( i, o ) -> o.set( i.getRealDouble() > threshold ? true : false ), new BitType() );

		mask = opService.morphology().fillHoles( mask );

		return mask;
	}



}
