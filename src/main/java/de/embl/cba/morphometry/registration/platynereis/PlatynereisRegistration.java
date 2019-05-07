package de.embl.cba.morphometry.registration.platynereis;

import de.embl.cba.morphometry.*;
import de.embl.cba.morphometry.geometry.CoordinatesAndValues;
import de.embl.cba.morphometry.geometry.CurveAnalysis;
import de.embl.cba.morphometry.geometry.ellipsoids.EllipsoidMLJ;
import de.embl.cba.morphometry.geometry.ellipsoids.EllipsoidVectors;
import de.embl.cba.morphometry.geometry.ellipsoids.Ellipsoids3DImageSuite;
import de.embl.cba.morphometry.geometry.ellipsoids.EllipsoidsMLJ;
import de.embl.cba.morphometry.regions.Regions;
import de.embl.cba.transforms.utils.Scalings;
import de.embl.cba.transforms.utils.Transforms;
import ij.IJ;
import ij.ImagePlus;
import net.imagej.ops.OpService;
import net.imglib2.*;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.interpolation.randomaccess.NearestNeighborInterpolatorFactory;
import net.imglib2.realtransform.AffineTransform3D;
import net.imglib2.roi.labeling.ImgLabeling;
import net.imglib2.type.NativeType;
import net.imglib2.type.logic.BitType;
import net.imglib2.type.numeric.RealType;
import net.imglib2.type.numeric.integer.IntType;
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.util.Intervals;

import java.util.List;

import static de.embl.cba.morphometry.Constants.*;
import static de.embl.cba.morphometry.viewing.BdvViewer.show;
import static de.embl.cba.transforms.utils.Scalings.createResampledArrayImg;
import static de.embl.cba.transforms.utils.Scalings.createRescaledArrayImg;
import static de.embl.cba.transforms.utils.Transforms.getScalingFactors;
import static java.lang.Math.toRadians;


public class PlatynereisRegistration< T extends RealType< T > & NativeType< T > >
{
	final PlatynereisRegistrationSettings settings;
	final OpService opService;

	private AffineTransform3D transformAtRegistrationResolution;
	private Img< IntType > watershedLabelImg;
	private AffineTransform3D registration;
	private double[] registrationCalibration;
	private RandomAccessibleInterval< BitType > yawAndOrientationAlignedIntensity;
	private RandomAccessibleInterval< T > isotropic;
	private RandomAccessibleInterval< T > image;
	private RandomAccessibleInterval< BitType > mask;
	private RandomAccessibleInterval< BitType > yawAlignedMask;
	private RandomAccessibleInterval yawAlignedIntensity;
	private CoordinateAndValue axialEmbryoCenter;
	private double[] inputCalibration;

	public PlatynereisRegistration(
			final PlatynereisRegistrationSettings settings,
			final OpService opService )
	{
		this.settings = settings;
		this.opService = opService;
	}

	public boolean run( RandomAccessibleInterval< T > image, double[] inputCalibration )
	{
		this.inputCalibration = inputCalibration;
		this.image = image;

//		if ( settings.showIntermediateResults )
//			show( image, "input image", null, inputCalibration, false );

		registration = new AffineTransform3D();

		isotropic = downSampleToRegistrationResolution( image );

//		isotropic = downSampleUsingImageJ1( image );

		if ( settings.invertImage ) invertImage();

		if ( ! segmentPlaty() ) return false;

		final AffineTransform3D ellipsoidFittingBasedTransformation = computeEllipsoidFittingBasedTransformation( this.mask );

		applyAlignmentToImageAndMask( ellipsoidFittingBasedTransformation );

		orientLongAxis();

		rollTransform();

		transformAtRegistrationResolution = registration;

		return true;

	}

	public RandomAccessibleInterval< T >  downSampleUsingImageJ1( RandomAccessibleInterval< T > image )
	{
		Logger.log( "Down-sampling to registration resolution: " + settings.registrationResolution + " micrometer" );
		final double[] scalingFactors = getScalingFactors( inputCalibration, settings.registrationResolution );
		Logger.log( "Scaling factor X: " + scalingFactors[ 0 ] );
		Logger.log( "Scaling factor Y: " + scalingFactors[ 1 ] );
		Logger.log( "Scaling factor Z: " + scalingFactors[ 2 ] );

		final ImagePlus imp = ImageJFunctions.wrap( image, "" );

		IJ.run( imp, "Scale...",
				"x=" + scalingFactors[ 0 ] +
						"y=" + scalingFactors[ 1 ] +
						"z=" + scalingFactors[ 2 ] + " interpolation=Bilinear average process create" );

		imp.close();

		final ImagePlus downScaled = IJ.getImage();

		return ImageJFunctions.wrapReal( downScaled );
	}

	public void invertImage()
	{
		isotropic = Utils.copyAsArrayImg( Utils.invertedView( isotropic ) );
		if ( settings.showIntermediateResults )
			show( isotropic,
					"inverted image",
					null, registrationCalibration, false );
	}

	private void applyAlignmentToImageAndMask( AffineTransform3D alignmentTransform )
	{
		Logger.log( "Transform image using Ellipsoid parameters..." );

		registration.preConcatenate( alignmentTransform );

		final RandomAccessibleInterval transformedView =
				Transforms.createTransformedView(
						mask, registration, new NearestNeighborInterpolatorFactory() );

		yawAlignedMask = Utils.copyAsArrayImg( transformedView );

		yawAlignedIntensity = Utils.copyAsArrayImg(
				Transforms.createTransformedView( isotropic, registration ) );

		if ( settings.showIntermediateResults )
			show( yawAlignedMask, "aligned mask",
					Transforms.origin(), registrationCalibration, false );

		if ( settings.showIntermediateResults )
			show( yawAlignedIntensity, "aligned intensities",
					Transforms.origin(), registrationCalibration, false );
	}

	private AffineTransform3D computeEllipsoidFittingBasedTransformation( RandomAccessibleInterval< BitType > mask )
	{
		/**
		 * Compute ellipsoid (probably mainly yaw) alignment
		 * - https://en.wikipedia.org/wiki/Euler_angles
		 */

		Logger.log( "Fit ellipsoid..." );

		final ImagePlus imagePlus = ImageJFunctions.wrap( Utils.get3DRaiAs5DRaiWithImagePlusDimensionOrder( mask  ), "mask" );
		final EllipsoidVectors ellipsoidVectors = Ellipsoids3DImageSuite.fitEllipsoid( imagePlus );
		final AffineTransform3D imageSuiteTransform = Ellipsoids3DImageSuite.createAlignmentTransform( ellipsoidVectors );

		final EllipsoidMLJ ellipsoidMLJ = EllipsoidsMLJ.computeParametersFromBinaryImage( mask );
		final AffineTransform3D mljTransform = EllipsoidsMLJ.createAlignmentTransform( ellipsoidMLJ );

		return imageSuiteTransform;
	}

	private void rollTransform()
	{

		final AffineTransform3D rollTransform =
				computeIntensityBasedRollTransform( yawAndOrientationAlignedIntensity );

		// changes whether the found structure should be at the top or bottom
		rollTransform.rotate( X, Math.PI );

		registration = registration.preConcatenate( rollTransform  );

		if ( settings.showIntermediateResults )
			show( Transforms.createTransformedView( isotropic, registration ),
					"aligned at registration resolution",
					Transforms.origin(), registrationCalibration, false );
	}

	private void orientLongAxis()
	{
		Logger.log( "Computing long axis orientation..." );

		final AffineTransform3D flippingTransform =
				computeFlippingTransform(
						yawAlignedIntensity,
						yawAlignedMask,
						settings.registrationResolution );

		registration = registration.preConcatenate( flippingTransform );

		yawAndOrientationAlignedIntensity = Utils.copyAsArrayImg(
				Transforms.createTransformedView(
						yawAlignedIntensity,
						flippingTransform,
						new NearestNeighborInterpolatorFactory() ) );

		if ( settings.showIntermediateResults )
			show( yawAndOrientationAlignedIntensity, "long axis aligned and oriented",
					Transforms.origin(), registrationCalibration, false );


	}

	private boolean segmentPlaty()
	{
		createMask();

		return true;
	}

	private RandomAccessibleInterval< DoubleType > distanceTransform()
	{
		/**
		 * Distance transform
		 * - Note: EUCLIDIAN distances are returned as squared distances
		 */

		Logger.log( "Distance transform..." );

		final RandomAccessibleInterval< DoubleType > distances = Algorithms.computeSquaredDistances( mask );

		if ( settings.showIntermediateResults )
			show( distances, "squared distances", null,
					registrationCalibration, false );
		return distances;
	}

	private void createMask()
	{
		mask = binarise();

		open( mask );

		removeSmallRegions( mask );

//		fillHolesInMask();

	}

	private void open( RandomAccessibleInterval< BitType > mask )
	{
		this.mask = Algorithms.erode( mask, 1 );
		this.mask = Algorithms.dilate( this.mask, 1 );
	}

	private RandomAccessibleInterval< BitType > binarise()
	{
		/**
		 *  Compute threshold
		 */

		double threshold = Algorithms.thresholdYen( isotropic );
		Logger.log( "Threshold: " + threshold );

		/**
		 * Create mask
		 */

		mask = Algorithms.createMask( isotropic, threshold );

		if ( settings.showIntermediateResults )
			show( Utils.copyAsArrayImg( mask ), "binary mask", null,
					registrationCalibration, false );

		return mask;
	}

	private void removeSmallRegions( RandomAccessibleInterval< BitType > mask )
	{
		Regions.removeSmallRegionsInMask(
				mask,
				settings.minimalObjectSize,
				settings.registrationResolution );

		if ( settings.showIntermediateResults )
			show( Utils.copyAsArrayImg( this.mask ), "small regions removed", null,
					registrationCalibration, false );
	}

	private void fillHolesInMask()
	{
		Logger.log( "Fill holes..." );
		for ( int d = 0; d < 3; ++d )
			mask = Algorithms.fillHoles3Din2D( mask, d, opService );

		if ( settings.showIntermediateResults )
			show( mask, "small regions removed and holes closed", null,
					registrationCalibration, false );
	}

	private RandomAccessibleInterval< T > downSampleToRegistrationResolution( RandomAccessibleInterval< T > image )
	{
		/**
		 *  Down-sampling to registration resolution
		 */

		Logger.log( "Down-sampling to registration resolution: " + settings.registrationResolution + " micrometer" );

		final double[] scalingFactors = getScalingFactors( inputCalibration, settings.registrationResolution );

		Logger.log( "Scaling factor X: " + scalingFactors[ 0 ] );
		Logger.log( "Scaling factor Y: " + scalingFactors[ 1 ] );
		Logger.log( "Scaling factor Z: " + scalingFactors[ 2 ] );

		isotropic = Scalings.createRescaledCellImg( image, scalingFactors );

		registrationCalibration = Utils.as3dDoubleArray( settings.registrationResolution );

		if ( settings.showIntermediateResults )
			show( isotropic,
					"isotropic sampled at registration resolution",
					null, registrationCalibration, false );

		return isotropic;
	}


	private ImgLabeling< Integer, IntType > computeWatershed(
			RandomAccessibleInterval< BitType > mask,
			RandomAccessibleInterval< DoubleType > distances,
			ImgLabeling< Integer, IntType > seedsLabelImg )
	{

		Logger.log( "Watershed..." );

		watershedLabelImg = ArrayImgs.ints( Intervals.dimensionsAsLongArray( mask ) );
		final ImgLabeling< Integer, IntType > watershedLabeling =
				new ImgLabeling<>( watershedLabelImg );

		if ( settings.showIntermediateResults )
			show( watershedLabelImg, "watershed",
					null, registrationCalibration, false );

		opService.image().watershed(
				watershedLabeling,
				Utils.invertedView( distances ),
				seedsLabelImg,
				false,
				false );

		Utils.applyMask( watershedLabelImg, mask );

		return watershedLabeling;
	}


	public AffineTransform3D getRegistrationTransform( double[] inputCalibration, double outputResolution )
	{
		final AffineTransform3D transform =
				Transforms.getScalingTransform( inputCalibration, settings.registrationResolution )
						.preConcatenate( transformAtRegistrationResolution.copy() )
						.preConcatenate( Transforms.getScalingTransform( settings.registrationResolution, outputResolution ) );

		return transform;
	}

	public RandomAccessibleInterval< BitType >
	getAlignedMask( double resolution, FinalInterval interval )
	{

		/**
		 * - TODO: using the mask just like this was cutting away signal from embryo..
		 * 		   the issue might be that during the rotations the voxels do not end up
		 * 		   precisely where they should be? Currently, I simple dilate "a bit".
		 * 		   Feels kind of messy...better way?
		 */

		Logger.log( "Creating aligned mask..." );

		final RandomAccessibleInterval< BitType > dilatedMask = Algorithms.dilate( mask, 2 );

		AffineTransform3D transform = transformAtRegistrationResolution.copy()
				.preConcatenate( Transforms.getScalingTransform( settings.registrationResolution, resolution ) );

		RandomAccessibleInterval< BitType > alignedMask =
				Utils.copyAsArrayImg(
					Transforms.createTransformedView(
							dilatedMask,
							transform,
							interval, // after the transform we need to specify where we want to "crop"
							new NearestNeighborInterpolatorFactory() // binary image => do not interpolate linearly!
					)
				);

		if ( settings.showIntermediateResults )
			show( alignedMask, "aligned mask at output resolution",
					Transforms.origin(), resolution );

		return alignedMask;
	}

	private  < T extends RealType< T > & NativeType< T > >
	AffineTransform3D computeIntensityBasedRollTransform(
			RandomAccessibleInterval< T > image )
	{
		Logger.log( "Computing intensity based roll transform" );

		final AffineTransform3D intensityBasedRollTransform =
				computeIntensityBasedRollTransform(
					image,
					settings.projectionXMin,
					settings.projectionXMax,
					settings.projectionBlurSigma,
					registrationCalibration );

		return intensityBasedRollTransform;
	}


	private AffineTransform3D computeFlippingTransform(
			RandomAccessibleInterval yawAlignedIntensities,
			RandomAccessibleInterval< BitType > yawAlignedMask,
			double calibration )
	{
//
//		final CoordinatesAndValues coordinatesAndValues =
//				Utils.computeAverageIntensitiesAlongAxisWithinMask(
//						yawAlignedIntensities,
//						yawAlignedMask,
//						X,
//						calibration );

		final CoordinatesAndValues coordinatesAndValues =
				Utils.computeAverageIntensitiesAlongAxis(
						yawAlignedIntensities,
						100.0,
						X,
						calibration );

		if ( settings.showIntermediateResults )
			Plots.plot(
					coordinatesAndValues.coordinates,
					coordinatesAndValues.values,
					"x",
					"average intensity" );

		CoordinateAndValue maximum =
				CurveAnalysis.maximum( coordinatesAndValues, null );

		AffineTransform3D affineTransform3D = new AffineTransform3D();

		if ( maximum.coordinate > 0 ) affineTransform3D.rotate( Z, toRadians( 180.0D ) );

		return affineTransform3D;
	}

	private < T extends RealType< T > & NativeType< T > >
	AffineTransform3D computeIntensityBasedRollTransform(
			RandomAccessibleInterval rai,
			double xMin,
			double xMax,
			double blurSigma,
			double[] registrationCalibration )
	{
		final RandomAccessibleInterval< T > longAxisProjection =
				Utils.createAverageProjectionAlongAxis(
					rai,
					X,
					xMin,
					xMax,
					settings.registrationResolution );

		if ( settings.showIntermediateResults )
			show( longAxisProjection, "long axis projection", null,
					registrationCalibration, false );

		final RandomAccessibleInterval< T > blurred = Utils.createBlurredRai(
				longAxisProjection,
				blurSigma,
				settings.registrationResolution );

		final Point maximum = Algorithms.getMaximumLocation( blurred, Utils.as2dDoubleArray( settings.registrationResolution ));
		final List< RealPoint > realPoints = Utils.asRealPointList( maximum );
		realPoints.add( new RealPoint( new double[]{ 0, 0 } ) );

		if ( settings.showIntermediateResults )
			show( blurred, "perpendicular projection - blurred ",
					realPoints, Utils.as2dDoubleArray( settings.registrationResolution ), false );

		// take the position of the maximum and create a roll transform such that this
		// maximum is oriented in the yz plane
//		final AffineTransform3D xAxisRollTransform = createXAxisRollTransform( maximum );

		// just flip depending on where the maximum is along the y-axis
		AffineTransform3D xAxisRollTransform = new AffineTransform3D();
		if ( maximum.getDoublePosition( Y ) > 0 ) xAxisRollTransform.rotate( X, toRadians( 180.0D ) );

		return xAxisRollTransform;
	}


	private static
	AffineTransform3D createXAxisRollTransform( Point maximum2DinYZPlane )
	{
		double angleToZAxisInDegrees = Angles.angle2DToCoordinateSystemsAxisInDegrees( maximum2DinYZPlane );
		AffineTransform3D rollTransform = new AffineTransform3D();

		Logger.log( "Roll angle: " + angleToZAxisInDegrees );
		rollTransform.rotate( X, toRadians( angleToZAxisInDegrees ) );

		return rollTransform;
	}


}
