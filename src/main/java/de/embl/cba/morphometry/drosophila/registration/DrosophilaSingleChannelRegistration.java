package de.embl.cba.morphometry.drosophila.registration;

import de.embl.cba.morphometry.*;
import de.embl.cba.morphometry.geometry.CentroidsParameters;
import de.embl.cba.morphometry.geometry.CoordinatesAndValues;
import de.embl.cba.morphometry.geometry.CurveAnalysis;
import de.embl.cba.morphometry.geometry.ellipsoids.EllipsoidMLJ;
import de.embl.cba.morphometry.geometry.ellipsoids.EllipsoidsMLJ;
import de.embl.cba.morphometry.refractiveindexmismatch.RefractiveIndexMismatchCorrectionSettings;
import de.embl.cba.morphometry.refractiveindexmismatch.RefractiveIndexMismatchCorrections;
import de.embl.cba.morphometry.regions.Regions;
import de.embl.cba.transforms.utils.Transforms;
import net.imagej.ops.OpService;
import net.imglib2.*;
import net.imglib2.algorithm.labeling.ConnectedComponents;
import net.imglib2.algorithm.neighborhood.HyperSphereShape;
import net.imglib2.histogram.Histogram1d;
import net.imglib2.img.Img;
import net.imglib2.img.array.ArrayImgs;
import net.imglib2.interpolation.randomaccess.NearestNeighborInterpolatorFactory;
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
import net.imglib2.type.numeric.real.DoubleType;
import net.imglib2.util.Intervals;
import net.imglib2.view.Views;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Set;

import static de.embl.cba.morphometry.Constants.*;
import static de.embl.cba.morphometry.viewing.BdvViewer.show;
import static de.embl.cba.transforms.utils.Scalings.createRescaledArrayImg;
import static de.embl.cba.transforms.utils.Transforms.getScalingFactors;
import static java.lang.Math.toRadians;


public class DrosophilaSingleChannelRegistration< T extends RealType< T > & NativeType< T > >
{
	final DrosophilaRegistrationSettings settings;
	final OpService opService;

	private RandomAccessibleInterval< BitType > embryoMask;
	private double coverslipPosition;
	private AffineTransform3D transformAtRegistrationResolution;
	private Img< IntType > watershedLabelImg;
	private double[] correctedCalibration;
	private AffineTransform3D registration;
	private double[] registrationCalibration;
	private RandomAccessibleInterval< BitType > yawAndOrientationAlignedIntensity;
	private RandomAccessibleInterval< T > isotropic;
	private RandomAccessibleInterval< T > intensityCorrected;
	private RandomAccessibleInterval< BitType > mask;
	private RandomAccessibleInterval< BitType > yawAlignedMask;
	private RandomAccessibleInterval yawAlignedIntensities;
	private double axialEmbryoCenterPosition;
	private EllipsoidMLJ ellipsoidParameters;

	public DrosophilaSingleChannelRegistration(
			final DrosophilaRegistrationSettings settings,
			final OpService opService )
	{
		this.settings = settings;
		this.opService = opService;
	}

	public void run( RandomAccessibleInterval< T > image, double[] inputCalibration )
	{
		if ( settings.showIntermediateResults )
			show( image, "input image", null, inputCalibration, false );

		registration = new AffineTransform3D();

		refractiveIndexScalingCorrection( image, inputCalibration );

		createIsotropicImage( image );

		refractiveIndexIntensityCorrection();

		if ( ! segmentEmbryo() ) return;

		computeEllipsoidParameters();

		if ( settings.onlyComputeEllipsoidParameters ) return;

		applyYawAlignmentToImageAndMask();

		orientLongAxis();

		rollTransform();

		transformAtRegistrationResolution = registration;

	}

	public EllipsoidMLJ getEllipsoidParameters()
	{
		return ellipsoidParameters;
	}

	private void applyYawAlignmentToImageAndMask()
	{
		registration.preConcatenate(
				EllipsoidsMLJ.createAlignmentTransform( ellipsoidParameters ) );

		yawAlignedMask = Utils.copyAsArrayImg( Transforms.createTransformedView( embryoMask, registration, new NearestNeighborInterpolatorFactory() ) );

		yawAlignedIntensities = Utils.copyAsArrayImg( Transforms.createTransformedView( isotropic, registration ) );
	}

	private void computeEllipsoidParameters()
	{
		/**
		 * Compute ellipsoid (probably mainly yaw) alignment
		 * - https://en.wikipedia.org/wiki/Euler_angles
		 */

		Logger.log( "Fit ellipsoid..." );

		ellipsoidParameters = EllipsoidsMLJ.computeParametersFromBinaryImage( embryoMask );
	}

	private void rollTransform()
	{
		/**
		 *  Roll transform
		 */

		final AffineTransform3D rollTransform = computeIntensityBasedRollTransform( yawAndOrientationAlignedIntensity );

		rollTransform.rotate( X, Math.PI ); // this changes whether the found structure should be at the top or bottom

		registration = registration.preConcatenate( rollTransform  );

		if ( settings.showIntermediateResults )
			show( Transforms.createTransformedView( intensityCorrected, registration ), "aligned at registration resolution", Transforms.origin(), registrationCalibration, false );
	}

	private void orientLongAxis()
	{
		/**
		 *  Long axis orientation
		 */

		Logger.log( "Computing long axis orientation..." );

		final AffineTransform3D orientationTransform = computeFlippingTransform( yawAlignedMask, yawAlignedIntensities, settings.registrationResolution );

		registration = registration.preConcatenate( orientationTransform );

		yawAndOrientationAlignedIntensity = Utils.copyAsArrayImg( Transforms.createTransformedView( yawAlignedIntensities, registration, new NearestNeighborInterpolatorFactory() ) );

		if ( settings.showIntermediateResults ) show( yawAndOrientationAlignedIntensity, "long axis aligned and oriented", Transforms.origin(), registrationCalibration, false );


	}

	private boolean segmentEmbryo()
	{

		createMask();

		/**
		 * Distance transform
		 * - Note: EUCLIDIAN distances are returned as squared distances
		 */

		Logger.log( "Distance transform..." );

		final RandomAccessibleInterval< DoubleType > distances = Algorithms.computeSquaredDistances( mask );

		if ( settings.showIntermediateResults )
			show( distances, "squared distances", null, registrationCalibration, false );


		/**
		 * Watershed seeds
		 * - if local maxima are defined as strictly larger (>) one misses them in case two
		 *   neighboring pixels in the centre of an object have the same distance value
		 * - if local maxima are >= and the search radius is only 1 pixel (four-connected) one gets false maxima at the corners objects
		 *   thus, the search radius should be always >= 2 pixels
		 * - triangle shaped appendices are an issue because they do not have a
		 *   maximum in the distance map
		 * - due to the elongated shape of the embryos there might not be a clear maximum => use also a global threshold
		 */

		final ImgLabeling< Integer, IntType > seedsLabelImg = createWatershedSeeds( distances );

		/**
		 * Watershed
		 */

		final ImgLabeling< Integer, IntType > labeling = computeWatershed( mask, distances, seedsLabelImg );

		if ( settings.showIntermediateResults )
			show( watershedLabelImg, "watershed", null, registrationCalibration, false );

		/**
		 * Get main embryo
		 */

		Logger.log( "Extract central embryo..." );

		final double[] center = getApproximateEmbryoCenter( labeling );

		final Set< LabelRegion< Integer > > centralRegions =
				Regions.getCentralRegions( labeling, center, (int) ( settings.centralRegionDistance / settings.registrationResolution ) );

		if ( centralRegions.size() == 0 ) return false;

		embryoMask = Algorithms.createMaskFromLabelRegions(
				centralRegions,
				Intervals.dimensionsAsLongArray( labeling ) );

		if ( settings.showIntermediateResults )
			show( Utils.copyAsArrayImg( embryoMask ), "embryo mask", null, registrationCalibration, false );

		/**
		 * Process main embryo mask
		 * - TODO: put the currently hard-coded values into settings
		 */

		embryoMask = Algorithms.close( embryoMask, ( int ) ( 20.0 / settings.registrationResolution ) );

		// embryoMask = Algorithms.open( embryoMask, ( int ) ( 20.0 / settings.registrationResolution ) );

		if ( settings.showIntermediateResults )
			show( embryoMask, "embryo mask - processed", null, registrationCalibration, false );

		return true;
	}

	private double[] getApproximateEmbryoCenter( ImgLabeling< Integer, IntType > labeling )
	{
		return new double[]{
					labeling.dimension( 0 ) / 2.0 ,
					labeling.dimension( 1 ) / 2.0,
					axialEmbryoCenterPosition / settings.registrationResolution };
	}

	private void createMask()
	{
		/**
		 *  Compute threshold
		 *  - TODO: How does Huang work?
		 *  - TODO: find some more scientific method to determine threshold...
		 */

		final Histogram1d< T > histogram = opService.image().histogram( Views.iterable( intensityCorrected ) );
		final double huang = opService.threshold().huang( histogram ).getRealDouble();
//		final double otsu = opService.threshold().otsu( histogram ).getRealDouble();
//		final double yen = opService.threshold().yen( histogram ).getRealDouble();

		double thresholdAfterIntensityCorrection = huang;

		Logger.log( "Threshold (after intensity correction): " + thresholdAfterIntensityCorrection );


		/**
		 * Create mask
		 */

		mask = Algorithms.createMask( intensityCorrected, thresholdAfterIntensityCorrection );

		if ( settings.showIntermediateResults ) show( mask, "binary mask", null, registrationCalibration, false );


		/**
		 * Process mask
		 * - remove small objects
		 * - close holes
		 */

		Regions.removeSmallRegionsInMask( mask, settings.minimalObjectSize, settings.registrationResolution );

		for ( int d = 0; d < 3; ++d )
		{
			mask = Algorithms.fillHoles3Din2D( mask, d, opService );
		}

		if ( settings.showIntermediateResults ) show( mask, "small objects removed and holes closed", null, registrationCalibration, false );

	}

	private void refractiveIndexIntensityCorrection()
	{
		/**
		 *  Compute intensity offset (for refractive index mismatch corrections)
		 */

		Logger.log( "Offset and threshold..." );

		final IntensityHistogram histogram = new IntensityHistogram( isotropic, 65535.0, 5.0 );

		CoordinateAndValue histogramMode = histogram.getMode();

		Logger.log( "Intensity offset: " + histogramMode.coordinate );


		/**
		 *  Compute approximate axial embryo center and coverslip coordinate
		 */

		final CoordinatesAndValues averageIntensitiesAlongZ =
				Utils.computeAverageIntensitiesAlongAxis( isotropic, 2, settings.registrationResolution );

		if ( settings.showIntermediateResults )
			Plots.plot( averageIntensitiesAlongZ.coordinates, averageIntensitiesAlongZ.values, "z [um]", "average intensities" );

		axialEmbryoCenterPosition = CurveAnalysis.maxLoc( averageIntensitiesAlongZ );
		coverslipPosition = axialEmbryoCenterPosition - DrosophilaRegistrationSettings.drosophilaWidth / 2.0;

		Logger.log( "Approximate coverslip coordinate [um]: " + coverslipPosition );
		Logger.log( "Approximate axial embryo center coordinate [um]: " + axialEmbryoCenterPosition );

		/**
		 *  Refractive index corrections
		 */

		Logger.log( "Refractive index intensity correction..." );

		final RefractiveIndexMismatchCorrectionSettings correctionSettings = new RefractiveIndexMismatchCorrectionSettings();
		correctionSettings.intensityOffset = histogramMode.coordinate;
		correctionSettings.intensityDecayLengthMicrometer = settings.refractiveIndexIntensityCorrectionDecayLength;
		correctionSettings.coverslipPositionMicrometer = coverslipPosition;
		correctionSettings.pixelCalibrationMicrometer = settings.registrationResolution;

		intensityCorrected = Utils.copyAsArrayImg( isotropic );
		RefractiveIndexMismatchCorrections.correctIntensity( intensityCorrected, correctionSettings );

		if ( settings.showIntermediateResults ) show( intensityCorrected, "intensity corrected channel 1", null, registrationCalibration, false );
	}

	private void createIsotropicImage( RandomAccessibleInterval< T > image )
	{
		/**
		 *  Down-sampling to registration resolution
		 *  - Speeds up calculations ( pow(3) effect in 3D! )
		 *  - Reduces noise
		 *  - Fills "holes" in staining
		 *  - TODO: bug: during down-sampling saturated pixels become zero
		 */

		Logger.log( "Down-sampling to registration resolution..." );

		isotropic = createRescaledArrayImg( image, getScalingFactors( correctedCalibration, settings.registrationResolution ) );
		registrationCalibration = Utils.as3dDoubleArray( settings.registrationResolution );

		if ( settings.showIntermediateResults )
			show( isotropic, "isotropic sampled at registration resolution", null, registrationCalibration, false );
	}


	private < T extends RealType< T > & NativeType< T > >
	void refractiveIndexScalingCorrection(
			RandomAccessibleInterval< T > image, double[] inputCalibration )
	{
		/**
		 *  Axial calibration correction due to refractive index mismatch
		 *  - We are using an NA 0.8 air lens imaging into water/embryo (1.33 < NA < 1.51)
		 *  - Complicated topic: http://www.bio.brandeis.edu/marderlab/axial%20shift.pdf
		 *  - We assume axial compression by factor ~1.6
		 */

		Logger.log( "Refractive index scaling correction..." );

		correctedCalibration = RefractiveIndexMismatchCorrections.getAxiallyCorrectedCalibration(
				inputCalibration, settings.refractiveIndexAxialCalibrationCorrectionFactor );

		if ( settings.showIntermediateResults )
			show( image, "corrected calibration", null, correctedCalibration, false );
	}

	private ImgLabeling< Integer, IntType > computeWatershed( RandomAccessibleInterval< BitType > mask, RandomAccessibleInterval< DoubleType > distances, ImgLabeling< Integer, IntType > seedsLabelImg )
	{
		Logger.log( "Watershed..." );

		// prepare result label image
		watershedLabelImg = ArrayImgs.ints( Intervals.dimensionsAsLongArray( mask ) );
		final ImgLabeling< Integer, IntType > watershedLabeling = new ImgLabeling<>( watershedLabelImg );

		opService.image().watershed(
				watershedLabeling,
				Utils.invertedView( distances ),
				seedsLabelImg,
				false,
				false );

		Utils.applyMask( watershedLabelImg, mask );
		return watershedLabeling;
	}


	public double[] getCorrectedCalibration()
	{
		return correctedCalibration;
	}

	public double getCoverslipPosition()
	{
		return coverslipPosition;
	}

	public Img< IntType > getWatershedLabelImg()
	{
		return watershedLabelImg;
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

		final RandomAccessibleInterval< BitType > dilatedMask = Algorithms.dilate( embryoMask, 2 );

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

		if ( settings.showIntermediateResults ) show( alignedMask, "aligned mask at output resolution", Transforms.origin(), resolution );

		return alignedMask;
	}

	private  < T extends RealType< T > & NativeType< T > >
	AffineTransform3D computeIntensityBasedRollTransform(
			RandomAccessibleInterval< T > image )
	{
		Logger.log( "Computing intensity based roll transform" );

		final AffineTransform3D intensityBasedRollTransform = computeIntensityBasedRollTransform(
				image,
				settings.projectionXMin,
				settings.projectionXMax,
				settings.projectionBlurSigma,
				registrationCalibration );

		return intensityBasedRollTransform;
	}

	private < T extends RealType< T > & NativeType< T > >
	double getThreshold( RandomAccessibleInterval< T > downscaled )
	{

		double threshold = 0;

		if ( settings.thresholdModality.equals( DrosophilaRegistrationSettings.HUANG_AUTO_THRESHOLD ) )
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

	private ImgLabeling< Integer, IntType > createWatershedSeeds( RandomAccessibleInterval< DoubleType > distance )
	{
		Logger.log( "Seeds for watershed...");

		double globalDistanceThreshold = Math.pow( settings.watershedSeedsGlobalDistanceThreshold / settings.registrationResolution, 2 );
		double localMaximaDistanceThreshold = Math.pow( settings.watershedSeedsLocalMaximaDistanceThreshold / settings.registrationResolution, 2 );
		int localMaximaSearchRadius = (int) ( settings.watershedSeedsLocalMaximaSearchRadius / settings.registrationResolution );

		final RandomAccessibleInterval< BitType >  seeds = Algorithms.createWatershedSeeds(
				distance,
				new HyperSphereShape( localMaximaSearchRadius ),
				globalDistanceThreshold,
				localMaximaDistanceThreshold );

		final ImgLabeling< Integer, IntType > seedsLabelImg = Utils.asImgLabeling( seeds, ConnectedComponents.StructuringElement.FOUR_CONNECTED );

		if ( settings.showIntermediateResults ) show( Utils.asIntImg( seedsLabelImg ), "watershed seeds", null, Utils.as3dDoubleArray( settings.registrationResolution ), false );

		return seedsLabelImg;
	}

	private AffineTransform3D computeFlippingTransform( RandomAccessibleInterval yawAlignedMask, RandomAccessibleInterval yawAlignedIntensities, double calibration )
	{
		final CoordinatesAndValues coordinatesAndValues = Utils.computeAverageIntensitiesAlongAxisWithinMask( yawAlignedIntensities, yawAlignedMask, X, calibration );

		if ( settings.showIntermediateResults ) Plots.plot( coordinatesAndValues.coordinates, coordinatesAndValues.values, "x", "average intensity" );

		double maxLoc = CurveAnalysis.maxLoc( coordinatesAndValues, null );

		AffineTransform3D affineTransform3D = new AffineTransform3D();

		if ( maxLoc > 0 ) affineTransform3D.rotate( Z, toRadians( 180.0D ) );

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
		final RandomAccessibleInterval< T > longAxisProjection = Utils.createAverageProjectionAlongAxis(
				rai,
				X,
				xMin,
				xMax,
				settings.registrationResolution );

		if ( settings.showIntermediateResults ) show( longAxisProjection, "channel2 projection", null, registrationCalibration, false );

		final RandomAccessibleInterval< T > blurred = Utils.createBlurredRai(
				longAxisProjection,
				blurSigma,
				settings.registrationResolution );

		final Point maximum = Algorithms.getMaximumLocation( blurred, Utils.as2dDoubleArray( settings.registrationResolution ));
		final List< RealPoint > realPoints = Utils.asRealPointList( maximum );
		realPoints.add( new RealPoint( new double[]{ 0, 0 } ) );

		if ( settings.showIntermediateResults ) show( blurred, "perpendicular projection - blurred ", realPoints, Utils.as2dDoubleArray( settings.registrationResolution ), false );

		final AffineTransform3D xAxisRollTransform = createXAxisRollTransform( maximum );

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
