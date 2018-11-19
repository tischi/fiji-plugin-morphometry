package de.embl.cba.morphometry.drosophila.shavenbaby;

import net.imglib2.FinalInterval;

public class ShavenBabyRegistrationSettings
{
	public static final String MANUAL_THRESHOLD = "Manual threshold";
	public static final String HUANG_AUTO_THRESHOLD = "Huang auto threshold";
	public static final String CENTROID_SHAPE_BASED_ROLL_TRANSFORM = "Shape - Centroids";
	public static final String INTENSITY_BASED_ROLL_TRANSFORM = "Intensity - Other Channel";
	public static final String PROJECTION_SHAPE_BASED_ROLL_TRANSFORM = "Shape - Projection";

	// all spatial values are in micrometer
	// morphometry length: 420
	// morphometry width: 160

	public static double drosophilaLength = 420;
	public static double drosophilaWidth = 160;

	public int svbChannelIndexOneBased = 2;
	public int otherChannelIndexOneBased = 1;

	public boolean showIntermediateResults = false;

	public double refractiveIndexAxialCalibrationCorrectionFactor = 1.6;
	public double refractiveIndexIntensityCorrectionDecayLength = 170; //170;

	public double registrationResolution = 6.0;
	public double outputResolution = 0.7;

	public double rollAngleMinDistanceToAxis = 0;
	public double rollAngleMinDistanceToCenter = drosophilaLength / 2.0 * 0.5;
	public double rollAngleMaxDistanceToCenter = drosophilaLength / 2.0 - 10.0;

	public double watershedSeedsGlobalDistanceThreshold = drosophilaWidth / 3.0;
	public double watershedSeedsLocalMaximaDistanceThreshold = 0.0; //drosophilaWidth / 6.0; // at least 3 pixels

	public String thresholdModality = MANUAL_THRESHOLD;
	public double thresholdInUnitsOfBackgroundPeakHalfWidth = 5.0;
	public double closingRadius = 10;

	public static double[] outputImageSize =
			new double[]{ drosophilaLength * 1.3,
					drosophilaWidth * 1.3,
					drosophilaWidth * 1.3};

	public double minimalObjectSize = drosophilaWidth * drosophilaWidth * drosophilaWidth;

	public double ch2ProjectionXMin = +20.0;
	public double ch2ProjectionXMax = +80.0;
	public double ch2ProjectionBlurSigma = 20.0;
	public double finalProjectionMinDistanceToCenter = 60;
	public String rollAngleComputationMethod = CENTROID_SHAPE_BASED_ROLL_TRANSFORM;
	public double watershedSeedsLocalMaximaSearchRadius = 2 * registrationResolution;

	public FinalInterval getOutputImageInterval()
	{
		final long[] min = new long[ 3 ];
		final long[] max = new long[ 3 ];

		for ( int d = 0; d < 3; ++d )
		{
			min[ d ] = - ( long ) ( outputImageSize[ d ] / 2.0 / outputResolution );
			max[ d ] = -1 * min[ d ];
		}

		return new FinalInterval( min, max );
	}
}
