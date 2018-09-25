package de.embl.cba.morphometry.drosophila.shavenbaby;

import net.imglib2.FinalInterval;

import static de.embl.cba.morphometry.Constants.X;
import static de.embl.cba.morphometry.Constants.Y;
import static de.embl.cba.morphometry.Constants.Z;

public class ShavenBabyRegistrationSettings
{
	public static final String MANUAL_THRESHOLD = "Manual threshold";
	public static final String HUANG_AUTO_THRESHOLD = "Huang auto threshold";
	public static final String CENTROID_SHAPE = "Shape - Centroids";
	public static final String AMNIOSEROSA = "Intensity - Amnioserosa";
	public static final String PROJECTION_SHAPE = "Shape - Projection";

	// all spatial values are in micrometer
	// morphometry length: 420
	// morphometry width: 160

	public static double drosophilaLength = 420;
	public static double drosophilaWidth = 160;

	public int shavenbabyChannelIndexOneBased = 1;
	public int amnioserosaChannelIndexOneBased = 2;

	public boolean showIntermediateResults = false;

	public double refractiveIndexScalingCorrectionFactor = 1.6;
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

	public double amaProjectionXMin = +20.0;
	public double amaProjectionXMax = +80.0;
	public double amaProjectionBlurSigma = 20.0;
	public double finalProjectionMinDistanceToCenter = 60;
	public String rollAngleComputationMethod = CENTROID_SHAPE;
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
