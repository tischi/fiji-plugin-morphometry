package de.embl.cba.morphometry.drosophila.registration;

import net.imglib2.FinalInterval;

public class DrosophilaRegistrationSettings
{
	public static final String MANUAL_THRESHOLD = "Manual threshold";
	public static final String HUANG_AUTO_THRESHOLD = "Huang auto threshold";
	public static final String CENTROID_SHAPE_BASED_ROLL_TRANSFORM = "Shape - Centroids";
	public static final String PROJECTION_SHAPE_BASED_ROLL_TRANSFORM = "Shape - Projection";
	public static final String MOMENTS = "Moments";
	public static final String INTENSITY =  "Intensity";

	// all spatial values are in micrometer

	public static double drosophilaLength = 420;
	public static double drosophilaWidth = 160;

	public int alignmentChannelIndexOneBased = 2;
	public int secondaryChannelIndexOneBased = 2;

	public boolean showIntermediateResults = false;

	public double refractiveIndexAxialCalibrationCorrectionFactor = 1.6;
	public double refractiveIndexIntensityCorrectionDecayLength = 170; //170;

	public double registrationResolution = 6.0;
	public double outputResolution = 0.7;

	public double rollAngleMinDistanceToAxis = 0;
	public double rollAngleMinDistanceToCenter = drosophilaLength / 2.0 * 0.5;
	public double rollAngleMaxDistanceToCenter = drosophilaLength / 2.0 - 10.0;

	public double watershedSeedsGlobalDistanceThreshold = drosophilaWidth / 3.0;
	public double watershedSeedsLocalMaximaDistanceThreshold = 0.0;

	public String thresholdModality = MANUAL_THRESHOLD;
	public double thresholdInUnitsOfBackgroundPeakHalfWidth = 5.0;

	public static double[] outputImageSize =
			new double[]{ drosophilaLength * 1.5,
					drosophilaWidth * 1.5,
					drosophilaWidth * 1.5};

	public double minimalObjectSize = drosophilaWidth * drosophilaWidth * drosophilaWidth;

	public double projectionXMin = +20.0;
	public double projectionXMax = +80.0;
	public double projectionBlurSigma = 20.0;
	public double finalProjectionMinDistanceToCenter = 60;
	public String rollAngleComputationMethod = INTENSITY;
	public double watershedSeedsLocalMaximaSearchRadius = 2 * registrationResolution;
	public String yawTransformComputationMethod;
	public double centralRegionDistance = drosophilaWidth * 0.5;
	public boolean onlyComputeEllipsoidParameters = false;

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
