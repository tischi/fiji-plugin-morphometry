package de.embl.cba.morphometry.drosophila.shavenbaby;

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
	public double refractiveIndexIntensityCorrectionDecayLength = 250; //170;

	public double registrationResolution = 6.0;
	public double outputResolution = 2.0;

	public double rollAngleMinDistanceToAxis = 0;
	public double rollAngleMinDistanceToCenter = drosophilaLength / 2.0 * 0.5;
	public double rollAngleMaxDistanceToCenter = drosophilaLength / 2.0 - 10.0;

	public double watershedSeedsGlobalDistanceThreshold = drosophilaWidth / 3.0;
	public double watershedSeedsLocalMaximaDistanceThreshold = 3 * registrationResolution; // at least 3 pixels

	public String thresholdModality = MANUAL_THRESHOLD;
	public double thresholdInUnitsOfBackgroundPeakHalfWidth = 5.0;
	public double closingRadius = 10;

	public double outputImageSizeX = 500;
	public double outputImageSizeY = 250;
	public double outputImageSizeZ = 250;

	public double minimalObjectSize = drosophilaWidth * drosophilaWidth * drosophilaWidth;

	public double amaProjectionXMin = +20.0;
	public double amaProjectionXMax = +80.0;
	public double amaProjectionBlurSigma = 20.0;
	public double finalProjectionMinDistanceToCenter = 60;
	public String rollAngleComputationMethod = CENTROID_SHAPE;
}
