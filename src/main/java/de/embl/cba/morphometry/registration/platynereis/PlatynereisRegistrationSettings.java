package de.embl.cba.morphometry.registration.platynereis;

import net.imglib2.FinalInterval;

public class PlatynereisRegistrationSettings
{

	// all spatial values are in micrometer

	public static double platyWidth = 140;
	public static double platyHeight = 80;
	public static double platyLength = 280;

	public boolean showIntermediateResults = false;

	public double registrationResolution = 8.0;
	public double outputResolution = 1.0;

	public double thresholdInUnitsOfBackgroundPeakHalfWidth = 5.0;

	public static double[] outputImageSize =
			new double[]{
					platyLength * 1.5,
					platyWidth * 2.0,
					platyHeight * 2.0};

	public double minimalObjectSize = 0.5 * platyWidth * platyHeight * platyLength;

	public double projectionXMin = -8.0;
	public double projectionXMax = +8.0;
	public double projectionBlurSigma = 8.0;
	public boolean invertImage = false;

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
