package de.embl.cba.morphometry.microglia;

import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import java.io.File;

public class MicrogliaMorphometrySettings<T extends RealType<T> & NativeType< T > >
{
	public static final String MANUAL_THRESHOLD = "Manual threshold";
	public static final String HUANG_AUTO_THRESHOLD = "Huang auto threshold";

	// all spatial values are in micrometer
	// morphometry length: 420
	// morphometry width: 160

	public OpService opService;

	public boolean showIntermediateResults = false;
	public double workingVoxelSize = 6.0;
	public double outputResolution = 2.0;


	public double watershedSeedsGlobalDistanceThreshold = Double.MAX_VALUE;
	public double watershedSeedsLocalMaximaDistanceThreshold = 3 * workingVoxelSize; // at least 3 pixels

	public String thresholdModality = MANUAL_THRESHOLD;
	public double thresholdInUnitsOfBackgroundPeakHalfWidth = 5.0;
	public double closingRadius = 3.0;


	public double[] inputCalibration;
	public RandomAccessibleInterval<T> image;
	public RandomAccessibleInterval<T> tubulin;

	public double maxPossibleValueInDataSet;
	public double maxShortAxisDist;
	public double interestPointsRadius;
	public File outputDirectory;
	public String inputDataSetName;
	public boolean returnEarly;
	public double minimalObjectSize;
	public boolean splitTouchingObjects = false;
	public double skeletonMaxLength;
	public double minimalObjectCenterDistance;
	public double maximalWatershedLength;
	public double minimalOverlapFraction = 0.05;
	public double minimalSumIntensityRatio = 0.5;
	public double maximalSumIntensityRatio = 1.5;
}
