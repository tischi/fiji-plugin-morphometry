package de.embl.cba.morphometry.spindle;

import ij.measure.Calibration;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import java.io.File;

public class SpindleMorphometrySettings <T extends RealType<T> & NativeType< T > >
{
	public static final String MANUAL_THRESHOLD = "Manual threshold";
	public static final String HUANG_AUTO_THRESHOLD = "Huang auto threshold";

	// all spatial values are in micrometer
	// morphometry length: 420
	// morphometry width: 160

	public static double drosophilaLength = 420;
	public static double drosophilaWidth = 160;

	public boolean showIntermediateResults = false;
	public double workingVoxelSize = 6.0;
	public double outputResolution = 2.0;

	public double rollAngleMinDistanceToAxis = 0;
	public double rollAngleMinDistanceToCenter = drosophilaLength / 2.0 * 0.5;
	public double rollAngleMaxDistanceToCenter = drosophilaLength / 2.0 - 10.0;

	public double watershedSeedsGlobalDistanceThreshold = Double.MAX_VALUE;
	public double watershedSeedsLocalMaximaDistanceThreshold = 3 * workingVoxelSize; // at least 3 pixels

	public String thresholdModality = MANUAL_THRESHOLD;
	public double thresholdInUnitsOfBackgroundPeakHalfWidth = 5.0;
	public double closingRadius = 1.0;

	public double outputImageSizeX = 500;
	public double outputImageSizeY = 250;
	public double outputImageSizeZ = 250;


	public double[] inputCalibration;
	public RandomAccessibleInterval<T> dnaImage;
	public RandomAccessibleInterval<T> tubulinImage;

	public double maxPossibleValueInDataSet;
	public double maxShortAxisDist;
	public double interestPointsRadius;
	public File outputDirectory;
	public String inputDataSetName;
	public double derivativeDelta;

	public double minimalMetaphasePlateVolumeInCalibratedUnits = 5; // um^3
	public double centralObjectRegionToleranceInCalibratedUnits = 3; // um
	public double erosionOfDnaMaskInCalibratedUnits = 1.0; // um
	public Calibration imagePlusCalibration;
}
