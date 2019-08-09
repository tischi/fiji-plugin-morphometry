package de.embl.cba.morphometry.spindle;

import ij.measure.Calibration;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import java.io.File;

public class SpindleMorphometrySettings <T extends RealType<T> & NativeType< T > >
{
	public boolean showIntermediateResults = false;
	public double workingVoxelSize = 6.0;

	public double watershedSeedsGlobalDistanceThreshold = Double.MAX_VALUE;
	public double watershedSeedsLocalMaximaDistanceThreshold = 3 * workingVoxelSize; // at least 3 pixels
	public double thresholdInUnitsOfBackgroundPeakHalfWidth = 5.0;

	public double[] inputCalibration;
	public RandomAccessibleInterval<T> dnaImage;
	public RandomAccessibleInterval<T> tubulinImage;

	public double maxDnaAxisDist;
	public double interestPointsRadius;
	public File outputDirectory;
	public String inputDataSetName;
	public double derivativeDelta;

	public double minimalDnaFragementsVolume = 5; // um^3
	public double maxCentralObjectRegionsDistance = 7; // um
	public double erosionOfDnaMaskInCalibratedUnits = 1.0; // um
	public Calibration imagePlusCalibration;
	public double maxSpindlePoleRefinementDistance;
	public double spindleDerivativeDelta = 1;
	public double dnaThresholdFactor = 1.0;
	public double dnaThresholdResolution = 1.5;
}
