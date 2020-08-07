package de.embl.cba.morphometry.spindle;

import ij.measure.Calibration;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

import java.io.File;

import static de.embl.cba.morphometry.spindle.SpindleMorphometrySettings.CellCenterDetectionMethod.BlurredDnaImage;
import static de.embl.cba.morphometry.spindle.SpindleMorphometrySettings.CellCenterDetectionMethod.BlurredTubulinImage;

public class SpindleMorphometrySettings <T extends RealType<T> & NativeType< T > >
{
	public boolean showIntermediateResults = false;
	public double workingVoxelSize = 0.25;

	public double watershedSeedsGlobalDistanceThreshold = Double.MAX_VALUE;
	public double watershedSeedsLocalMaximaDistanceThreshold = 3 * workingVoxelSize; // at least 3 pixels
	public double thresholdInUnitsOfBackgroundPeakHalfWidth = 5.0;

	public double[] inputCalibration;

	public double interestPointsRadius;
	public File outputDirectory;
	public String inputDataSetName;
	public double derivativeDelta;

	/**
	 * Spatial
	 */
	// TODO: make all micrometer everything relative to something
	public double maxDnaLateralRadius = 6; // um
	public double minimalDnaFragmentsVolume = 5; // um^3
	public double maxCentralObjectRegionsDistance = 7; // um
	public double cellRadius = 6.0; // um
	public double erosionOfDnaMaskInCalibratedUnits = 1.0; // um
	public double maxSpindlePoleRefinementDistance = 3.0;
	public double spindleDerivativeDelta = 1;

	/**
	 * Intensity
	 */
	public double dnaThresholdFactor = 0.5;
	public double dnaThresholdResolution = 1.5;
	public int minimalDynamicRange = 7;

	/**
	 * Other
	 */
	public String version;
	public Calibration imagePlusCalibration;
	public long dnaChannelIndex;
	public long tubulinChannelIndex;
	public boolean showOutputImage = false;
	public boolean showMetaphaseClassification = false;
	public boolean useCATS = false;
	public File classifier;
	public CellCenterDetectionMethod cellCenterDetectionMethod;

	public enum CellCenterDetectionMethod
	{
		None,
		BlurredDnaImage,
		BlurredTubulinImage
	}

	public static final String CCDM_NONE = "None";
	public static final String CCDM_DNA = "BlurredDnaImage";
	public static final String CCDM_TUBULIN = "BlurredTubulinImage";

	public String toString()
	{
		String settings = new String();

		settings += "\n";
		settings += "## Spindle Morphometry Settings\n";
		settings += "workingVoxelSize: " + workingVoxelSize + "\n";
		settings += "dnaThresholdResolution: " + dnaThresholdResolution + "\n";
		settings += "dnaThresholdFactor: " + dnaThresholdFactor + "\n";
		settings += "spindleDerivativeDelta: " + spindleDerivativeDelta + "\n";
		settings += "minimalDynamicRange: " + minimalDynamicRange + "\n";
		settings += "minimalDnaFragmentsVolume: " + minimalDnaFragmentsVolume + "\n";
		settings += "maxCentralObjectRegionsDistance: " + maxCentralObjectRegionsDistance + "\n";
		settings += "erosionOfDnaMaskInCalibratedUnits: " + erosionOfDnaMaskInCalibratedUnits + "\n";

		return settings;
	}

}
