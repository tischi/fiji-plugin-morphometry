package de.embl.cba.morphometry.spindle;

import ij.ImagePlus;
import ij.measure.Calibration;
import net.imagej.DatasetService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;
import org.ilastik.ilastik4ij.ui.IlastikOptions;
import org.scijava.app.StatusService;
import org.scijava.display.DisplayService;
import org.scijava.log.LogService;
import org.scijava.ui.UIService;

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

	public double maxDnaAxisDist;
	public double interestPointsRadius;
	public File outputDirectory;
	public String inputDataSetName;
	public double derivativeDelta;

	public double minimalDnaFragmentsVolume = 5; // um^3
	public double maxCentralObjectRegionsDistance = 7; // um
	public double cellRadius = 6.0; // um
	public double erosionOfDnaMaskInCalibratedUnits = 1.0; // um
	public Calibration imagePlusCalibration;
	public double maxSpindlePoleRefinementDistance;
	public double spindleDerivativeDelta = 1;
	public double dnaThresholdFactor = 0.5;
	public double dnaThresholdResolution = 1.5;
	public int minimalDynamicRange = 7;
	public String version;
	public long dnaChannelIndex;
	public long tubulinChannelIndex;
	public boolean showOutputImage = false;
	public boolean showMetaphaseClassification = false;
	public String classifier;
	public File classifierFile;
	public File classifierExecutionFile;
	public CellCenterDetectionMethod cellCenterDetectionMethod;
	// SciJava
	public IlastikOptions ilastikOptions;
	public LogService logService;
	public StatusService statusService;
	public DatasetService datasetService;
	public ImagePlus imagePlus;
	public DisplayService displayService;
	public UIService uiService;


	public enum CellCenterDetectionMethod
	{
		None,
		BlurredDnaImage,
		BlurredTubulinImage
	}

	// TODO: this could be removed once the Command's can have enums as choices
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
