package de.embl.cba.morphometry.spindle;

import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.measurements.Measurements;

import java.util.HashMap;
import java.util.Map;

public class SpindleMeasurements
{
	public static final String DNA_AXIAL_EXTEND = "DNA_Width";
	public static final String DNA_LATERAL_EXTEND = "DNA_Length";
	public static final String DNA_VOLUME = "DNA_Volume";
	public static final String SPINDLE_VOLUME = "Spindle_Volume";
	public static final String DNA_RELATIVE_CENTRAL_INTENSITY = "DNA_Normalised_Central_Intensity";
	public static final String DNA_SPINDLE_CENTER_DISTANCE = "Dna_Center_To_Spindle_Center_Distance";
	public static final String SPINDLE_AXIS_TO_COVERSLIP_PLANE_ANGLE_DEGREES = "Spindle_Axis_To_Coverslip_Plane_Angle_Degrees";
	public static final String LENGTH_UNIT = "um";
	public static final String VOLUME_UNIT = "um3";
	public static final int ALIGNED_DNA_AXIS = 2;
	public static final String COMMENT = "Comment";
	public static final String SPINDLE_POLE_REFINEMENT_DISTANCE = "Spindle_Pole_Refinement_Distance";
	public static final String SPINDLE_AXIAL_EXTEND = "Spindle_Length";
	public static final String SPINDLE_LATERAL_EXTEND = "Spindle_Width";
	public static final String ANALYSIS_INTERRUPTED_LOW_DYNAMIC_DNA =
			"Analysis interrupted: Too low dynamic range in DNA image";
	public static final String ANALYSIS_INTERRUPTED_LOW_DYNAMIC_TUBULIN =
			"Analysis interrupted: Too low dynamic range in tubulin image";
	public static final String ANALYSIS_FINISHED = "Analysis finished.";
	public static final String SEP = "_";
	public static final String SPINDLE_COV = "Spindle_CoV";

	public Double dnaLateralExtend = Double.NaN;
	public Double dnaAxialExtend = Double.NaN;
	public Double dnaVolumeCalibrated = Double.NaN;
	public Double dnaRelativeCentralIntensity = Double.NaN;
	public Double spindleAxialExtend = Double.NaN;
	public Double spindlePoleARefinementDistance = Double.NaN;
	public Double spindlePoleBRefinementDistance = Double.NaN;
	public Double spindleThreshold = Double.NaN;
	public Double spindleVolume = Double.NaN;
	public Double spindleWidthMin = Double.NaN;
	public Double spindleWidthMax = Double.NaN;
	public Double dnaCenterToSpindleCenterDistance = Double.NaN;
	public Double angleSpindleAxisToCoverslipPlaneInDegrees = Double.NaN;
	public Double dnaVolumeThreshold = Double.NaN;
	public String log = ANALYSIS_FINISHED;
	public Double dnaInitialThreshold = Double.NaN;
	public String version;
	public Double spindleCoefficientOfVariation;

	private HashMap< Integer, Map< String, Object > > objectMeasurements;

	public SpindleMeasurements( HashMap< Integer, Map< String, Object > > objectMeasurements )
	{
		this.objectMeasurements = objectMeasurements;
	}


	public void setObjectMeasurements( )
	{
		addMeasurement( "Version", version );

		addMeasurement( "DNA_Initial_Threshold", dnaInitialThreshold );

		addMeasurement( "DNA_Volume_Threshold", dnaVolumeThreshold );

		addMeasurement(
				DNA_AXIAL_EXTEND + SEP + LENGTH_UNIT,
				dnaAxialExtend );

		addMeasurement(
				SpindleMeasurements.DNA_LATERAL_EXTEND + SEP + SpindleMeasurements.LENGTH_UNIT,
				dnaLateralExtend );

		addMeasurement(
				SpindleMeasurements.DNA_VOLUME + SEP + SpindleMeasurements.VOLUME_UNIT,
				dnaVolumeCalibrated );

		addMeasurement(
				SpindleMeasurements.DNA_RELATIVE_CENTRAL_INTENSITY,
				dnaRelativeCentralIntensity );

		addMeasurement(
				SPINDLE_POLE_REFINEMENT_DISTANCE + SEP
						+ "PoleA" + SEP + SpindleMeasurements.LENGTH_UNIT,
				spindlePoleARefinementDistance );

		addMeasurement(
				SPINDLE_POLE_REFINEMENT_DISTANCE + SEP
						+ "PoleB" + SEP + SpindleMeasurements.LENGTH_UNIT,
				spindlePoleBRefinementDistance );

		addMeasurement( SPINDLE_AXIAL_EXTEND + SEP + SpindleMeasurements.LENGTH_UNIT,
				spindleAxialExtend );

		addMeasurement( "Spindle_Intensity_Threshold",  spindleThreshold );

		addMeasurement( getSpindleVolumeKey(), spindleVolume );

		addMeasurement( SpindleMeasurements.getSpindleWidthMinKey(), spindleWidthMin );

		addMeasurement( SpindleMeasurements.getSpindleWidthMaxKey(), spindleWidthMax );

		addMeasurement(
				SpindleMeasurements.DNA_SPINDLE_CENTER_DISTANCE + SEP + SpindleMeasurements.LENGTH_UNIT,
				dnaCenterToSpindleCenterDistance );

		addMeasurement(
				SpindleMeasurements.SPINDLE_AXIS_TO_COVERSLIP_PLANE_ANGLE_DEGREES,
				angleSpindleAxisToCoverslipPlaneInDegrees );

		addMeasurement(
				SPINDLE_COV,
				spindleCoefficientOfVariation );

		Measurements.addMeasurement(
				objectMeasurements,
				0,
				COMMENT,
				log );
	}

	private void addMeasurement( String name, Object value )
	{
		Logger.log( name + ": " + value  );

		Measurements.addMeasurement(
				objectMeasurements,
				0,
				name,
				value );
	}

	public static String getSpindleVolumeKey()
	{
		return SpindleMeasurements.SPINDLE_VOLUME + SEP + SpindleMeasurements.VOLUME_UNIT;
	}

	public static String getSpindleWidthMaxKey()
	{
		return "Spindle_Width_Max" + SEP + LENGTH_UNIT;
	}

	public static String getSpindleWidthMinKey()
	{
		return "Spindle_Width_Min" + SEP + LENGTH_UNIT;
	}

	@Deprecated
	public static String getSpindleWidthKey()
	{
		return SPINDLE_LATERAL_EXTEND + SEP + LENGTH_UNIT;
	}

	public static String getDnaLateralExtendKey()
	{
		return DNA_LATERAL_EXTEND
				+ SEP + LENGTH_UNIT;
	}
}
