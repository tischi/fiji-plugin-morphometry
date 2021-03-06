package de.embl.cba.morphometry.commands;

import de.embl.cba.morphometry.ImageScience;
import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.measurements.Measurements;
import de.embl.cba.morphometry.spindle.SpindleMeasurements;
import de.embl.cba.morphometry.spindle.SpindleMorphometry;
import de.embl.cba.morphometry.spindle.SpindleMorphometrySettings;
import de.embl.cba.tables.Tables;
import ij.CompositeImage;
import ij.IJ;
import ij.ImagePlus;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import org.scijava.Context;
import org.scijava.ItemIO;
import org.scijava.ItemVisibility;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import javax.swing.*;
import java.io.File;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;

import static de.embl.cba.morphometry.spindle.SpindleMorphometrySettings.*;


@Plugin(type = Command.class, menuPath = "Plugins>Spindle3D>Spindle3D..." )
public class Spindle3DCommand< R extends RealType< R > > implements Command
{

	public SpindleMorphometrySettings settings = new SpindleMorphometrySettings();

	@Parameter
	private Context context;

	@Parameter
	public OpService opService;

	@Parameter ( label = "Input Image File" )
	public File inputImageFile;

	//	@Parameter ( label = "Input Image Files Parent Directory", style = "directory" )
	public File inputImageFilesParentDirectory = new File("/" ); // to store a relative directory

	@Parameter ( label = "Output Directory", style = "directory" )
	public File outputDirectory;

	// @Parameter ( label = "Voxel Size for Analysis" )
	public double voxelSpacingDuringAnalysis = settings.workingVoxelSize;

	//	@Parameter ( label = "DNA threshold factor" )
	public double dnaThresholdFactor = settings.dnaThresholdFactor;

	// @Parameter ( label = "Minimum Dynamic Range [segmentation threshold gray value]" )
	public int minimalDynamicRange = settings.minimalDynamicRange;

	@Parameter ( label = "DNA Channel [one-based index]" )
	public long dnaChannelIndexOneBased = 2;

	@Parameter ( label = "Spindle Channel [one-based index]" )
	public long spindleChannelIndexOneBased = 1;

	//	@Parameter ( label = "Initial Cell Center Detection Method", choices = { CCDM_NONE, CCDM_DNA, CCDM_TUBULIN } )
	public String cellCenterDetectionMethodChoice = CCDM_NONE;

	//	@Parameter ( label = "Use CATS for Metaphase Detection" )
	public boolean useCATS = false;

	//	@Parameter ( label = "CATS Classifier" )
	public File classifier;

	@Parameter ( label = "Show Intermediate Results" )
	public boolean showIntermediateResults = false;

	@Parameter( visibility = ItemVisibility.MESSAGE )
	private String version = "Spindle Morphometry Version: " + Spindle3DVersion.VERSION;

	@Parameter( type = ItemIO.OUTPUT )
	private double spindleVolume;

	public boolean saveResults = true;

	private String imageName;
	private HashMap< Integer, Map< String, Object > > objectMeasurements;

	public void run()
	{
		if ( ! ImageScience.isAvailable() ) return;
		setSettingsFromUI();
		processFile( inputImageFile );
	}

	private void setSettingsFromUI()
	{
		settings.showIntermediateResults = showIntermediateResults;
		settings.workingVoxelSize = voxelSpacingDuringAnalysis;
		settings.maxDnaLateralRadius = 6;
		settings.derivativeDelta = 3.0; // TODO: how to set this?
		settings.spindleDerivativeDelta = 1.0;
		settings.thresholdInUnitsOfBackgroundPeakHalfWidth = 5.0;
		settings.watershedSeedsLocalMaximaDistanceThreshold = 1.0;
		settings.watershedSeedsGlobalDistanceThreshold = 2.0;
		settings.maxSpindlePoleRefinementDistance = 3.0;
		settings.interestPointsRadius = 0.5;
		settings.outputDirectory = outputDirectory;
		settings.dnaThresholdFactor = dnaThresholdFactor;
		settings.minimalDynamicRange = minimalDynamicRange;
		settings.version = version;
		settings.useCATS = useCATS;
		settings.classifier  = classifier;
		settings.cellCenterDetectionMethod = SpindleMorphometrySettings.CellCenterDetectionMethod.valueOf( cellCenterDetectionMethodChoice );

		Logger.log( settings.toString() );
	}

	public HashMap< Integer, Map< String, Object > > getObjectMeasurements()
	{
		return objectMeasurements;
	}

	private void processFile( File file )
	{
		setImageName();

		logStart();

		final ImagePlus imagePlus = Utils.openWithBioFormats( file.toString() );
		setSettingsFromImagePlus( imagePlus );

		final RandomAccessibleInterval< R > raiXYCZ = ImageJFunctions.wrapReal( imagePlus );

		settings.dnaChannelIndex = dnaChannelIndexOneBased - 1;
		settings.tubulinChannelIndex = spindleChannelIndexOneBased - 1;

		//final OpService service = context.service( OpService.class );
		SpindleMorphometry morphometry = new SpindleMorphometry( settings, opService );
		final String log = morphometry.run( raiXYCZ );
		Logger.log( log );

		final SpindleMeasurements measurements =
				morphometry.getMeasurements();

		spindleVolume = measurements.spindleVolume;

		objectMeasurements = morphometry.getObjectMeasurements();

		addImagePathToMeasurements(
				inputImageFilesParentDirectory.toPath(),
				inputImageFile,
				objectMeasurements,
				"Path_InputImage" );

		if ( saveResults ) new File( getOutputDirectory() ).mkdirs();

		if ( log.equals( SpindleMeasurements.ANALYSIS_FINISHED ))
		{
			if ( settings.showOutputImage == true || saveResults )
			{
				final CompositeImage outputImage = morphometry.createOutputImage();

				if ( settings.showOutputImage == true )
					outputImage.show();

				if ( saveResults )
					saveOutputImageAndAddImagePathsToMeasurements( outputImage );
			}
		}

		if ( saveResults ) saveMeasurements( morphometry );

		logEnd();
	}

	private void setImageName()
	{
		imageName = inputImageFile.getName().replace( ".tif", "" );
		imageName = inputImageFile.getName().replace( ".ome", "" );
		imageName = inputImageFile.getName().replace( ".zip", "" );
	}

	private void logEnd()
	{
		Logger.log( "Done!" );
	}

	private void logStart()
	{
		Logger.log( "## Spindle Morphometry Measurement" );
		Logger.log( "Processing file " + imageName );
	}

	private void saveMeasurements( SpindleMorphometry morphometry )
	{
		final JTable jTable = Measurements.asTable( objectMeasurements );

		final File tableOutputFile = new File( getOutputDirectory() + "measurements.txt" );

		Logger.log( "Saving " + tableOutputFile );

		Tables.saveTable( jTable, tableOutputFile );
	}

	private String getOutputDirectory()
	{
		return outputDirectory
				+ File.separator
				+ imageName
				+ File.separator;
	}

	private void setSettingsFromImagePlus( ImagePlus imagePlus )
	{
		settings.inputCalibration = Utils.getCalibration( imagePlus );
		settings.imagePlusCalibration = imagePlus.getCalibration();
		settings.inputDataSetName = imagePlus.getTitle();
	}

	private void saveOutputImageAndAddImagePathsToMeasurements( ImagePlus imagePlus )
	{
		final Path parentPath = inputImageFilesParentDirectory.toPath();

		final File outputImageFile = new File( getOutputDirectory() + imageName + "-out.zip" );

		addImagePathToMeasurements( parentPath, outputImageFile, objectMeasurements, "Path_OutputImage" );

		Logger.log( "Saving: " + outputImageFile );
		IJ.saveAs( imagePlus, "ZIP", outputImageFile.toString() );
	}

	private static void addImagePathToMeasurements(
			Path parentPath,
			File inputImageFile,
			HashMap< Integer, Map< String, Object > > objectMeasurements,
			String path_inputImage )
	{
		final Path relativeInputImagePath = parentPath.relativize( inputImageFile.toPath() );

		Measurements.addMeasurement(
				objectMeasurements,
				0,
				path_inputImage,
				relativeInputImagePath );
	}

}
