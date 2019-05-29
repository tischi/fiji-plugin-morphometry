package de.embl.cba.morphometry.commands;

import de.embl.cba.morphometry.Logger;
import de.embl.cba.morphometry.Utils;
import de.embl.cba.morphometry.measurements.Measurements;
import de.embl.cba.morphometry.spindle.SpindleMorphometry;
import de.embl.cba.morphometry.spindle.SpindleMorphometrySettings;
import de.embl.cba.tables.Tables;
import ij.IJ;
import ij.ImagePlus;
import net.imagej.ops.OpService;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;
import org.scijava.command.Command;
import org.scijava.plugin.Parameter;
import org.scijava.plugin.Plugin;

import javax.swing.*;
import java.io.File;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Map;


@Plugin(type = Command.class, menuPath = "Plugins>Measure>Spindle Morphometry" )
public class SpindleMorphometryCommand< R extends RealType< R > > implements Command
{
	@Parameter
	public OpService opService;

	@Parameter ( label = "Input Image File" )
	public File inputImageFile;

	// TODO: how to make this more clear and easy?
	// (it is to remove to part of the path to only store a relative directory )
	@Parameter ( label = "Input Image Files Parent Directory", style = "directory" )
	public File inputImageFilesParentDirectory;

	@Parameter ( label = "Output Directory", style = "directory" )
	public File outputDirectory;

	@Parameter ( label = "Voxel Size for Analysis" )
	public double voxelSpacingDuringAnalysis = 0.25;

	@Parameter ( label = "DNA Channel [one-based index]" )
	public long dnaChannelIndexOneBased = 2;

	@Parameter ( label = "Spindle Channel [one-based index]" )
	public long spindleChannelIndexOneBased = 1;

	@Parameter ( label = "Show Intermediate Results" )
	public boolean showIntermediateResults = false;

	private String imageName;
	private SpindleMorphometrySettings settings = new SpindleMorphometrySettings();
	private HashMap< Integer, Map< String, Object > > objectMeasurements;

	public void run()
	{
		setSettingsFromUI();
		processFile( inputImageFile );
	}

	private void setSettingsFromUI()
	{
		settings.showIntermediateResults = showIntermediateResults;
		settings.workingVoxelSize = voxelSpacingDuringAnalysis;
		settings.maxDnaAxisDist = 6;
		settings.derivativeDelta = 3.0; // TODO: how to set this?
		settings.spindleDerivativeDelta = 1.0;
		settings.thresholdInUnitsOfBackgroundPeakHalfWidth = 5.0;
		settings.watershedSeedsLocalMaximaDistanceThreshold = 1.0;
		settings.watershedSeedsGlobalDistanceThreshold = 2.0;
		settings.maxSpindlePoleRefinementDistance = 3.0;
		settings.interestPointsRadius = 0.5;
		settings.outputDirectory = outputDirectory;
	}

	private void processFile( File file )
	{
		imageName = inputImageFile.getName().replace( ".tif", "" );

		logStart();

		final ImagePlus imagePlus = IJ.openImage( file.toString() );
		setSettingsFromImagePlus( imagePlus );

		final RandomAccessibleInterval< R > rai = ImageJFunctions.wrapReal( imagePlus );

		final RandomAccessibleInterval< R > dapi =
				Views.hyperSlice( rai, 2, dnaChannelIndexOneBased - 1 );

		final RandomAccessibleInterval< R > tubulin =
				Views.hyperSlice( rai, 2, spindleChannelIndexOneBased - 1 );

		settings.dnaImage = dapi;
		settings.tubulinImage = tubulin;

		SpindleMorphometry morphometry = new SpindleMorphometry( settings, opService );
		morphometry.run();

		objectMeasurements = morphometry.getObjectMeasurements();

		saveOutputImageAndAddImagePathsToMeasurements( morphometry.getOutputImage() );

		saveMeasurements( morphometry );

		logEnd();

	}

	private void logEnd()
	{
		Logger.log( "Done!" );
	}

	private void logStart()
	{
		Logger.log( " " );
		Logger.log( "# Spindle Morphometry" );
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

		final File outputImageFile =
				new File( getOutputDirectory() + imageName + "-out.tif" );

		outputImageFile.getParentFile().mkdirs();

		addImagePathToMeasurements( parentPath, inputImageFile, objectMeasurements, "Path_InputImage" );

		addImagePathToMeasurements( parentPath, outputImageFile, objectMeasurements, "Path_OutputImage" );

		Logger.log( "Saving: " + outputImageFile );
		IJ.saveAsTiff( imagePlus, outputImageFile.toString() );
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
